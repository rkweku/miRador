#!/usr/loca/bin/python3

# miRunner is the miRNA prediction script within the miRador package.
# It functions by first identifying inverted repeats within the genome.
# Then, it maps all small RNAs to the identified inverted repeats to create
# a list of candidate miRNAs.

#### Written by Reza Hammond
#### rkweku@udel.edu 

######## IMPORT ################
import csv,sys,os,re,subprocess,multiprocessing,datetime,time
import statistics
import string
from multiprocessing import Process, Queue, Pool, current_process
from multiprocessing.connection import wait
from itertools import product, repeat

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from Bio.Blast.Applications import NcbiblastnCommandline
from datetime import datetime

########################## EXECUTION VARIABLES ###############################
# Name of the genome file
genomeFilename = 'genome/Maize/Zea_mays.AGPv4.dna.toplevel.fa'
#genomeFilename = 'genome/Arabidopsis/all.fa'
#genomeFilename = 'genome/fake/fakeGenome_c.fa'
# List of library file names
libFilenamesList = ['libs/MAIZE_sbsCSHL_sRNA/2421_chopped.txt', 'libs/MAIZE_sbsCSHL_sRNA/2432_chopped.txt', 'libs/MAIZE_sbsCSHL_sRNA/Mz_dcl1_chopped.txt', 'libs/MAIZE_sbsCSHL_sRNA/1682_chopped.txt']
#libFilenamesList = ['libs/fake/fakeLibs_c.txt']#, 'libs/AT_pub2_sRNA/100_chopped.txt']
#libFilenamesList = ['libs/AT_pub2_sRNA/99_chopped.txt',  'libs/AT_pub2_sRNA/100_chopped.txt']#, 'libs/AT_pub2_sRNA/724_chopped.txt', 'libs/AT_pub2_sRNA/2202_chopped.txt']
#libFilenamesList = ['libs/AT_DAS1_sRNA/2730_chopped.txt', 'libs/AT_DAS1_sRNA/2731_chopped.txt', 'libs/AT_DAS1_sRNA/2732_chopped.txt', 'libs/AT_DAS1_sRNA/2733_chopped.txt', 'libs/AT_DAS1_sRNA/2734_chopped.txt', 'libs/AT_DAS1_sRNA/2735_chopped.txt', 'libs/AT_DAS1_sRNA/2736_chopped.txt', 'libs/AT_DAS1_sRNA/2737_chopped.txt', 'libs/AT_DAS1_sRNA/3205_chopped.txt', 'libs/AT_DAS1_sRNA/3206_chopped.txt', 'libs/AT_DAS1_sRNA/3207_chopped.txt', 'libs/AT_DAS1_sRNA/3208_chopped.txt', 'libs/AT_DAS1_sRNA/3271_chopped.txt', 'libs/AT_DAS1_sRNA/3382_chopped.txt']
# Optional parameter of merged map files for each library. If these are set,
# mapping the libFilenamesList will be bypassed
mapFilenames = []
#mapFilenames = ['libs/AT_pub2_sRNA/99_chopped.map']#
# Penalty for a score
gap = 12            ## Default:12 (from emboss website)
# Score for a match
match = 3           ## Default:3
# Score for a mismatch
mismatch = -4       ## Default:-4
# Minimum score for a candidate miRNA precursor
threshold = 40      ## Default:50 - This is the score threshold for cutoff
# Maximum inverted repeat length that can be identified as a candidate
# Default is 300 based on max suggested by Axtell and Meyers (TPC 2018)
maxRepLen = 300
# Flag to run einverted for this genome file
runEInvertedFlag = 1
# First letter of genus and first 2 letters of species. This is important
# for us to properly annotate miRNAs as novel or not when referencing 
# mirBase. This is the same format as would be found in mirBase.
threeLetterIdentifier = "zma"
# Flag to utilize parallelization
parallel = 1        ## Default: 1 (Yes)
nthreads = 8

bowtieBuildPath = os.path.expanduser('~/tools/bowtie-1.2.2/bowtie-build')
bowtiePath = os.path.expanduser('~/tools/bowtie-1.2.2/bowtie')

# The name of the query miRNA File in FASTA format
queryMirnasFilename = 'output.fa'
# The name of the subject miRNA file in FASTA format
subjectSequencesFilename = 'mirBase/mirBaseMirnas.fa'
# The name of the BLAST database
dbFilename = 'mirBase/mirBaseMirnas.db'

directory = './'


### Steps #####

##############################################################################
class Genome:
    """
    Class grouping of genome functions and data structures
    """

    def __init__(self, filename):
        self.filename = filename
        self.chrDict = {}
        self.IRFastaFilename = 'invertedRepeats/Inverted_Seqs.fa' 
        self.IRAlignmentFilename = 'invertedRepeats/Inverted_Alings.inv'
        
        self.genomeSeqList = self.readFasta(self.filename)

        # Create an empty dictionary in IRDictByChr for as many chromosomes
        # exist
        self.IRDictByChr = [{} for i in range(len(self.chrDict))]

        # Build the bowtie index if it does not exist yet
        self.indexFilename = self.buildBowtieIndex(genomeFilename)

    def readFasta(self, filename):
        """
        Read a genome FASTA file into memory as a list

        Args:
            filename: Name of the file to be read
        Returns:
            Variable wholeFile containing the entire file

        """

        # Create an empty list to hold full genome FASTA file
        fastaList = []
        count = 0
        emptyCount = 0

        # loop through file using csv.reader and store in wholeFile
        f = open(filename, 'r')
        wholeFile = f.read()

        # Loop through the FASTA file and split by sequence name in case
        # there are some sequence names with no sequence on the next line   
        for entry in wholeFile.split('>')[1:]:
            chromoInfo = entry.partition('\n')

            # We will call the sequence identifier of the FASTA file
            # the chrName. Technically, this does not have to be a chromosome
            # if the user's file does not contain full chromosomes, but
            # the variable name will at least indicate that
            chrName = chromoInfo[0].split()[0]
            # Add seqID to chrDict for indexing of chromosomes
            self.chrDict[chrName] = count

            # Strip the newline character from the sequence
            sequence = chromoInfo[2].replace('\n','').strip()

            # Increment the entry counter
            count += 1
            # If the sequence exists, append the sequence ID and the
            # matching sequence to fastaList 
            if(sequence):
                fastaList.append((chrName, sequence))

            else:
                emptyCount += 1
        
        print("Total entries in File: %s | Total empty entries in file: %s" %
            (count, emptyCount))
        print("Total entries in the genome FASTA list: %s" % (len(fastaList)))

        return(fastaList)

    def buildBowtieIndex(self, filename):
        """
        Code to create a bowtie index for the inverited repeats file.

        Args:
            filename: Path and name of the inverted repeats to generate 
                bowtie index for
        Returns:
            Path of bowtie index

        """

        # Set the index filename. Remove any file extension and folders
        # from the filename path to ensure the index file is written
        # to the correct folder that is hardcoded here
        filenameStripped = ''.join(os.path.splitext(
            filename.split('/')[-1])[0])
        indexFilename = "genome/bowtieIndex/%s" % (filenameStripped)

        if(self.checkBowtieNeedsUpdate(indexFilename)):
            print("Building a bowtie index for %s" % (filename))
            with open("genome/bowtieIndex/%s_bowtiebuild.log" %\
                    filenameStripped, 'w') as logFile:
                subprocess.call([bowtieBuildPath, filename, indexFilename],
                    stdout = logFile)

            logFile.close()

        return(indexFilename)

    def checkBowtieNeedsUpdate(self, indexFilename):
        """
        Check if the index filename already exists. Return true if it
        does not and needs to be updated, but return false if doesn't
        need an update.

        Args:
            indexFilename: Name of the index file that would be created
                by bowtie-build for the genome file

        Returns:
            True if no update is needed, false if an update is needed

        """

        # Check if the bowtieIndex file exists in the genome folder.
        # Create it if it does not
        if(not os.path.isdir("genome/bowtieIndex")):
            os.mkdir('genome/bowtieIndex')

        pattern = re.compile("%s(.*?).ebwt" % os.path.basename(indexFilename))

        # Check if a bowtie index file exists, and if it does, return false,
        # otherwise, return true so that bowtie-build can be run
        for filepath in os.listdir('genome/bowtieIndex'):
            if(pattern.match(filepath)):
                return(False)

        return(True) 

    def einverted(self, chrAndSeq):
        """
        Fuunction to run einverted for a sequence
        
        Args:
            chrAndSeq: A tuple where the firstt element of the tuple is 
            the chr (or sequence name) and the second is the sequence
        """

        outputFastaFilenamesList = []
        outputAlignmentFilenamesList = []

        # Open FNULL to suppress the output of einverted becuase we do not
        # really need to know it is running for each proc
        FNULL = open(os.devnull, 'w')

        # Separate the name and sequence from the tuple
        name = chrAndSeq[0]
        seq = chrAndSeq[1]

        # Each sequence needs to be in its own fasta file, so we will
        # first create temp input files for each sequence prior to calling
        # einverted
        tempInput = "%s.fa" % name
        f_out = open(tempInput,'w')
        f_out.write('>%s\n%s\n' % (name,seq))
        f_out.close()

        # Names of temporary output files to store results prior to merging
        outputFastaFilename = "invertedRepeats/%s.fa.temp" % name
        outputAlignmentFilename = "invertedRepeats/%s.alignment.temp" % \
            name

        # Call einverted utilizing this current sequence with the user
        # defined arguments from the config file.
        retcode = subprocess.call(["einverted", "-sequence", tempInput,
            "-gap", str(gap), "-threshold", str(threshold), "-match",
            str(match), "-mismatch", str(mismatch), "-maxrepeat",
            str(maxRepLen), "-outfile", outputAlignmentFilename,
            "-outseq", outputFastaFilename], stdout=FNULL,
            stderr=subprocess.STDOUT)

        # If a return code of anything but 0 is returned, it means there
        # was a problem and it should be investigated. Temp files wiill
        # remain from the run to assist in the debugging process
        if(retcode != 0):
            print("Something wrong happened while running einverted for "\
                "%s" % (tempInput))
            sys.exit()

        ## Cleanup entry specifc FASTA file
        if os.path.exists(tempInput):
            os.remove(tempInput)

        # Close FNULL
        FNULL.close()

        return(outputFastaFilename, outputAlignmentFilename)

    def combineIRTempFiles(self, IRFastaFilenamesList,
            IRAlignmentFilenamesList):
        """
        This function combines the temporary einverted files into one file
        for final analysis. However, if the user has opted to not run 
        einverted due to a previous run alredy existing for this genome,
        this function will bypass the merging steps and only add the IR
        data to the IR dictionary

        Args:
            IRFastaFilenamesList: List of the inverted repeat FASTA files
            IRAlignmentFilenamesList: List of inverted repeat alignment files
            
        """

        # If einverted was run, combine the temp FASTA files
        if(runEInvertedFlag):
            ## Combine inverted repeats sequences FASTA file

            fasta_out = open(self.IRFastaFilename ,'w')

            # Loop through all FASTA files and merge into one file
            for filename in IRFastaFilenamesList:
                with open(filename) as fastaFile:
                    for line in fastaFile:
                        fasta_out.write(line)

            fasta_out.close()

        # Set a counter to process each inverted repeat by line number
        counter = 0
        IRCounter = 0

        # If einverted was run, open the output alignment file to write
        # the results to
        if(runEInvertedFlag):
            ## Combine inverted repeats alignments results
            align_out = open(self.IRAlignmentFilename ,'w')

        # Loop through all alignment files and merge into one file
        for filename in IRAlignmentFilenamesList:
            with open(filename) as alignmentFile:
                # Loop through the alignment files to add them to the
                # merged file and add them to IRDictByChr

                failFlag = 0
                toWriteList = []
                ### *******************
                # I've realized that I can do some filtering here before
                # adding the IRs to the combined file.
                # One such filter is to prevent complex duplex with no
                # secondary stems or large internal loops.
                # I can include a processing step when counter % 5 is 3 
                # to search for more than 5 '|' in a row
                for line in alignmentFile:
                    toWriteList.append(line)

                    # Split the entire line on spaces for parsing.
                    # Remove preceeding whitespace with lstrip first
                    parsedLine = line.lstrip().rstrip().split(' ')

                    # If the current line counter % 5 is 1, there will be
                    # a lot of useful information here. Begin to parse
                    # this data into specific variables
                    # Note that if einverted changes the output format,
                    # these lines here can fail and would need to be
                    # readjusted
                    if(counter % 5 == 1):

                        chrName = parsedLine[0].split(':')[0]
                        score = int(parsedLine[2].split(':')[0])
                        matches, totalBases = map(int,
                            parsedLine[3].split('/'))
                        percMatch = round(float(matches)/totalBases, 3)

                        if(percMatch == 1):
                            gaps = int(parsedLine[6])
                        elif(percMatch < .1):
                            gaps = int(parsedLine[8])
                        else:
                            gaps = int(parsedLine[7])

                    # If the current line counter % 5 is 2, the 5'
                    # repeat start and end coordinates will be contained
                    # within this line
                    elif(counter % 5 == 2):
                        start5 = int(parsedLine[0])
                        hairpin5 = parsedLine[1].upper()
                        end5 = int(parsedLine[2])

                    # The alignment between the two strands is given in
                    # the 3rd line of the alignment.
                    elif(counter % 5 == 3):
                        alignmentIndicators = line.lstrip().rstrip()

                    # If the current line counter % 5 is 4, the 3'
                    # repeat start and end coordinates will be contained
                    # within this line. Also, calculate the loop length
                    elif(counter % 5 == 4):
                        start3 = int(parsedLine[2])
                        hairpin3 = parsedLine[1].upper()
                        end3 = int(parsedLine[0])
                        loop = int(start3) - int(end5) - 1

                        # Get the index of the chromosome to add the
                        # inverted repeat to
                        index = self.chrDict[chrName]

                        # Add the inverted repeat to the appropriate
                        # list within IRDictByChr
                        IRName = "mir%s" % IRCounter
                        self.IRDictByChr[index][IRName] = (start5, end5,
                            start3, end3, loop, 'w', hairpin5,
                            alignmentIndicators, hairpin3)

                        IRCounter += 1
                        IRName = "mir%s" % IRCounter
                        self.IRDictByChr[index][IRName] = (start5, end5,
                            start3, end3, loop, 'c', hairpin5,
                            alignmentIndicators, hairpin3)

                        if(not failFlag and runEInvertedFlag):
                            for entry in toWriteList:
                                align_out.write(entry)
                        toWriteList = []
                        failFlag = 0
                        IRCounter += 1

                    # Increment the counter
                    counter += 1 

        # If einverted was run and the temp files were merged, close
        # the output file and delete the temp files
        if(runEInvertedFlag):
            align_out.close()

            # Delete individual inverted files and fasta files
            print("Combined files '%s and %s'\nDeleting temp files" %\
                (self.IRAlignmentFilename, self.IRFastaFilename))

            # Combine the inverted repeats FASTA and alignmenet filenames
            # lists to delete all of these temp files
            garbage = IRFastaFilenamesList + IRAlignmentFilenamesList
            for file in garbage:
                os.remove(file)

class Library:
    """
    Class grouping of library functions and data structures
    """

    def __init__(self, filename, chrDict):
        self.filename = filename
        self.fastaFilename = "".join(self.filename.split('.')[:-1]) + '.fa'
        self.libDict = self.readTagCount()
        # Create the mapped filename for this library
        self.mapFilename = "%s.map" % ("".join(os.path.splitext(
            self.filename)[0]))

        # Using chrDict from the genome file, create a tuple of multiple
        # dictionaries. Each dictionary in the tuple will represent a 
        # chromosome. Within these dictionaries are keys 'w' and 'c' for
        # each strand of the chromosome. These will then contain empty
        # ditionaries as values which will be populated by positions and
        # the sequences of the reads that map there.
        self.mappedList = [{'w': {}, 'c': {}} for i in range(len(chrDict))]

    def readTagCount(self):
        """
        Read a library in tag count format into memory. Because we are
        performing miRNA identification with these sequences, we know 
        that we will not need all sequences. Additionally, as each 
        sequence is read, we will write this to a FASTA output file for
        bowtie

        Return:
            Dictionary of library read in

        """

        # Start timer for function
        funcStart = time.time()

        # Create an empty dictionary to store the full library
        libDict = {}
        libSize = 0

        # Open a file to write all FASTA sequences to
        g = open(self.fastaFilename, 'w')

        # Open the file and loop through line by line to store the library into
        # a dictionary. Each tag will be a key and the abundance will be the
        # value
        with open(self.filename) as f:
            for line in f:
                # Store the sequence and count into variables, then add them
                # to the dictionary. There should not be any duplicate
                # sequences in this format, however, 
                tag = line.split('\t')[0]
                count = line.split('\t')[1].strip()
                libDict[tag] = [int(count), 0]

                g.write(">s_%s\n%s\n" % (libSize, tag))
                libSize += 1

            print("Total entries in library %s: %s" % (self.filename, libSize))
            g.close()

        # Stop timer for function
        funcEnd = time.time()

        # Calculate the execution time and print it to the user
        execTime = round(funcEnd - funcStart, 2)
        print("Time to read library %s: %s seconds" % (self.filename,
            execTime))

        return(libDict)

    def mapper(self, indexFilename):
        """
        Map small RNAs to the provided index file

        Args:
            indexFilename: Path and name of the index for the genome. 
                This will be a fragment if fragFasta was run
        Returns:
            Filename of mapped data

        """

        # Strip the filename of its folders and create the output map
        # name with that stripped filename in the libs folder
        indexNameStripped = os.path.basename(indexFilename)

        if(not os.path.isdir("libs")):
            os.mkdir("libs")

        outputFilename = '%s.map' % "".join(os.path.splitext(
            self.fastaFilename)[:-1])

        logFilename = "%s_bowtie.log" % os.path.splitext(outputFilename)[:-1] 

        # Run bowtie
        print("Mapping small RNAs to the genome files for %s" %\
            (self.fastaFilename))

        with open(logFilename, 'w') as logFile:
            # Run bowtie with the following options:
            # -a to report all valid alignments as we want multihits
            # -m 20 to suppress all alignments with more than 20 matches
            # to the genome. We expect few multi-matches to the genome
            # --best and --strata ensures only the best alignments are reported
            # and so that less optimum but passable alignments do not appear
            # -v 0 Allow no mismatch
            # --sam-nohead removes the header from the SAM file. This is useful
            # because we have to merge the fragment alignments for parallel
            # runs
            # --no-unal suppresses sequences with no alignemnt. This helps to
            # keep the map file manageable and filter out these sequences
            # earlier for efficiency
            ### Note that the output of bowtie is send to stderr for some 
            ### which is why this log flie goes there
            subprocess.call([bowtiePath, indexFilename, "-f",
                self.fastaFilename, "-a", "-m 20", "--best", "--strata", "-p",
                str(nthreads), "-v 0", "-S", outputFilename, "--sam-nohead",
                "--no-unal"], stderr = logFile)

        logFile.close()

        return(outputFilename, logFilename)

    def createMappedList(self, chrDict):
        """
        Create a dictionary to hold the mapped sRNAs for simple querying

        Args:
            chrDict: This is the genome chromosome dictionary with the
                index position for the chromosome so that mappedList
                and List can be linked by index positions

        Returns:
            Dictionary of the mapped file with chromosomes as a key and 
            subdictionary for each chromosome. This subdictionary contains
            another subdictionary with strand as the key. These
            subdictionaries have a position as a key and then a list as
            values. That list contains all tags that map to that same
            position

        """

        with open(self.mapFilename, 'r') as inFile:
            for line in inFile:
                flag = line.split('\t')[1]
                chrName = line.split('\t')[2]
                position = int(line.split('\t')[3])
                sequence = line.split('\t')[9]

                # If the flag is a star, it means that it did not align
                # and thus must we break from the loop here
                if(chrName == '*'):
                    continue

                # If the 16th bit flag is not set, the strand is w
                elif(not int(flag) & 0b10000):
                    strand = 'w'

                # If the 16th bit of the flag is set, the strand is c. We 
                # thus need to make the sequence the reverse complement to
                # ensure we can find the right sequence later
                elif(int(flag) & 0b10000):
                    strand = 'c'
                    sequence = sequence.translate(str.maketrans(
                        "ACGT","TGCA"))[::-1]

                # Get the index of the chromosome to add to the appropriate
                # index of mappedList
                chrIndex = chrDict[chrName]

                # If the position is not yet in this dictionary, create an
                # entry for it
                if(position not in self.mappedList[chrIndex][strand]):
                    self.mappedList[chrIndex][strand][position] = [sequence]

                # If the position already exists in this dictionary,
                # just append the sequence to the position list
                else:
                    self.mappedList[chrIndex][strand][position].append(sequence)

                # Increment the hits variable for this sequence in libDict
                # to indicate its number of mapped locations found
                self.libDict[sequence][1] += 1

def findOverlappingSeqsAndAbuns(position, subMappedDict, IREnd, libDict,
        IRMappedTags):

    """
    This function will serve a few purposes. First, we need to track the
    total abundance of sRNAs that map to each arm of an inverted repeat.
    It will then add the sequence and abundance of that sequence to a
    variable, provided the sequence length is between 18 and 26 (only
    20-26 are considered candidates, but we need these for variant
    abundance calculations

    Args:
        position: The position of the tags that are present at this location
        subMappedDict: A dictionary of mapped tags and their mapping strands
            at each position. This is the subdictionary of mappedList
            to avoid passing chrName and IRStrand
        IREnd: The last position of the inverted repeat
        libDict: The entire library dictionary to be queried for abundances
        IRMappedTags: Dictionary of tags mapping to the inverted repeat.
            This will be updated with a sequence and hits normalizd abundance
            if the tag the appropriate length and does not extend beyond the
            arm of the inverted repeat

        Returns:
            The total abundance of tags that map to the IR arm, regardless of
            tag length

    """

    totalAbun = 0

    # Pull all tags that map to this position 
    tagList = [mappedTags for mappedTags in \
        subMappedDict[position]]

    # Loop through the list of sequences that map to
    # this precursor
    for sequence in tagList:
        abundance = libDict[sequence][0]
        hits = libDict[sequence][1]
        hitsNormalizedAbundance = abundance/hits

        # We need a quick check to determine if the sequence extends beyond
        # the IR Arm. If it does, we cannot add it as a mapping tag. This
        # should only be true if there are two tags that map to a position
        # where one sequence is longer than the other(s) and extends into
        # the loop or outside the precursor
        if(position + len(sequence) - 1 > IREnd):
            continue

        # Add the HNA to the total abundance
        totalAbun += hitsNormalizedAbundance

        # If the tag is beween 18 and 26 nt in length, add
        # it its sequence and HNA to IRMappedTags.
        # Note that whil we only want 20 and 24, 1 nt variants
        # (tailed or truncated on both 5' and 3' ends) call for 18
        # and 26 nt sRNAs to be recorded in the event that a 20nt or
        # 24 nt sRNA are candidates. We will eliminate tags from being
        # miRNA candidates early in the filtration step
        if(len(sequence) >= 18 and len(sequence) <= 26):
            IRMappedTags[position].append((
                sequence, hitsNormalizedAbundance))

    return(totalAbun)

def mapSRNAsToIRs(IRDict, mappedDict, libDict):
    """
    Map small RNAs to the inverted repeats. This will first read in
    the map file (merged mapfile if run in parallel) and identify any
    small RNAs which overlap with inverted repeats. Any small RNAs that
    do map will be added to the IRDict with the coordinates of overlap.
    Do not count small RNAs that map to the loop of the precursor

    Args:
        IRDict: Dictionary of the inverted repeats in one chromosome
        mappedDict: Dictionary of positions and sequences that map to those
            positions in one chromosome
        libDict: Dictionary of all library tags and their abundances

    """

    mirnasCount = 0
    listOfMirnas = []

    # Initialize a dictionary to hold the precursors
    precursorsDict = {}

    # Loop through the inverted repeats to find all tags within our
    # length cutoffs that map to them
    for IRName, invertedRepeat in IRDict.items():
        # Get the strand of the inverted repeat and various other
        # information about the inverted repeat
        start5 = invertedRepeat[0]
        end5 = invertedRepeat[1]
        start3 = invertedRepeat[2]
        end3 = invertedRepeat[3]
        loopStart = end5 + 1
        loopEnd = start3 - 1
        IRStrand = invertedRepeat[5]

        # Initialize dictionaries to store postions of sRNAs that
        # map to inverted repeats
        IRMappedTags5 = {}
        IRMappedTags3 = {}

        # Initialize a list to store the positions that tags map to
        # within the loop so that we can get their abundance
        mappedLoopPositions = []

        # Initialize variables to store the total abundance of tags
        # that map to the precursor
        totalAbun5 = 0
        totalAbun3 = 0
        loopAbun = 0

        # Loop through all postiions that map to the inverted 
        # repeats and create an entry in a dictionary on the
        # 5' and 3' variables, respectively. Initialize values
        # (which will be abundance) to an empty list. We use a list
        # as it is possible two tags map to the same location, so we
        # need to leave the option open that this will happen
        for position in range(start5, end3):
            # We are looping through this dictionary in key order,
            # so if position has passed the range of this inverted
            # repeat, there is no point in investigating the other
            # positions as none will fall in this range
            if(position in mappedDict[IRStrand]):
                # Store the value of this dictionary in mappedTagList
                mappedTagList = mappedDict[IRStrand][position]
                for mappedTag in mappedTagList:
                    # If the current tag maps to the 5' arm of the
                    # inverted repeat, add the position to the 
                    # IRMappedTags5 dictionary
                    if(position >= start5 and position+len(mappedTag) - 1 
                            <= end5):
                        IRMappedTags5[position] = []

                    # If the current tag maps to the 3' arm of the
                    # inverted repeat, add the position to the 
                    # IRMappedTags3 dictionary
                    elif(position >= start3 and position+len(mappedTag) - 1
                            <= end3):
                        IRMappedTags3[position] = []

                    # If the tag maps to the loop of the IR, we
                    # only need to 
                    elif(position + len(mappedTag) - 1 >= loopStart and 
                            position <= loopEnd):
                        mappedLoopPositions.append(position)

        # If there are tags that map to both stems of the
        # inverted repeat, the inverted repeat will need to be
        # saved as a candidate precursor miRNA
        if(IRMappedTags5 and IRMappedTags3):
            # Get the positions of tags that map to the 5' and 3' ends
            # of the candidate precursors
            position5List = IRMappedTags5.keys()
            position3List = IRMappedTags3.keys()

            # Initialize a list to track positions that are not miRNA
            # or miRNA* candidates
            positionsToPop = []

            # Loop through the positions that have tags mapping to the
            # 5' arm of an inverted repeat to get the sequences and
            # abundances
            for position in position5List:
                # Find the sequences and abundances of all tags that
                # map to an inverted repeat
                totalAbun5 += findOverlappingSeqsAndAbuns(position,
                    mappedDict[IRStrand], end5, libDict,
                    IRMappedTags5)

                # If no tags were added to the IRMappedTags5
                # dictionary at this position, tag the entry for removal
                # to eliminate entries with no mapping sRNAs
                if(not IRMappedTags5[position]):
                    positionsToPop.append(position)

            # Iterate through positionsToPop and remove those positions
            # from IRMappedTags5
            for position in positionsToPop:
                IRMappedTags5.pop(position, None)

            # If all positions were popped from IRMappedTags5,
            # there is no point in investigating the 3' arm because
            # there can't possibly be a miRNA duplex.
            if(not IRMappedTags5):
                continue

            # Reinitialize positionsToPop back to an empty list
            positionsToPop = []

            # Loop through the positions that have tags mapping to the
            # 3' arm of an inverted repeat to get the sequences and
            # abundances
            for position in position3List:
                # Find the sequences and abundances of all tags that
                # map to an inverted repeat
                totalAbun3 += findOverlappingSeqsAndAbuns(position,
                    mappedDict[IRStrand], end3, libDict,
                    IRMappedTags3)

                # If no tags were addded to the IRMappedTags3
                # dictionary at this position, remove the position from 
                # the dictionary
                if(not IRMappedTags3[position]):
                    positionsToPop.append(position)

            # Iterate through positionsToPop and remove those positions
            # from IRMappedTags3
            for position in positionsToPop:
                IRMappedTags3.pop(position, None)

        # If after all of our checks, there are definitely tags
        # that map to both ends of the IR repeat, we will add
        # the precursor to findOverlappingSeqsAndAbuns
        if(IRMappedTags5 and IRMappedTags3):
            # Create a list of tags that map to both arms of the
            # inverted repeat to ensure no double counting occurs
            # when counting loop abundance. This is necessary for
            # when two tags map to the same position, but one is 
            # longer than the other and extends into the loop
            IRMappedSequences = []
            for position, mappedSeqList in IRMappedTags5.items():
                for sequenceAbun in mappedSeqList:
                    IRMappedSequences.append(sequenceAbun[0])
            for position, mappedSeqList in IRMappedTags3.items():
                for sequenceAbun in mappedSeqList:
                    IRMappedSequences.append(sequenceAbun[0])

            # Before finishing with this precursor, we need to
            # determine the abundance of tags that map to the loop
            for position in mappedLoopPositions:
                tagList = [mappedTags for mappedTags in \
                mappedDict[IRStrand][position]]

                # We can't just use findOverlappingSeqsAndAbuns
                # because it is is mapping to an IR arm, so we will
                # just use a component of that here
                for sequence in tagList:
                    # If the sequence is in IRMappedSequences, skip it
                    if(sequence not in IRMappedSequences):
                        abundance = libDict[sequence][0]
                        hits = libDict[sequence][1]
                        hitsNormalizedAbundance = abundance/hits
                        loopAbun += hitsNormalizedAbundance

            # Add the mapping information to the precursor dictionary
            precursorsDict[IRName] = (IRMappedTags5, IRMappedTags3,
                totalAbun5, totalAbun3, loopAbun)

    return(precursorsDict)

def writePrecursors(filename, chrDict, IRDictByChr, mappedTagsToPrecursors):
    """
    Write the predicted precursors to a file with the designated file name

    Args:
        filename: The name of the output file
        chrDict: Dictionary of chromosomes and their corresponding indices
            in the precursorsList
        IRDictByChr: List of dictionaries with the inverted repeat information
        mappedTagsToPrecursors: List of dictionaries containing all tags 
            that map to each chromosome of the genome. Each list is ui
            dictionary 

    """

    # Open the output file
    with open(filename, 'w') as f:
        for chrName in sorted(chrDict.keys()):
            chrIndex = chrDict[chrName]

            mappedTagsDict = mappedTagsToPrecursors[chrIndex]

            # Within each chromosome, loop through each precursor
            # and extract information to write to the file
            for precursorName, precursor in mappedTagsDict.items():
                # Coordinates are in the 0th element of the tuple,
                # and the dictioniary of tags mapping to the 5' and 3'
                # arms are in the 1st and 2nd element of the tuple,
                # respectively
                coordinates = IRDictByChr[chrIndex][precursorName]
                arm5 = precursor[0]
                arm3 = precursor[1]

                # Write the name and thechromosome number
                # followed by the coordinates of the precursor
                for i in range(len(coordinates)):
                    if(i == 1):
                        f.write('%s,' % chrName)
                    if(i == 6 or i == 7 or i == 8 or i == 9):
                        f.write('%s\n' % coordinates[i])
                    else:
                        f.write('%s,' % coordinates[i])

                # Loop through the positions in sorted order to write
                # the mapping tags on the 5' arm to the file
                for position in sorted(arm5.keys()):
                    # Because more than 1 tag can map to the same
                    # postiion, loop through the list of entries to
                    # write them all to the file
                    for entry in arm5[position]:
                        tag = entry[0]
                        abundance = entry[1]
                        f.write('%s,%s,%s\n' % (tag, position, abundance))

                # Loop through the positions in sorted order to write
                # the mapping tags on the 3' arm to the file
                for position in sorted(arm3.keys()):
                    # Because more than 1 tag can map to the same
                    # postiion, loop through the list of entries to
                    # write them all to the file
                    for entry in arm3[position]:
                        tag = entry[0]
                        abundance = entry[1]
                        f.write('%s,%s,%s\n' % (tag, position, abundance))

def findSequenceInIR(sequence, IRArm, tagLength):
    """  
    If a sequence cannot be found in a simple search in an inverted repeat
    arm, it is because there are gaps that interrupt the sequence. Thus, we
    must do some work to find the tag, WITH gaps, as wel as the start and
    end positions of the tag with gaps to identify potential duplexes

    Args:
        sequence: Small RNA to be idenfitied in the inverted repeat
        IRArm: The aligning arm of the inverted repeat with gaps as hyphens
        tagLength: The length of the sequence
        
    Returns:
        A tuple of the sequence with gaps, its start position, and its
        end position on the inverted repeat arm

    """

    # Get the start position of the sequence on the IR arm without any
    # gap sequences
    localStart = IRArm.replace('-','').find(sequence)

    # Initialize offset and sequence gap count to 0
    offset = 0
    gapCount = 0

    # Initialize offset update flag to True
    offsetUpdate = True

    # Loop through until no update is made to the offset through
    # successive iterations
    while(offsetUpdate):
        # Store the previous offset value for update determination
        oldOffset = offset

        # Set the offset to the number of gaps prior to the 5' start
        # position + the current offset (important because if 
        # is a gap after the previously found offset)
        offset = IRArm[:localStart + \
            offset].count('-')

        # If the new offset is the same as the
        # old (number of gaps), set the
        # offsetUpdate flag to False to break
        if(offset == oldOffset):
            offsetUpdate = False

    # Initialize non gap base counts
    baseCount = 0

    # Initialize a variable to hold the sequence with gaps
    sequenceWithGaps = ''

    # Initialize a counter to 0 so that we can
    # fill the new sequence variables one by one
    counter = 0

    # Loop through the arm of the IR to fill
    # sequence5 with its sequence, including
    # gaps
    while(baseCount < tagLength):
        # Store the current base as the current
        # nucleotide at this calculated position
        base = IRArm[localStart + offset + counter]

        # Add the base to the new sequence
        sequenceWithGaps += base

        # If the base is not a gap, increment the
        # baseCount variable to ensure we only
        # record as many bases to replicate the
        # mapped sRNA
        if(base != '-'):
            baseCount += 1

        else:
            gapCount += 1

        # Increment the counter to ensure we
        # progress through the arm
        counter += 1

    localStart += offset
    localEnd = localStart + tagLength + gapCount - 1

    return(sequenceWithGaps, localStart, localEnd)

def getAlign(base1, base2):
    """
    This function gets the alignment between two bases

    Args:
        base1: The first nucleotide to be compared
        base2: The second nucleotide to be compared
    Returns:
        Match, mismatch, or gap

    """

    # Block to check if the alignment is a gap
    if(base1 == '-' or base2 == '-'):
        return("gap")

    # Block to check for simple mismatches (ie anything that cannot
    # possibly be a G-U wobble
    if(base1 == 'A' and base2 != 'T'):
        return("mismatch")
    elif(base1 == 'C' and base2 != 'G'):
        return("mismatch")

    # Block to check for a G-C match or G-U wobble
    elif(base1 == 'G'):
        if(base2 == 'C'):
            return("match")
        elif(base2 == 'T'):
            return("wobble")
        else:
            return("mismatch")

    # Block to check for U-A match or U-G wobble
    elif(base1 == 'T'):
        if(base2 == 'A'):
            return("match")
        elif(base2 == 'G'):
            return("wobble")
        else:
            return("mismatch")

    # If no other block was a hit, then we have a match
    else:
        return("match")

def getAlignment(arm5, arm3, alignStart, alignEnd):
    """
    This function gets the alignment that was agenerated by einverted
    for the candidate miRNA and miRNA* sequences. It will use the 5' and
    3' arms of the inverted repeats, and the previously determined 
    alignStart and alignEnd to determine which positions are the candidates

    Args:
        arm5: The 5' arm of the inverted repeat, generated by einverted
        arm3: The 3' arm of the inverted repeat, generated by einverted
        alignStart: The start position of the alignment
        alignEnd: The end position of the alignment

    """

    # Initialize variables to count the alignment informatoin
    matchCount = 0
    mismatchCount = 0
    wobbleCount = 0
    gapCount = 0

    # Loop through each position of the alignment to determine if the
    # position is a match, mismatch, wobble, or gap
    for position in range(alignStart, alignEnd):
        # Get the aligning information at this position from the
        # getAlign function
        align = getAlign(arm5[position], arm3[position])

        # If the alignment between the two bases is a match, increment
        # the number of matches and close the gap if it had been 
        # previously opened
        if(align == "match"):
            matchCount += 1

        # If the alignment between the two bases is a mismatch,
        # increment the number of mismatches 
        elif(align == "mismatch"):
            mismatchCount += 1

        elif(align == "wobble"):
            wobbleCount += 1

        else:
            gapCount += 1

    return(matchCount, mismatchCount, wobbleCount, gapCount)

def getVariantAbundance(mappedTagsDict, candidateSequence,
    candidatePosition, strand):
    """
    This function will serve two purposes. First, we need to get the sum of
    abundance of all eight 1 nt variants of the candidate sequence. Second,
    this function will return 0 and thus eliminate a candidate sequence if
    it is not the maximum of all 1 nt variants

    Args:
        mappedTagsDict: A dictionary with key of position and valuse of 
            tuples containing sequences and abundances of those sequences
        candidateTagAndAbundance: The mapped candidate sequence
        candidatePosition: The position on the mappin chromosome of the
            candidate
    Returns:
        A list of variant abundances for the input sequence

    """

    candidateSequenceLength = len(candidateSequence)
    variantAbundanceList = []

    if(strand == 'w'):
        # Check if any variants exist in which the 5' is tailed by 1 nt
        if(candidatePosition - 1 in mappedTagsDict):
            # Loop through all tags that map to this location to identify
            # any potential variants for this case
            for mappedTag in mappedTagsDict[candidatePosition - 1]:
                mappedSequence = mappedTag[0]

                # Check for a variant if the mapped sequence is extended on
                # both the 5' end on the 3' end
                if(len(mappedSequence) == candidateSequenceLength + 2):
                    if(mappedSequence[1:-1] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is extended on
                # the 5' end unchanged on the 3' end
                elif(len(mappedSequence) == candidateSequenceLength + 1):
                    if(mappedSequence[1:] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is extended on
                # the 5' end and shortened on the 3' end
                elif(len(mappedSequence) == candidateSequenceLength):
                    if(mappedSequence[1:] == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])

        # Check if any variants exist in which the 5' is unchanged
        if(candidatePosition in mappedTagsDict):
            # Loop through all tags that map to this location to identify
            # any potential variants for this case
            for mappedTag in mappedTagsDict[candidatePosition]:
                mappedSequence = mappedTag[0]

                # Check for a variant if the mapped sequence is unchanged on
                # the 5' end and extended on the 3' end
                if(len(mappedSequence) == candidateSequenceLength + 1):
                    if(mappedSequence[:-1] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])
            
                # Check for a variant if the mapped sequence is unchanged on
                # the 5' end and shortened on the 3' end
                if(len(mappedSequence) == candidateSequenceLength) - 1:
                    if(mappedSequence == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])

        # Check if any variants exist in which the 5' is shortened
        if(candidatePosition + 1 in mappedTagsDict):
            # Loop through all tags that map to this location to identify
            # any potential variants for this case
            for mappedTag in mappedTagsDict[candidatePosition + 1]:
                mappedSequence = mappedTag[0]

                # Check for a variant if the mapped sequence is shortened on
                # the 5' end and extended on the 3' end
                if(len(mappedSequence) == candidateSequenceLength):
                    if(mappedSequence[:-1] == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is shortened on
                # the 5' end and unchanged on the 3' end
                if(len(mappedSequence) == candidateSequenceLength - 1):
                    if(mappedSequence == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is shortened on
                # the 5' end and shortened on the 3' end
                if(len(mappedSequence) == candidateSequenceLength - 2):
                    if(mappedSequence == candidateSequence[1:-1]):
                        variantAbundanceList.append(mappedTag[1])

    # If the strand is c, mapping is reverse complelmented and thus 
    # the indexing for the variant calculations will basically be the
    # opposite of above
    else:
        # Check if any variants exist in which the 3' is tailed
        if(candidatePosition - 1 in mappedTagsDict):
            # Loop through all tags that map to this location to identify
            # any potential variants for this case
            for mappedTag in mappedTagsDict[candidatePosition - 1]:
                mappedSequence = mappedTag[0]

                # Check for a variant if the mapped sequence is extended on
                # both the 3' end on the 5' end
                if(len(mappedSequence) == candidateSequenceLength + 2):
                    if(mappedSequence[1:-1] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is extended on
                # the 3' end unchanged on the 5' end
                elif(len(mappedSequence) == candidateSequenceLength + 1):
                    if(mappedSequence[:-1] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is extended on
                # the 3' end and shortened on the 5' end
                elif(len(mappedSequence) == candidateSequenceLength):
                    if(mappedSequence[:-1] == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])

        # Check if any variants exist in which the 3' is unchanged
        if(candidatePosition in mappedTagsDict):
            # Loop through all tags that map to this location to identify
            # any potential variants for this case
            for mappedTag in mappedTagsDict[candidatePosition]:
                mappedSequence = mappedTag[0]

                # Check for a variant if the mapped sequence is unchanged on
                # the 3' end and extended on the 5' end
                if(len(mappedSequence) == candidateSequenceLength + 1):
                    if(mappedSequence[1:] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])
            
                # Check for a variant if the mapped sequence is unchanged on
                # the 3' end and shortened on the 5' end
                if(len(mappedSequence) == candidateSequenceLength) - 1:
                    if(mappedSequence == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])

        # Check if any variants exist in which the 3' is shortened
        if(candidatePosition + 1 in mappedTagsDict):
            # Loop through all tags that map to this location to identify
            # any potential variants for this case
            for mappedTag in mappedTagsDict[candidatePosition + 1]:
                mappedSequence = mappedTag[0]

                # Check for a variant if the mapped sequence is shortened on
                # the 3' end and extended on the 5' end
                if(len(mappedSequence) == candidateSequenceLength):
                    if(mappedSequence[1:] == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is shortened on
                # the 3' end and unchanged on the 5' end
                if(len(mappedSequence) == candidateSequenceLength - 1):
                    if(mappedSequence == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])

                # Check for a variant if the mapped sequence is shortened on
                # the 3' end and shortened on the 5' end
                if(len(mappedSequence) == candidateSequenceLength - 2):
                    if(mappedSequence == candidateSequence[1:-1]):
                        variantAbundanceList.append(mappedTag[1])

    if(not variantAbundanceList):
        return([0])

    else:
        return(variantAbundanceList)

def filterPrecursors(mappedTagsToPrecursors, IRDict, overhang):
    """
    This function will perform the sRNA mapping and abundance filters.
    It will first try to find a miRNA and miRNA* pair by identifying
    tags that map to opposite sides of the precursor. It will also create
    splits of the c and w strand if there are tags that map to both

    Args:
        mappedTagsToPrecursors: Dictionary of tag information mapping
            to the precursor, identified by the precursor name
        IRDict: Dictionary of the inverted repeats in one chromosome
        overhang: Maximum length of overhang that a duplex can have

    """

    # Though the intention is to remove all candidates that do not meet 
    # all filter requirements, it is necessary to report to the user why
    # certain precursors are removed. Thus, we create a variable here
    # to contain all precursors that have an identifiable miRNA:miRNA*
    # duplex prior to abundance filters. This means a duplex that
    # fails alignment filters will not be contained here
    precursorsWithDuplex = {}
    finalCandidates = {}
    uniqueCandidatesList = {}

    # Begin to loop though all of the candidate precursors for the
    # various filters. Each loop begins on the chromosome dictionary
    for precursorName, mappedTagsTuple in  mappedTagsToPrecursors.items():
        # Initialize a flag for if the 5' or 3' end of the precursor
        # contains a candidate miRNA
        is3Candidate = False
        is5Candidate = False

        precursor = IRDict[precursorName]
        # Store various elements of the precursor dictionary values
        # for quick accession
        start5 = precursor[0]
        end5 = precursor[1]
        start3 = precursor[2]
        end3 = precursor[3]
        strand = precursor[5]
        arm5 = precursor[6]
        alignmentIndicators = precursor[7]
        arm3 = precursor[8]

        # Store the various elements of the mapped tags tuple
        mappedTagsDict5 = mappedTagsTuple[0]
        mappedTagsDict3 = mappedTagsTuple[1]
        totalAbun5 = mappedTagsTuple[2]
        totalAbun3 = mappedTagsTuple[3]
        loopAbun = mappedTagsTuple[4]

        # Begin a series of loops to identify if there are any tags on
        # the 5' and 3' strands that overlap within a short, user defined
        # overhang
        for candidate5Pos, mappedTagList5 in mappedTagsDict5.items():
            for mapped5Tag in mappedTagList5:
                # Get the length of the 5' candidate tag so that we can
                # determine local positions on the precursors
                tag5Length = len(mapped5Tag[0])
                tag5Abun = mapped5Tag[1]

                # If the length of the tag is not between 20 and 24,
                # just move to the next tag. We do this here because
                # we have to store tags that are 1 nt variants of
                # candidate miRNA or miRNA* sequences
                if(tag5Length < 20 or tag5Length > 24):
                    continue

                # If the strand is w, the sequence will require no
                # modifications
                if(strand == 'w'):
                    sequence5 = mapped5Tag[0]

                # If the strand is c, we need to reverse complement
                # the mapped sequence so that we can find it on the
                # IR arm
                else:
                    sequence5 = mapped5Tag[0].translate(
                        str.maketrans("ACGT","TGCA"))[::-1]

                oldSequence5 = sequence5

                # If we are unable to find the sequences in the
                # IR arm, we know it is for one of two     
                # possibilities. Because we allow 1 mismatch in 
                # bowtie by default, we know that it is possilbe
                # that the sRNA that is mapped to this position
                # may not be exact with the IR. Also, there can
                # be a gap in the alignment, so we must identify
                # which case (if not both) it is before proceeding
                if(sequence5 not in arm5):
                    sequence5, local5Start, local5End = \
                        findSequenceInIR(sequence5, arm5, tag5Length)

                # If the sequence can be found, update the local positions
                # as they may be shifted due to gaps prior
                else:
                    local5Start = arm5.find(sequence5)
                    local5End = local5Start + tag5Length - 1

                # Check to confirm that the sequence with gaps is the
                # same sequence as before
                if(oldSequence5 != sequence5.replace('-','')):
                    print("findSequenceInIR messed up for %s. "\
                        "Contact Reza to debug" % oldSequence5)
                    print(precursorName, oldSequence5, sequence5, local5Start, local5End)
                    sys.exit()

                # Loop through all mapped tags in the 3' dictionary to
                # identify any candidate miRNA:miRNA* pairs with the
                # current 5' mapped tag
                for candidate3Pos, mappedTagList3 in mappedTagsDict3.items():
                    for mapped3Tag in mappedTagList3:
                        # Get the length of the 3' candidate tag so that
                        # we can determine local positions on the
                        # precursor for mapping comparisons. A candidate
                        # will be recorded if a miRNA:miRNA* pair can
                        # be identified
                        tag3Length = len(mapped3Tag[0])
                        tag3Abun = mapped3Tag[1]

                        # If the length of the tag is not between 20 and
                        # 24, just move to the next tag. We do this here
                        # because we have to store tags that are 1 nt
                        # variants of candidate miRNA or miRNA* sequences
                        if(tag3Length < 20 or tag3Length > 24):
                            continue

                        # If the strand is w, the sequence needs to be
                        # reversed because it is on the 3' arm of the IR
                        if(strand == 'w'):
                            sequence3 = mapped3Tag[0][::-1]

                        # If the strand is c, the sequence needs to be
                        # complemented (but not reversed) because it is
                        # on the 3' arm of the IR, but the reverse strand
                        # of the genome
                        else:
                            sequence3 = (mapped3Tag[0].translate(
                                str.maketrans("ACGT","TGCA")))

                        oldSequence3 = sequence3

                        # If we are unable to find the sequences in the
                        # IR arm, we need to find the alignment sequence,
                        # start, and end positions
                        if(sequence3 not in arm3):
                            sequence3, local3Start, local3End = \
                                findSequenceInIR(sequence3, arm3, tag3Length)

                        else:
                            local3Start = arm3.find(sequence3)
                            local3End = local3Start + tag3Length - 1

                        # Check to confirm that the sequence with gaps is
                        # the same sequence as before
                        if(oldSequence3 != sequence3.replace('-','')):
                            print("findSequenceInIR messed up for %s. "\
                                "Contact Reza to debug" % oldSequence3)
                            print(precursorName, oldSequence3, sequence3, local3Start, local3End)
                            sys.exit()

                        # If there is a 3' overhang on either the sequence,
                        # we have a candidate duplex and will investigate
                        # it further
                        if((local5End - local3End == overhang) and (
                                local5Start - local3Start == overhang)):
                            # Because we can have overhangs, the alignment
                            # should start and end at the postiions just
                            # prior to the overhang
                            alignStart = max(local5Start, local3Start)
                            alignEnd = min(local5End, local3End)

                            # Get the einverted alignment for the two
                            # sequences
                            matchCount, mismatchCount, wobbleCount,\
                                gapCount = getAlignment(arm5, arm3,
                                alignStart, alignEnd)
                            
                            # If the alignment between the overlapping
                            # regions of the miRNA and miRNA* do not
                            # exceed our alignment filters, add the
                            # duplex to the precursorsWithDuplex filter
                            if(gapCount + mismatchCount + (wobbleCount * .5) 
                                   <= 5 and gapCount <= 3):

                                ### Code for the abundance filter
                                #variant5Abun = totalAbun5
                                #variant3Abun = totalAbun3
                                variant5Abun = tag5Abun
                                variant3Abun = tag3Abun

                                # Get the abundance of all eight 1-nt
                                # variants of both 5' and 3' tags
                                variant5AbunList = getVariantAbundance(
                                    mappedTagsDict5, mapped5Tag[0],
                                    candidate5Pos, strand)
                                variant3AbunList = getVariantAbundance(
                                    mappedTagsDict3, mapped3Tag[0],
                                    candidate3Pos, strand)

                                if(tag5Abun < max(variant5AbunList) or
                                        tag3Abun < max(variant3AbunList)):
                                    continue

                                if(variant5Abun == -1 or 
                                        variant3Abun == -1):
                                    continue

                                variant5Abun += sum(variant5AbunList)
                                variant3Abun += sum(variant3AbunList)

                                # Get the proportion of reads coming from
                                # the miRNA duplex compred to the rest
                                # of the reads mapping to the duplex
                                proportion = (variant5Abun +
                                    variant3Abun) / (totalAbun5
                                    + totalAbun3 + loopAbun)

                                if(tag5Abun >= tag3Abun):
                                    duplex = ["5p", mapped3Tag[0],
                                        candidate5Pos, candidate3Pos, tag5Abun,
                                        tag3Abun, matchCount, mismatchCount,
                                        wobbleCount, gapCount, variant5Abun,
                                        variant3Abun, totalAbun5, totalAbun3,
                                        loopAbun, proportion]

                                    # Add the precursor name as a key to
                                    # precursorsWithDuplex if it does not
                                    # yet exist. The valu will be a list of 
                                    # duplexes found in the precursor, but the
                                    # first element will be the IR coordinates
                                    if(precursorName not in precursorsWithDuplex):
                                        precursorsWithDuplex[precursorName] = {}

                                    precursorsWithDuplex[precursorName][\
                                        mapped5Tag[0]] = duplex

                                    # If the sum of the two tags in the
                                    # make up more than 75% of the read
                                    # abundance in the entire precursor,
                                    # add the duplex to the candidates
                                    # dictionary
                                    if(proportion >= .75):
                                        # Add the precursor name as a key to
                                        # finalCandidates if it does not
                                        # yet exist. The valu will be a list of 
                                        # duplexes found in the precursor, but the
                                        # first element will be the IR coordinates
                                        if(precursorName not in finalCandidates):
                                            finalCandidates[precursorName] = {}

                                        finalCandidates[precursorName][\
                                            mapped5Tag[0]] = duplex

                                if(tag3Abun >= tag5Abun):
                                    duplex = ["3p", mapped5Tag[0],
                                        candidate3Pos, candidate5Pos, tag3Abun,
                                        tag5Abun, matchCount, mismatchCount,
                                        wobbleCount, gapCount, variant3Abun,
                                        variant5Abun, totalAbun3, totalAbun5,
                                        loopAbun, proportion]

                                    # Add the precursor name as a key to
                                    # precursorsWithDuplex if it does not
                                    # yet exist. The valu will be a list of 
                                    # duplexes found in the precursor, but the
                                    # first element will be the IR coordinates
                                    if(precursorName not in precursorsWithDuplex):
                                        precursorsWithDuplex[precursorName] = {}

                                    precursorsWithDuplex[precursorName][\
                                        mapped3Tag[0]] = duplex

                                    # If the sum of the two tags in the
                                    # make up more than 75% of the read
                                    # abundance in the entire precursor,
                                    # add the duplex to the candidates
                                    # dictionary
                                    if(proportion >= .75):
                                        # Add the precursor name as a key to
                                        # finalCandidates if it does not
                                        # yet exist. The valu will be a list of 
                                        # duplexes found in the precursor, but the
                                        # first element will be the IR coordinates
                                        if(precursorName not in finalCandidates):
                                            finalCandidates[precursorName] = {}

                                        finalCandidates[precursorName][
                                            mapped3Tag[0]] = duplex

    return(precursorsWithDuplex, finalCandidates)

def mergeLibSeqs(LibList):
    """
    Create a list of unique tags for the purpose of mapping unique
    sequences to the genome

    Args:
        libDictList: A list of dictionaries of each library

    Returns:
        List of unique sequences

    """

    # Create an empty list for the unique sequences
    allSeqs = []

    # Loop through each Library object and add all sequences to uniqueSeqs
    for LibClass in LibList:
        allSeqs.extend(list(LibClass.libDict.keys()))

    # Create a set of the allSeqs to get a list of unique sequences so
    # that they may be mapped with bowtie
    uniqueSeqs = list(set(allSeqs))

    # Create a libs folder if it does not exist yet
    if(not os.path.isdir("libs")):
        os.mkdir("libs")

    # Open a file to write the unique sequences to and create a counter
    # to name the sequences
    f = open("libs/uniqueSeqs.fa", "w")
    counter = 0

    # Write the unique sequences to a file so that they can be mapped
    # to the inverted repeats together
    for sequence in uniqueSeqs:
        f.write(">s_%s\n%s\n" % (counter, sequence))
        counter += 1

    return(uniqueSeqs)

def writeCandidates(filename, candidatesByLibDict, filteredPrecursorsDict,
        IRDictByChr, libFilenamesList, chrDict):
    """
    Write the candidate miRNAs to the provided filename

    Args:
        filename: Name of file to write data to (in both csv and fasta format)
        candidatesByLibDict: A dictionary of miRNA precursors that have been
            validated in any library as a key, and a list of library file 
            names as keys. This list must contain more than one file to be
            a candidate
        filtredPrecursorsDict: A dictionary of more detailed information
            on the precursors. Candidates are separated by chromosome.
        IRDictByChr: Dictionary of the inverted repeats and their coordinates
        libFilenamesList: A list of the library filenames used for predictions
        chrDict: Dictionary of chromosome names as the key and their index
            in a list for many variables

    """

    outputFilename = filename + '.csv'
    fastaFilename = filename + '.fa'

    # Open the output files
    with open(outputFilename, 'w') as f, open(fastaFilename, 'w') as g:
        # Write the column names
        f.write("miR Name,Chr,Strand,miR Position,miR Sequence,miR Length,"\
            "Star Position,Star Sequence,Star Length,")

        for libName in libFilenamesList:
            libNameNoFolders =  ''.join(os.path.splitext(
                libName.split('/')[-1])[0])
            f.write("miR Abun in {0},Star Abun in {0},1-nt Variants miR Abun "\
                "in {0},1-nt Variants Star Abun in {0},Total Precursor Abun "\
                "in {0},Proportion of reads from miR:miR* in {0}".format(
                libNameNoFolders))

            if(libName == libFilenamesList[-1]):
                f.write("\n")
            else:
                f.write(",")

        # Loop through candidates and identify any that have been
        # identified in more than one library. These should be written
        # to the output files
        for chrName, precursorDict in candidatesByLibDict.items():
            # Get the chrIndex to pull from the proper index of the
            # data structures separating chromosomes as a base list
            chrIndex = chrDict[chrName]

            for precursorName, candidatesDict in precursorDict.items():
                coordinates = IRDictByChr[chrIndex][precursorName]
                strand = coordinates[5]

                mirCount = 1

                for mirSeq, libList in candidatesDict.items():
                    # We need to begin writing all information that is not
                    # library specific, so we will pull this from the first
                    # library in the list before looping through all libraries
                    # and populating their respective columns in the output
                    # file
                    firstLib = libList[0]
                    duplexInfo = filteredPrecursorsDict[firstLib][chrName][
                        precursorName][mirSeq]
                    arm = duplexInfo[0]
                    starSeq= duplexInfo[1]
                    mirPos = duplexInfo[2]
                    starPos = duplexInfo[3]

                    # If there is one miRNA candidate in the duplex,
                    # we don't need to create a unique ID other than
                    # the name of the precursor
                    if(len(candidatesDict.keys()) == 1):
                        mirName = 'miR%s-%s' % (precursorName.split(
                            'mir')[1], arm)

                    # If there is more than one miRNA candidate in the
                    # duplex, append mirCount to the precursor number
                    # to create a unique identifier for each miRNA
                    # in the duplex
                    else:
                        mirName = 'miR%s_%s-%s' % (precursorName.split(
                            'mir')[1], mirCount, arm)
                        mirCount += 1

                    # Write the candidate name and the sequence to the
                    # fasta file
                    g.write('>%s\n' % mirName)
                    g.write("%s\n" % mirSeq)

                    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s," % (mirName,
                        chrName, strand, mirPos, mirSeq, len(mirSeq),
                        starPos, starSeq, len(starSeq)))

                    for libName in libFilenamesList:
                        libNameNoFolders =  ''.join(os.path.splitext(
                            libName.split('/')[-1])[0])
                        # If the sequence was predicted in the current
                        # library, fill its columns with the relevant
                        # data
                        if(libNameNoFolders in libList):
                            duplexInfo = filteredPrecursorsDict[\
                                libNameNoFolders][chrName][precursorName][\
                                mirSeq]
                            mirAbun = duplexInfo[4]
                            starAbun = duplexInfo[5]
                            variantMirAbun = duplexInfo[10]
                            variantStarAbun = duplexInfo[11]
                            totalAbun5 = duplexInfo[12]
                            totalAbun3 = duplexInfo[13]
                            loopAbun = duplexInfo[14]
                            variantAbun = variantMirAbun + variantStarAbun
                            totalAbun = totalAbun5 + totalAbun3 + loopAbun
                            proportion = duplexInfo[15]

                            f.write("%s,%s,%s,%s,%s,%s" % (mirAbun, starAbun,
                                variantMirAbun, variantStarAbun, totalAbun,
                                proportion))

                        # If the sequence was not predicted in the current
                        # library, write 0 for all columns
                        else:
                            f.write("0,0,0,0,0,0")

                        if(libName == libFilenamesList[-1]):
                            f.write("\n")
                        else:
                            f.write(",")

def getFastaDate(filename, directory):
    """Get the date of the input filename creation date

    Args:
        filename: The name of the subject sequence filaneme
        directory: The directory where the BLAST database would be stored

    Returns:
        The date of the blast DB update

    """

    # Get the time in the struct_time format
    updateTime = datetime.fromtimestamp(os.path.getmtime
        (filename))

    return(updateTime)

def getLocalDBDate(dbName, directory):
    """Get the date of the blast DB creation date

    Args:
        dbName: The name of the blast DB
        directory: The directory where the BLAST database would be stored

    Returns:
        The date of the blast DB update

    """

    # Only check the creation date if the file even exists
    if(os.path.isfile('%s.nhr' % dbName)):
        # Get the time in the struct_time format
        localDatabaseTime = datetime.fromtimestamp(os.path.getmtime
            ('%s.nhr' % dbName))

        return(localDatabaseTime)

    # File doesn't exist so return year 1 for comparison sake.
    # Note: datetime.MINYEAR threw error hence the hardcoded date
    else:
        print("There is no file with name %s. Creating BLAST database."\
            % dbName)
        return(datetime(1,1,1))

def checkNeedUpdateByDate(subjectSequencesFilename, dbName, directory):
    """Check to see if the BLAST database needs to be updated

    Args:
        : The name of the subject sequences file
        dbName: The name of the BLAST database
        directory: The directory where the BLAST database would be stored

    Returns:
        Bool value indicating if the BLAST DB needs to be updated

    """

    # Get the date that the FASTA file was last updated
    fastaFileDate = getFastaDate(subjectSequencesFilename, directory)

    # Get the date that the database was updated last
    localDatabaseDate = getLocalDBDate(dbName, directory)

    # If the expression database is newer than the BLAST database,
    # an update is needed
    if(fastaFileDate > localDatabaseDate):
        print("An update is needed to the BLAST database and will be created "\
        "now. This will take time so please be patient.")
        return(True)

    # If the expression database is older than the BLAST database,
    # no update is needed
    else:
        print("No update is needed to the BLAST database. Proceeding with "\
        "the BLAST.")
        return(False)


def createBlastDB(subjectSequencesFilename, dbFilename):
    """Create a local DB with the subject sequences for BLAST later

    Args:
        subjectSequencesFilename: The name of the file FASTA file holding 
            all the sequenes to be BLASTed against
        dbFilename: The name of the database file that will be written to

    """

    # Create a blast DB from the small RNAs
    os.system("makeblastdb -in %s -dbtype nucl -parse_seqids -out %s" %
        (subjectSequencesFilename, dbFilename))

def blastMirnas(subjectSequencesFilename, dbFilename,
        candidateSequencesFilename):
    """Create a subject database for the known miRNAs and BLAST
       all candidate miRNAs to those miRNAs to find the evolutionarily
       conserved miRNAS

    Args:
        flastFilename: The name of the FASTA file holding the subject
            sequences
        dbFilename: The name of the database file that will be written to
        candidateSequencesFilename: The name of the FASTA file holding the
            candidate sequences
        outputFilename: The name of the output file
    Returns:
        The BLAST results in plain text format

    """

    blastFilename = 'mirBase/blastResults.txt'

    # Run blastn-short, but set word size to 11 as 5 is too short IMO
    NcbiblastnCommandline(query=candidateSequencesFilename,
        db=dbFilename, task='blastn', word_size='15', outfmt=6,
        num_threads=8, strand='plus', out=blastFilename)()[0]

    return(blastFilename)

def addSequencesToOutput(querySequencesFilename, subjectSequencesFilename):
    """Add the full sequences to the XML file for each alignment

    Args:
        querySequencesFilename: Name of the file holding the query
            mirnas
        subjectSequencesFilename: Name of the file holding all subject 
            sequences and IDs

    """

    blastFilename = 'mirBase/blastResults.txt'

    # Parse query and subject sequences 
    querySequences = getSequencesFromFasta(querySequencesFilename)
    subjectSequences = getSequencesFromFasta(subjectSequencesFilename)

    # Read entire output file into to be traversed
    blastResults = readFile(blastFilename, '\t')

    # Reopen outputFile to add the sequences to it
    f_out = open(blastFilename, 'w')

    # Loop through the output file to find which row to add the sequence to
    for row in blastResults:
        # Get the query and subject IDs from the row
        queryID = row[0]
        subjectID = row[1]

        # Using the query and subject IDs, get the query and subject
        # sequences from the respective dictionaries
        querySeq = querySequences[queryID].replace('U','T')
        subjectSeq = subjectSequences[subjectID].replace('U','T')

        # Loop through the row to write the new row to the output file
        for entry in row:
            # As long as the entry is not the last entry, print the entry
            # with a tab after
            if entry != row[-1]:
                f_out.write('%s\t' % entry)

                # If the entry is either the queryID or subjectID, print the
                # corresponding sequence after it
                if(entry == queryID):
                    f_out.write('%s\t' % querySeq)

                if(entry == subjectID):
                    f_out.write('%s\t' % subjectSeq)

            # If the entry is the last element in the row, just write it to the
            # file without a tab
            else:
                f_out.write('%s\n' % entry)

    f_out.close()

def getSequencesFromFasta(filename):
    """Read a FASTA file into memory as a dictionary

    Args:
        filename: Name of the file to be read
    Returns:
        Variable wholeFile containing the entire file

    """

    # Create empty dictionary to hold full FASTA file
    queryMirnas = {}

    # loop through file using csv.reader and store in wholeFile
    f = open(filename, 'r')
    wholeFile = f.readlines()

    # Loop through the miRNAs and save the miRNA tags only
    for i in range(0,len(wholeFile),2):
        # Strip > and newline character from the ID
        seqID = wholeFile[i][1:-1].split(' ')[0].rstrip()
        # Strip the newline character from the sequence
        sequence = wholeFile[i+1]
        # Store the seqID in the dictionary with sequence as definition
        # Strip the newline character from the sequence
        queryMirnas[seqID] = sequence[:-1]

    return(queryMirnas)

def readFile(filename, separator):
    """Read a file into memory

    Args:
        filename: Name of the file to be read
        separator: Delimiter of the file

    Returns:
        Variable wholeFile containing the entire file

    """

    # Create empty array to hold full csv file
    wholeFile = []

    # loop through file using csv.reader and store in wholeFile
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter = separator)
        for row in reader:
            wholeFile.append(row)

    return(wholeFile)

def readBlastResults(filename):
    """Read in the results of the BLAST file after it has had
       its sequences added back to it

    Args:
        filename: The name of the blast file to be read in

    Returns:
        A dictionary of the BLAST results with the candidate miRNA
        name as the key and a dictionary of results as the value

    """

    blastDict = {}

    # Loop through the blast results fil and add the data to a dictionary
    with open(filename) as f:
        for line in f:
            # Remove the new lin character and split the line on tabs
            blastResult = line.split('\n')[0].split('\t')

            # Store each element of the blast result as a separate
            # variable
            queryName = blastResult[0]
            subjectName = blastResult[2]
            # Separate the subject organism identifiers from the miRNA
            # family
            subjectOrganismIdentifier = subjectName.split('-')[0]
            mirnaFamilyName = subjectName.split('-')[1]

            # We don't really care for subsets of a miRNA family, so
            # remove this if it is present
            if(not mirnaFamilyName[-1].isdigit()):
                counter = -1
                while(not mirnaFamilyName[counter].isdigit()):
                    counter -= 1 
                    trimmedMirnaFamilyName = mirnaFamilyName[0:counter+1]

                mirnaFamilyName = trimmedMirnaFamilyName

            # If the query name is not the blast dictionary, initialize
            # it as an empty dictionary
            if(queryName not in blastDict):
                blastDict[queryName] = {}

            # If the mirna family doesn't exist yet in the blastDict
            # subdictionary, initialize the list for organisms before
            # adding to it
            if(mirnaFamilyName not in blastDict[queryName]):
                blastDict[queryName][mirnaFamilyName] = []

            blastDict[queryName][mirnaFamilyName].append(blastResult)

    return(blastDict)

def createSimilarityDict(blastDict, threeLetterIdentifier):
    """Create a similarity dictionary to assist in the annotation of the
       candidate miRNAs that either sho w the equivalence of the candidate
       sequence with a known sequence, or that a candidate miRNA is
       highly similar to a known family in a number of species

    Args:
        blastDict: A dictionary of the blast results with the candidate
            miRNA name as the key and a dictionary of similar known
            miRNAs as the value. This dictionary has the miRNA family
            as a key and then a list of blast results as the value
        threeLetterIdentifier: First letter of genus and first 2 
            letters of species..

    Returns:
        A dictionary that showing the equivalence or near equivalence of
        a candidate miRNA to known miRNAs

    """

    similarityDict = {}

    for candidateMirna, candidateBlastDict in blastDict.items():
        breakFlag = 0
        for mirnaFamily, blastResultsList in candidateBlastDict.items():
            for blastResult in blastResultsList:
                candidateSeq = blastResult[1]
                subjectName = blastResult[2]
                subjectSeq = blastResult[3]
                subjectOrganismIdentifier = subjectName.split('-')[0]
                percIdentity = float(blastResult[4])
                alignLength = int(blastResult[5])
                numMismatch = int(blastResult[6])
                numGaps = int(blastResult[7])

                numMatches = alignLength - numMismatch - numGaps

                # Right now, we are going to classify anything that has
                # fewer than 5 total differences between the candidate
                # sequence and a miRNA known miRNA as a member of that family
                if(len(candidateSeq) - numMatches < 5):
                    # If candidateMirna or mirnaFamily haven't yet been
                    # initialized in similarityDict, initiailize them
                    if(candidateMirna not in similarityDict):
                        similarityDict[candidateMirna] = {}
                    if(mirnaFamily not in similarityDict[candidateMirna]):
                        similarityDict[candidateMirna][mirnaFamily] = []

                    # If the candidate sequence and subject sequences are the
                    # same (and of course the same organism), then set the
                    # the value of the similarityDict just to the miRNA name
                    # and skip the rest of the blast results for this sequence
                    if(candidateSeq == subjectSeq and threeLetterIdentifier ==
                            subjectOrganismIdentifier):
                        similarityDict[candidateMirna] = subjectName
                        # No need to continue with this candidate miRNA, so
                        # break from the loop to inspect the next candidate
                        breakFlag = 1
                        break

                    # If the exact sequence hasn't been found, add the
                    # three letter organism identifier to the similarity
                    # dictionary for this miRNA family
                    if(subjectOrganismIdentifier not in similarityDict[
                            candidateMirna][mirnaFamily]):
                        similarityDict[candidateMirna][mirnaFamily].append(
                            subjectOrganismIdentifier)

            # If an exact match was found in the same organism, break from
            # the loop for this candidateMirna as we no longer need to
            # search for a match
            if(breakFlag):
                breakFlag = 0
                break

    return(similarityDict)

def annotateCandidates(finalFilename, similarityDict, threeLetterIdentifier):
    """Using the similarityDict, annotate the data in finalFilename with
       the proper miRNA name, with the already existing family(ies) that the
       candidate miRNA likely belongs to, or the conservation of an existing
       family that has not yet been identified in this organism but is 
       present in another. If a candidate miRNA is novel, it MUST have been
       identified in more than one library to stay

    Args:
        finalFilename: Name of file to write data to (in both csv and
            fasta format)
        similarityDict: Dictionary of the candidate miRNAs and their
            potential annotations if there is no equivalence, or just
            simply its known name if it already exists
        threeLetterIdentifier: First letter of genus and first two
            of the species

    """

    outputFilename = finalFilename + '.csv'
    newOutputFilename = finalFilename + '_annotated.csv'

    finalCandidatesFile = readFile(outputFilename, ',')

    f = open(newOutputFilename, 'w')

    header = finalCandidatesFile[0]
    header.append("Classification Flag")

    # Get the index of the first library proportion. Useful for
    # when we need to validate novel tags in more than 1 library,
    # but we won't use it for conserved families already found
    # in this species
    startIterIndex = header.index("Star Length") + 6

    for entry in header:
        if(entry == header[-1]):
            f.write("%s\n" % entry)
        else:
            f.write("%s," % entry)

    for line in finalCandidatesFile[1:]:
        mirName = line[0]
        chrName = line[1]
        strand = line[2]
        libCount = 0
        similarFlag = 0

        # Check if the candidate miRNA name is in the similarityDict
        # to determine if it is completely novel or not
        if(mirName in similarityDict):
            # If the value in similarityDict of the candidate miRNA
            # is simply a string, that means that the candidate miRNA
            # is already present within this organism in mirBase, thus
            # replace the candidate name with the proper annotated name
            if(isinstance(similarityDict[mirName], str)):
                similarFlag = True
                line[0] = similarityDict[mirName]
                line.append("Known")

            # If the candidate miRNA is similar to sequences in mirBase, but
            # not identical to anything in this organism, we will need to
            # write this information to the output file.
            else:
                # Loop through each miRNA family that the candidate sequence
                # had fewer than 5 differences to
                for mirFamily, organismList in similarityDict[mirName].items():
                    # If the organism being studied has a very similar
                    # sequence to one that is already known for this organism
                    # in miRBase, we will tag it as a member of this family
                    if(threeLetterIdentifier in organismList):
                        similarFlag = True
                        line.append("New member of existing family")
                        line.append("%s-%s" % (threeLetterIdentifier,
                            mirFamily))

                    # If the sequence is only similar to miRNAs found in
                    # other organismList, don't tag as a member of any family.
                    # Rather, just write the miRNA family name
                    else:
                        for i in range(startIterIndex, len(line), 6):
                            if(float(line[i])):
                                libCount += 1
                        if(libCount > 1):
                            line.append("Conserved family of the following:")
                            line.append(mirFamily)

                            # Add the list of organismList with this same
                            # miRNA family that met the similarity
                            # requirement
                            toWrite = ""
                            for organism in sorted(organismList):
                                if(organism == organismList[-1]):
                                    toWrite += "%s" % organism
                                else:
                                    toWrite += "%s " % organism

                            line.append(toWrite)

        # If the candidate miRNA had no similar sequence, it is completely
        # novel by our tests and thus requires validation in more than one
        # library. Thus, here we will literate through the results and
        # remove candidate miRNAs that were identified in just one library
        else:
            similarFlag = False
            for i in range(startIterIndex, len(line), 6):
                if(float(line[i])):
                    libCount += 1

            if(libCount > 1):
                line.append("Novel")

        if(similarFlag or (not similarFlag and libCount > 1)):
            for i in range(len(line)):
                if i == len(line) - 1:
                    f.write("%s\n" % line[i])

                else:
                    f.write("%s," % line[i])

def main():
    progStart = time.time()
    LibList = []

    # Required overhang between top and bottom strands of miRNA duplex
    # Hardcoded to 2 here, but in such a way that could technically allow
    # modifications
    overhang = 2

    if(parallel):
        nproc = int(round(accel,1))
        pool = Pool(nproc)

    # Create a path for genome if it does not exist already
    if not os.path.isdir("genome"):
        os.mkdir('genome')
    if(not os.path.isdir("invertedRepeats")):
        os.mkdir("invertedRepeats")

    # Create genome object
    GenomeClass = Genome(genomeFilename)

    ##########################################################################

    ############### Find inverted repeats in genome file #####################

    ##########################################################################


    # If not set to run in parallel, run the sequential function
    # for einverted
    if(runEInvertedFlag):
        # Create an empty list for both the inverted repeat FASTA files
        # and alignment files
        IRFastaFilenamesList = []
        IRAlignmentFilenamesList = []

        # If parallel is set, run einverted using the parallel version
        if(parallel):
            print("Running einverted in parallel")

            res = pool.starmap_async(GenomeClass.einverted,
                zip(GenomeClass.genomeSeqList))

            results = res.get()

            # Loop through the results and add the inverted repeat filenames
            # to their respective lists
            for result in results:
                IRFastaFilenamesList.append(result[0])
                IRAlignmentFilenamesList.append(result[1])

        else:
            print("Running einverted sequentially")

            for entry in GenomeClass.genomeSeqList:
                IRName, IRSeq = GenomeClass.einverted(
                    entry)

                IRFastaFilenamesList.append(IRName)
                IRAlignmentFilenamesList.append(IRSeq)

    # If einverted was not run, set the temp file lists to be just the list
    # of the merged final file so that we can create the IRDictByChr using the
    # previously merged file
    else:
        IRFastaFilenamesList = [GenomeClass.IRFastaFilename]
        IRAlignmentFilenamesList = [GenomeClass.IRAlignmentFilename]

    # Combine the inverted repeat temp files of the einverted runs
    # into one file
    GenomeClass.combineIRTempFiles(IRFastaFilenamesList,
        IRAlignmentFilenamesList)

    #########################################################################

    ######################## Map small RNAs to genome #######################
 
    #########################################################################

    precursorsWithDuplexDictByLib = {}
    filteredPrecursorsDict = {}
    candidatesByLibDict = {}

    # Populate candidatesByLibDict chromosomes with empty dictionaries
    for chrName, chrIndex in GenomeClass.chrDict.items():
        candidatesByLibDict[chrName] = {}

    for libraryFilename in libFilenamesList:
        libNameNoFolders =  ''.join(os.path.splitext(
            libraryFilename.split('/')[-1])[0])
        
        print("Beginning to process %s." % libraryFilename)
        Lib = Library(libraryFilename, GenomeClass.chrDict)

        precursorsWithDuplexDictByLib[libNameNoFolders] = {}
        filteredPrecursorsDict[libNameNoFolders] = {}

        for chrName in sorted(GenomeClass.chrDict.keys()):
            precursorsWithDuplexDictByLib[libNameNoFolders][chrName] = {}
            filteredPrecursorsDict[libNameNoFolders][chrName] = {}


        # Map small RNAs to the genome

        # If the parameter mapFilenames is not set, map the sRNA libraries
        # to the genome 
        if(not mapFilenames):
            mapFilename, logFilename = Lib.mapper(
                GenomeClass.indexFilename)

        else:
            print("Mapping will not be performed for file %s" % Lib.filename)

        print("Creating the mapped list for %s" % Lib.filename)
        funcStart = time.time()

        # Create a dictionary with the sequence of all tags that
        # map to a position on every chromosome
        Lib.createMappedList(GenomeClass.chrDict)

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        print("Time to create the mappedList: %s seconds" % (execTime))

    
        #######################################################################

        ################# Map small RNAs to inverted repeats ##################

        #######################################################################

        print("Mapping sRNAs to the inverted repeats")

        funcStart = time.time()

        mappedTagsToPrecursors = []

        if(parallel):
            # Run mapSRNAsToIRs in parallel
            res = pool.starmap_async(mapSRNAsToIRs, zip(GenomeClass.IRDictByChr,
                Lib.mappedList, repeat(Lib.libDict)))

            mappedTagsToPrecursors = res.get()

        else:
            for i in range(len(GenomeClass.chrDict)):

                mappedTagsToPrecursors.append(mapSRNAsToIRs(
                    GenomeClass.IRDictByChr[i], Lib.mappedList[i], Lib.libDict))

        unfilteredFilename = os.path.splitext(Lib.filename)[0] +\
            '_all_precursors.txt'

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        print("Time to map sRNAs to inverted repeats: %s seconds" % (execTime))

        print("Writing precursors to a file")

        writePrecursors(unfilteredFilename, GenomeClass.chrDict,
            GenomeClass.IRDictByChr, mappedTagsToPrecursors)

        #######################################################################

        ################### Filter precursor candidates #######################

        #######################################################################

        print("Filtering candidate precursors")

        funcStart = time.time()

        # If we are running in parallel, run filterPrecursors in parallel.
        # Parallelize by chromosomes
        if(parallel):
            res = pool.starmap_async(filterPrecursors, zip(mappedTagsToPrecursors,
                GenomeClass.IRDictByChr, repeat(overhang)))

            results = res.get()

            for chrName in sorted(GenomeClass.chrDict.keys()):
                chrIndex = GenomeClass.chrDict[chrName]

                # Get the index of chrDict[chrName]
                precursorsWithDuplexDictByLib[libNameNoFolders][chrName] \
                     = results[chrIndex][0]
                filteredPrecursorsDict[libNameNoFolders][chrName] = \
                    results[chrIndex][1]

        # If running sequentially, run filterPrecursors sequentially
        else:
            for chrName in sorted(GenomeClass.chrDict.keys()):
                # Get the index of each chromosome that will be processed
                # sequentially
                chrIndex = GenomeClass.chrDict[chrName]
                precursorList = mappedTagsToPrecursors[chrIndex]
                IRDict = GenomeClass.IRDictByChr[chrIndex]

                precursorsWithDuplexDictByLib[libNameNoFolders][chrIndex], \
                    filteredPrecursorsDict[libNameNoFolders][chrName] = \
                    filterPrecursors(precursorList, IRDict, overhang)

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        print("Time to map filter inverted repeats: %s seconds" % (execTime))

        ###
        # Prior to writing this library's results, add its miRNAs and
        # corresponding precursors to a dictionary

        # Loop through each chromosome of the final candidates dictionary
        # for this library
        for chrName, subFilteredPrecursorsDict in \
                filteredPrecursorsDict[libNameNoFolders].items():

            # Loop through each duplex in the precursor and add it to the 
            # dictionary tracking which library it has been found in
            for precursorName, duplexDict in subFilteredPrecursorsDict.items():
                if(precursorName not in candidatesByLibDict[chrName]):
                    candidatesByLibDict[chrName][precursorName] = {} 

                for mirCandidate in duplexDict.keys():
                    if(mirCandidate not in candidatesByLibDict[chrName]\
                            [precursorName]):
                        candidatesByLibDict[chrName][precursorName]\
                            [mirCandidate] = []

                    candidatesByLibDict[chrName][precursorName]\
                        [mirCandidate].append(libNameNoFolders)

    # Prior to exiting, we need to write the candidates to a file that
    # are actually identified in more than one library. These are the final
    # candidates
    finalFilename = 'output'

    writeCandidates(finalFilename, candidatesByLibDict, filteredPrecursorsDict,
        GenomeClass.IRDictByChr, libFilenamesList, GenomeClass.chrDict)

    ##########################################################################

    ################### Annotate candidate miRNAs ############################

    ##########################################################################

    # Check to see if the BLAST database needs to be updated
    updateFlag = checkNeedUpdateByDate(subjectSequencesFilename, dbFilename,
        directory)

    # Update the blast database if updateFlag is true
    if(updateFlag):
        # Create the database from the file holding the subject sequences
        localStartTime = time.time()
        dbName = createBlastDB(subjectSequencesFilename, dbFilename)

    # BLAST query miRNAs to known miRNAs
    localStartTime = time.time()
    print("Performing BLAST")
    blastFilename = blastMirnas(subjectSequencesFilename, dbFilename,
        queryMirnasFilename)

    # Add field for the subject and query sequences in the BLAST output
    # because these sequences are not within by default
    print("Adding sequences to output file")
    addSequencesToOutput(queryMirnasFilename, subjectSequencesFilename)

    # Read the BLAST results into a dictionary for quick querying
    blastDict = readBlastResults(blastFilename)

    # Create a dictionary of either the known miRNAs that the candidate
    # miRNAs are equal to, or the miRNA families and species that the
    # candidate miRNA is highly similar to
    similarityDict = createSimilarityDict(blastDict, threeLetterIdentifier)

    # Properly annotate the candidate miRNAs with the data in similarityDict
    annotateCandidates(finalFilename, similarityDict, threeLetterIdentifier)

    progEnd = time.time()
    execTime = round(progEnd - progStart, 2)

    print("Total runtime was %s seconds" % execTime)

def test():
    # Check to see if the BLAST database needs to be updated
    updateFlag = checkNeedUpdateByDate(subjectSequencesFilename, dbFilename, directory)

    # Update the blast database if updateFlag is true
    if(updateFlag):
        # Create the database from the file holding the subject sequences
        localStartTime = time.time()
        dbName = createBlastDB(subjectSequencesFilename, dbFilename)

    # BLAST query miRNAs to known miRNAs
    localStartTime = time.time()
    print("Performing BLAST")
    blastFilename = blastMirnas(subjectSequencesFilename, dbFilename,
        queryMirnasFilename)

    # Add field for the subject and query sequences in the XML output file
    print("Adding sequences to output file")
    addSequencesToOutput(queryMirnasFilename, subjectSequencesFilename)

    blastDict = readBlastResults(blastFilename)

    similarityDict = createSimilarityDict(blastDict, threeLetterIdentifier)

    finalFilename = 'output'

    annotateCandidates(finalFilename, similarityDict, threeLetterIdentifier)

if __name__ == '__main__':
    if parallel == 1:
        accel = int(multiprocessing.cpu_count()*0.50)
    else:
        pass

    main()
