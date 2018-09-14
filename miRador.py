#!/usr/loca/bin/python3

# miRunner is the miRNA prediction script within the miRador package.
# It functions by first identifying inverted repeats within the genome.
# Then, it maps all small RNAs to the identified inverted repeats to create
# a list of candidate miRNAs.

#### Written by Reza Hammond
#### rkweku@udel.edu 

#### Things to do still #####
# Perform check on fragmentation, particularly when there are multiple sequences but very short sequences

######## IMPORT ################
import sys,os,re,subprocess,multiprocessing,datetime,time
import statistics
import string
from multiprocessing import Process, Queue, Pool, current_process
from multiprocessing.connection import wait
from itertools import product, repeat

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

########################## EXECUTION VARIABLES ###############################
# Name of the genome file
genomeFilename = 'genome/Maize/Zea_mays.AGPv4.dna.toplevel.fa'
#genomeFilename = 'genome/fake/fakeGenome_c.fa'
# List of library file names
libFilenames = ['libs/MAIZE_sbsCSHL_sRNA/2421_chopped.txt', 'libs/MAIZE_sbsCSHL_sRNA/2432_chopped.txt', 'libs/MAIZE_sbsCSHL_sRNA/Mz_dcl1_chopped.txt', 'libs/MAIZE_sbsCSHL_sRNA/1682_chopped.txt']
#libFilenames = ['libs/fake/fakeLibs_c.txt']#, 'libs/AT_pub2_sRNA/100_chopped.txt']
#libFilenames = ['libs/AT_pub2_sRNA/724_chopped.txt',  'libs/AT_pub2_sRNA/2202_chopped.txt']
# Optional parameter of merged map files for each library. If these are set,
# mapping the libFilenames will be bypassed
mapFilenames = []
#mapFilenames = ['libs/AT_pub2_sRNA/99_chopped.map']#
# Maximum length over overhang
overhang = 2
# Penalty for a score
gap = 12            ## Default:12 (from emboss website)
# Score for a match
match = 3           ## Default:3
# Score for a mismatch
mismatch = -4       ## Default:-4
# Minimum score for a candidate miRNA precursor
threshold = 50      ## Default:50 - This is the score threshold for cutoff
# Maximum inverted repeat length that can be identified as a candidate
# Default is 300 based on max suggested by Axtell and Meyers (TPC 2018)
maxRepLen = 300
# Flag to run einverted for this genome file
runEInvertedFlag = 1
# Flag to utilize parallelization
parallel = 1        ## Default: 1 (Yes)
nthreads = 8
debug = 0

bowtieBuildPath = os.path.expanduser('~/tools/bowtie-1.2.2/bowtie-build')
bowtiePath = os.path.expanduser('~/tools/bowtie-1.2.2/bowtie')

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

        # Create an empty list in IRList for as many chromosomes
        # exist
        self.IRList = [[] for i in range(len(self.chrDict))]

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

        # If einverted was run, open the output alignment file to write
        # the results to
        if(runEInvertedFlag):
            ## Combine inverted repeats alignments results
            align_out = open(self.IRAlignmentFilename ,'w')

        # Loop through all alignment files and merge into one file
        for filename in IRAlignmentFilenamesList:
            with open(filename) as alignmentFile:
                # Loop through the alignment files to add them to the
                # merged file and add them to IRList

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
                        matches, totalBases = map(int, parsedLine[3].split('/'))
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
                        # list within IRList
                        self.IRList[index].append((start5, end5,
                            start3, end3, loop, 'w', hairpin5,
                            alignmentIndicators, hairpin3))
                        self.IRList[index].append((start5, end5,
                            start3, end3, loop, 'c', hairpin5,
                            alignmentIndicators, hairpin3))

                        if(not failFlag and runEInvertedFlag):
                            for entry in toWriteList:
                                align_out.write(entry)
                        toWriteList = []
                        failFlag = 0

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
                # to the dictionary. There should not be any duplicate sequences
                # in this format, however, 
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

def mapSRNAsToIRs(IRList, mappedDict, libDict, mirnasList):
    """
    Map small RNAs to the inverted repeats. This will first read in
    the map file (merged mapfile if run in parallel) and identify any
    small RNAs which overlap with inverted repeats. Any small RNAs that
    do map will be added to the IRDict with the coordinates of overlap.
    Do not count small RNAs that map to the loop of the precursor

    Args:
        IRList: List of the inverted repeats in one chromosome
        mappedDict: Dictionary of positions and sequences that map to those
            positions in one chromosome
        libDict: Dictionary of all library tags and their abundances

    """

    mirnasCount = 0
    listOfMirnas = []

    # Initialize a dictionary to hold the precursors
    precursorsList = []

    # Loop through the inverted repeats to find all tags within our
    # length cutoffs that map to them
    for invertedRepeat in IRList:
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
                        #if(mappedTag in mirnasList):
                        #    mirnasCount += 1
                        #    listOfMirnas.append(mappedTag)
                        #    print(mappedTag)

                    # If the current tag maps to the 3' arm of the
                    # inverted repeat, add the position to the 
                    # IRMappedTags3 dictionary
                    elif(position >= start3 and position+len(mappedTag) - 1
                            <= end3):
                        IRMappedTags3[position] = []
                        #if(mappedTag in mirnasList):
                        #    mirnasCount += 1
                        #    listOfMirnas.append(mappedTag)
                        #    print(mappedTag)

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

            # Add the precursor to the precursorsList
            precursorsList.append([invertedRepeat,
                IRMappedTags5, IRMappedTags3,
                totalAbun5, totalAbun3, loopAbun])

    return(precursorsList)

def writePrecursors(filename, chrDict, precursorsList):
    """
    Write the predicted precursors to a file with the designated file name

    Args:
        filename: The name of the output file
        chrDict: Dictionary of chromosomes and their corresponding indices
            in the precursorsList
        precursorsList: List containing all candidate precursors

    """

    # We want to easily count the number of precursors, so we will
    # assign them a name with a number here to ensure single core
    # processing and thus no duplicate numbers
    precursorCount = 0

    # Open the output file
    with open(filename, 'w') as f:
        for chrName in sorted(chrDict.keys()):
            chrIndex = chrDict[chrName]

            precursorSubList = precursorsList[chrIndex]

            # Within each chromosome, loop through each precursor
            # and extract information to write to the file
            for precursor in precursorSubList:
                # Coordinates are in the 0th element of the tuple,
                # and the dictioniary of tags mapping to the 5' and 3'
                # arms are in the 1st and 2nd element of the tuple,
                # respectively
                coordinates = precursor[0]
                arm5 = precursor[1]
                arm3 = precursor[2]

                name = "mir%05d" % precursorCount
                precursorCount += 1
                # I want the names of the precursors to have 5 digits,
                # but if we end up having more than 99,999 precursors,
                # we need to change the script a bit. Inform the user
                # and exit
                if(precursorCount > 99999):
                    print("Problem: The number of precursors exceeds "\
                        "the number of precursors Reza thought was "\
                        "possible. Extend name from %05d to %06d")
                    sys.exit()

                # Add the name to the precursor
                precursor.insert(0, name)

                # Write the name and the  chromosome number
                # followed by the coordinates of the precursor
                f.write('%s,%s,' % (name, chrName))
                for i in range(len(coordinates)):
                    if(i == 5 or i == 6 or i == 7 or i == 8):
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
        try:
            base = IRArm[localStart + offset + \
                counter]
        except IndexError:
            print(sequence, localStart, offset, counter, baseCount, tagLength, localStart, localEnd)
            print(IRArm)
            sys.exit()
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
    localEnd = localStart + tagLength + gapCount

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
#                        print(mappedSequence, mappedTag[1], 1)

                # Check for a variant if the mapped sequence is extended on
                # the 5' end unchanged on the 3' end
                elif(len(mappedSequence) == candidateSequenceLength + 1):
                    if(mappedSequence[1:] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 2)

                # Check for a variant if the mapped sequence is extended on
                # the 5' end and shortened on the 3' end
                elif(len(mappedSequence) == candidateSequenceLength):
                    if(mappedSequence[1:] == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 3)

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
#                        print(mappedSequence, mappedTag[1], 4)
            
                # Check for a variant if the mapped sequence is unchanged on
                # the 5' end and shortened on the 3' end
                if(len(mappedSequence) == candidateSequenceLength) - 1:
                    if(mappedSequence == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 5)

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
#                        print(mappedSequence, mappedTag[1], 6)

                # Check for a variant if the mapped sequence is shortened on
                # the 5' end and unchanged on the 3' end
                if(len(mappedSequence) == candidateSequenceLength - 1):
                    if(mappedSequence == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 7)

                # Check for a variant if the mapped sequence is shortened on
                # the 5' end and shortened on the 3' end
                if(len(mappedSequence) == candidateSequenceLength - 2):
                    if(mappedSequence == candidateSequence[1:-1]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 8)

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
#                        print(mappedSequence, mappedTag[1], 1)

                # Check for a variant if the mapped sequence is extended on
                # the 3' end unchanged on the 5' end
                elif(len(mappedSequence) == candidateSequenceLength + 1):
                    if(mappedSequence[:-1] == candidateSequence):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 2)

                # Check for a variant if the mapped sequence is extended on
                # the 3' end and shortened on the 5' end
                elif(len(mappedSequence) == candidateSequenceLength):
                    if(mappedSequence[:-1] == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 3)

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
#                        print(mappedSequence, mappedTag[1], 4)
            
                # Check for a variant if the mapped sequence is unchanged on
                # the 3' end and shortened on the 5' end
                if(len(mappedSequence) == candidateSequenceLength) - 1:
                    if(mappedSequence == candidateSequence[1:]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 5)

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
#                        print(mappedSequence, mappedTag[1], 6)

                # Check for a variant if the mapped sequence is shortened on
                # the 3' end and unchanged on the 5' end
                if(len(mappedSequence) == candidateSequenceLength - 1):
                    if(mappedSequence == candidateSequence[:-1]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 7)

                # Check for a variant if the mapped sequence is shortened on
                # the 3' end and shortened on the 5' end
                if(len(mappedSequence) == candidateSequenceLength - 2):
                    if(mappedSequence == candidateSequence[1:-1]):
                        variantAbundanceList.append(mappedTag[1])
#                        print(mappedSequence, mappedTag[1], 8)

    if(not variantAbundanceList):
        return([0])

    else:
        return(variantAbundanceList)

def filterPrecursors(precursorsList, overhang):
    """
    This function will perform the sRNA mapping and abundance filters.
    It will first try to find a miRNA and miRNA* pair by identifying
    tags that map to opposite sides of the precursor. It will also create
    splits of the c and w strand if there are tags that map to both

    Args:
        precursorsList: List of precursors with the sRNAs that
            map to the precursors
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
    for precursor in precursorsList:
        # Initialize a flag for if the 5' or 3' end of the precursor
        # contains a candidate miRNA
        is3Candidate = False
        is5Candidate = False

        # Store various elements of the precursor dictionary values
        # for quick accession
        name = precursor[0]
        start5 = precursor[1][0]
        end5 = precursor[1][1]
        start3 = precursor[1][2]
        end3 = precursor[1][3]
        strand = precursor[1][5]
        arm5 = precursor[1][6]
        alignmentIndicators = precursor[1][7]
        arm3 = precursor[1][8]
        totalAbun5 = precursor[4]
        totalAbun3 = precursor[5]
        loopAbun = precursor[6]

        # Begin a series of loops to identify if there are any tags on
        # the 5' and 3' strands that overlap within a short, user defined
        # overhang
        for candidate5Pos, mappedTagList5 in precursor[2].items():
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
                    print(name, oldSequence5, sequence5, local5Start, local5End)
                    sys.exit()

                # Loop through all mapped tags in the 3' dictionary to
                # identify any candidate miRNA:miRNA* pairs with the
                # current 5' mapped tag
                for candidate3Pos, mappedTagList3 in precursor[3].items():
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
                            print(name, oldSequence3, sequence3, local3Start, local3End)
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
                                    precursor[2], mapped5Tag[0],
                                    candidate5Pos, strand)
                                variant3AbunList = getVariantAbundance(
                                    precursor[3], mapped3Tag[0],
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

                                duplex = [mapped5Tag[0], mapped3Tag[0],
                                    candidate5Pos, candidate3Pos, tag5Abun,
                                    tag3Abun, matchCount, mismatchCount,
                                    wobbleCount, gapCount,
                                    variant5Abun, variant3Abun, totalAbun5,
                                    totalAbun3, loopAbun, proportion]

                                # Add the precursor name as a key to
                                # precursorsWithDuplex if it does not
                                # yet exist. The valu will be a list of 
                                # duplexes found in the precursor, but the
                                # first element will be the IR coordinates
                                if(name not in precursorsWithDuplex):
                                    precursorsWithDuplex[name] = []

                                precursorsWithDuplex[name].append(duplex)

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
                                    if(name not in finalCandidates):
                                        finalCandidates[name] = []

                                    finalCandidates[name].append(duplex)

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

def writeCandidates(filename, chrDict, precursorsListByChr):
    """
    Write the precursors with their duplex to the output file

    Args:
        filename: The name of the output file
        chrDict: Dictionary of chromosomes and their corresponding
            positions
        precursorsListByChr: List containing candidate precursors and
            their miRNA:miRNA* duplex for each chromosome

    """

    precursorFilename = filename + '_precursors.txt'
    fastaFilename = filename + '.fa'

    # Open the output file
    with open(precursorFilename, 'w') as f, open(fastaFilename, 'w') as g:
        # loop through the chromosomes, sorted in numerical order
        for chrName in sorted(chrDict.keys()):
            chrIndex = chrDict[chrName]

            precursorsDict = precursorsListByChr[chrIndex]
            # Within each chromosome, loop through each precursor
            # and extract information to write to the file
            for precursorName, duplexList in precursorsDict.items():
                precursor = duplexList[0]
                duplexCount = 0

                # Write the precursor (its name, coordinates, and alignment)
                # to the precursor file
                f.write('%s,%s,' % (precursorName, chrName))
                for i in range(len(precursor)):
                    if(i == 5 or i == 6 or i == 7 or i == 8):
                        f.write('%s\n' % precursor[i])
                    else:
                        f.write('%s,' % precursor[i]) 

                # Loop through all duplexes that were identified for a
                # precursor and write the stored information to the file
                for duplex in duplexList:
                    candidate5Seq = duplex[0]
                    candidate3Seq = duplex[1]
                    candidate5Pos = duplex[2]
                    candidate3Pos = duplex[3]
                    candidate5Abun = duplex[4]
                    candidate3Abun = duplex[5]
                    matchCount = duplex[6]
                    mismatchCount = duplex[7]
                    wobbleCount = duplex[8]
                    gapCount = duplex[9]
                    variant5Abun = duplex[10]
                    variant3Abun = duplex[11]
                    totalAbun5 = duplex[12]
                    totalAbun3 = duplex[13]
                    loopAbun = duplex[14]
                    variantAbun = variant5Abun + variant3Abun
                    totalAbun = totalAbun5 + totalAbun3 + loopAbun
                    proportion = duplex[15]

                    # In the event that there is more than one duplex for
                    # this precursor, we need a unique identifier for the
                    # miRNA and miRNA*
                    if(len(duplexList) > 1):
                        mirName = "miR%s_%s" % (precursorName.split(
                            'mir')[1], duplexCount)
                        duplexCount += 1
                    else:
                        mirName = "miR%s" % precursorName.split('mir')[1]

                    # Our official candidate of this precursor will
                    # be the sequence with greatest abundance. If the
                    # sequences have the same abundance, we will go with
                    # the mapped 5 sequence. Maybe this should be
                    # different though...
                    if(candidate3Abun > candidate5Abun):
                        mirSeq = candidate3Seq
                        mirPos = candidate3Pos
                        mirAbun = candidate3Abun
                        starSeq = candidate5Seq
                        starPos = candidate5Pos
                        starAbun = candidate5Abun

                    else:
                        mirSeq = candidate5Seq
                        mirPos = candidate5Pos
                        mirAbun = candidate5Abun
                        starSeq = candidate3Seq
                        starPos = candidate3Pos
                        starAbun = candidate3Abun

                    f.write("%s-5p\tSequence: %s\tPosition:%s\t"\
                        "Abundance:%s\n" % (mirName, mirSeq, mirPos,
                        mirAbun))

                    f.write("%s-3p\tSequence: %s\tPosition:%s\t"\
                        "Abundance:%s\n" % (mirName, starSeq, starPos,
                        starAbun))

                    f.write("Match:%s, Mismatch:%s, Wobble:%s, Gap:%s, "\
                        "1-nt Variant Abundance:%s, Total sRNA Precursor "\
                        "Abundance:%s, Proportion:%s\n\n" % (matchCount,
                        mismatchCount, wobbleCount, gapCount, variantAbun,
                        totalAbun, proportion))

                    g.write(">%s\n%s\n" % (mirName, mirSeq))

def main():
    #findSequenceInIR('ATATTTCAGTTGCCTTTCTA', 'AAGGA-ACTAAAGT--AAGA-TATA-T-T-TCAGT-T-GC-CT-T-TCTATATGTAT-A-TATCAAAAGAAGGCTAA-GATCCCAAAAG-T-ATAAAGGAGATTTAAAAG', 17, 37, 20)
    #findSequenceInIR('TCCAACACGTTCTTCATCTTC', 'TATGGATTATTTATTG-TAATATCTTCTTCTGAACATATCGGAGTT-ATTGGA-GTC-CAACACGTTCTTCATCTTCTTTTCGGCCAGA-ACAT', 52, 73, 21)
    #findSequenceInIR('TGCCCGTAATTATAAAACATA', 'TAATTATAAAACATATTTAATTTAACATATGAAAC--A-ATACAAA-CAATACAAACAGTTT-TA', 35, 14, 21)
    #sys.exit()
    progStart = time.time()
    LibList = []

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
    # of the merged final file so that we can create the IRList using the
    # previously merged file
    else:
        IRFastaFilenamesList = [GenomeClass.IRFastaFilename]
        IRAlignmentFilenamesList = [GenomeClass.IRAlignmentFilename]

    # Combine the inverted repeat temp files of the einverted runs
    # into one file
    GenomeClass.combineIRTempFiles(IRFastaFilenamesList,
        IRAlignmentFilenamesList)

    # If the debug statement is set on, count the number of miRNAs from
    # the mirna input file that have an identified inverted repeat at its
    # position
    mirnasList = []
    if(debug): 
        mirnasList = findMirnasInIRs(GenomeClass.IRList, GenomeClass.chrDict,
            "Arabidopsis_mirnas.txt")

    #########################################################################

    ######################## Map small RNAs to genome #######################
 
    #########################################################################

    #I want to remove this loop so that we do not store multiple library files at once. 
    #This should help reduce memory footprint
    #Also, why am I getting 75% of precursors all of a sudden?

    precursorsWithDuplexDictByLib = {}
    finalCandidatesDictByLib = {}
    candidatesByLibDict = {}

    for libraryFilename in libFilenames:
        print("Beginning to process %s." % libraryFilename)
        Lib = Library(libraryFilename, GenomeClass.chrDict)

        precursorsWithDuplexDictByLib[libraryFilename] = {}
        finalCandidatesDictByLib[libraryFilename] = {}

        for chrName in sorted(GenomeClass.chrDict.keys()):
            chrIndex = GenomeClass.chrDict[chrName]
            
            precursorsWithDuplexDictByLib[libraryFilename][chrIndex] = {}
            finalCandidatesDictByLib[libraryFilename][chrIndex] = {}


        # Map small RNAs to the genome

        # If the parameter mapFilenames is not set, map the sRNA libraries
        # to the genome 
        if(not mapFilenames):
            mapFilename, logFilename = Lib.mapper(
                GenomeClass.indexFilename)

        else:
            print("Mapping will not be performed for file %s" % Lib.filename)

        #######################################################################

        ################# Map small RNAs to inverted repeats ##################

        #######################################################################

        # Create a dictionary with the sequence of all tags that
        # map to a position on every chromosome
        print("Creating the mapped list for %s" % Lib.filename)
        Lib.createMappedList(GenomeClass.chrDict)

        precursorsListByChr = []

        print("Mapping small RNAs to inverted repeats")

        if(parallel):
            # Run mapSRNAsToIRs in parallel
            res = pool.starmap_async(mapSRNAsToIRs, zip(GenomeClass.IRList,
                Lib.mappedList, repeat(Lib.libDict), repeat(mirnasList)))

            precursorsListByChr = res.get()

        else:
            for i in range(len(GenomeClass.chrDict)):
                #funcStart = time.time()

                precursorsListByChr.append(mapSRNAsToIRs(GenomeClass.IRList[i],
                    Lib.mappedList[i], Lib.libDict, mirnasList))

                #funcEnd = time.time()
                #execTime = round(funcEnd - funcStart, 2)
                #print("Time to map to chromosome: %s seconds" % (execTime))

        unfilteredFilename = os.path.splitext(Lib.filename)[0] +\
            '_all_precursors.txt'

        print("Writing precursors to a file")
    
        writePrecursors(unfilteredFilename, GenomeClass.chrDict,
            precursorsListByChr)

        # If the debug flag is set, we will want to count the number of
        # miRNAs in the miRNA input file that were sequenced in this
        # library
        if(debug):
            mirnasPresent = 0
            for mirna in mirnasList:
                if(mirna in Lib.libDict):
                    #print(mirna, Lib.libDict[mirna])
                    mirnasPresent += 1

            print("%s miRNAs found in %s." % (mirnasPresent, Lib.filename))

            # Find the number of miRNAs that have a reportable abundance
            # in the unfiltered precursor list (precursors that have at
            # least one small RNA that map to the 5' and 3' arm of the IR
            findMirnasInAllPrecursors(precursorsListByChr, mirnasList)

        #######################################################################

        ################### Filter precursor candidates #######################

        #######################################################################

        print("Filtering candidate precursors")

        # If we are running in parallel, run filterPrecursors in parallel.
        # Parallelize by chromosomes
        if(parallel):
            res = pool.starmap_async(filterPrecursors, zip(precursorsListByChr,
                repeat(overhang)))

            results = res.get()

            for chrName in sorted(GenomeClass.chrDict.keys()):
                chrIndex = GenomeClass.chrDict[chrName]
                # Get the index of chrDict[chrName]
                precursorsWithDuplexDictByLib[libraryFilename][chrIndex] \
                     = results[chrIndex][0]
                finalCandidatesDictByLib[libraryFilename][chrIndex] = \
                    results[chrIndex][1]

        # If running sequentially, run filterPrecursors sequentially
        else:
            for chrName in sorted(GenomeClass.chrDict.keys()):
                # Get the index of each chromosome that will be processed
                # sequentially
                chrIndex = GenomeClass.chrDict[chrName]
                precursorList = precursorsListByChr[chrIndex]

                precursorsWithDuplexDictByLib[libraryFilename][chrIndex], \
                    finalCandidatesDictByLib[libraryFilename][chrIndex] = \
                    filterPrecursors(precursorList, overhang)

        # Prior to writing this library's results, add its miRNAs and
        # corresponding precursors to a dictionary because we need to
        # identify a miRNA in more than 1 library to validate it
        for chrName, libraryFinalCandidatesDict in \
                finalCandidatesDictByLib[libraryFilename].items():
            for precursorName, duplexList in libraryFinalCandidatesDict.items():
                for duplex in duplexList:
                    if(precursorName in candidatesByLibDict):
                        candidatesByLibDict[precursorName].append(
                            libraryFilename)
                    else:
                        candidatesByLibDict[precursorName] = [libraryFilename]

        # Create filenames for the candidate files
        unfilteredPrecursorFilename = os.path.splitext(Lib.filename)[0] +\
            '_nearly_final_candidates'
        finalCandidateFilename = os.path.splitext(Lib.filename)[0] +\
            '_final_candidates'

        # Write the candidate miRNAs from this library to a file
        writeCandidates(unfilteredPrecursorFilename, GenomeClass.chrDict,
            precursorsWithDuplexDictByLib[libraryFilename])
        writeCandidates(finalCandidateFilename, GenomeClass.chrDict,
            finalCandidatesDictByLib[libraryFilename])

        progEnd = time.time()
        execTime = round(progEnd - progStart, 2)

    for precursorName, precursorList in candidatesByLibDict.items():
        if(len(precursorList) > 1):
            print(precursorName, precursorList)

#    if(debug): 
#        mirnasList = findMirnasInAllPrecursors(precursorsWithDuplexList, GenomeClass.chrDict,
#            "Arabidopsis_mirnas.txt")
#        mirnasList = findMirnasInAllPrecursors(finalCandidatesList, GenomeClass.chrDict,
#            "Arabidopsis_mirnas.txt")

    print("Total runtime was %s seconds" % execTime)

def findMirnasInIRs(IRList, chrDict, mirnaFilename):
    """
    This function counts the number of validated miRNAs that are
    found within the coordinates of an inverted repeat

    """

    countFound = 0
    countMirnas = 0
    countIRs = 0
    mirnasWithPrecursors = []
    countPrecursors = 0

    with open(mirnaFilename, 'r') as mirnaFile:
        for line in mirnaFile:
            mirnaFoundFlag = 0
            countMirnas += 1
            parsedLine = line.strip().split('\t')
            mirName = parsedLine[0]
            mirSeq = parsedLine[1]
            mirChr = parsedLine[2]
            mirPos = int(parsedLine[3])

            index = chrDict[mirChr]
            for repeat in IRList[index][::2]:
                countIRs += 1
                start5 = repeat[0]
                end5 = repeat[1]
                start3 = repeat[2]
                end3 = repeat[3]

                if(mirPos >= start5 and mirPos + (len(mirSeq)) <= end5 or
                        mirPos >= start3 and mirPos <= end3 + len(mirSeq)):
                    mirnasWithPrecursors.append(mirSeq)

                    if(not mirnaFoundFlag):
                        countFound += 1
                        mirnaFoundFlag = 1
                    countPrecursors += 1


    print("Found %s total precursors" % countPrecursors)
    print("Found %s out of %s. %s%%" % (countFound, countMirnas, countFound/countMirnas*100))

    return(mirnasWithPrecursors)

def findMirnasInAllPrecursors(precursorsListByChr, mirnasList):
    """
    Find the number of miRNAs that are present in the all_precursor file

    """

    countMirnasFound = 0
    mappedList = []
    flag = 0
    for precursorsList in precursorsListByChr:
        for precursor in precursorsList:
            mapped5Tags = precursor[2]
            mapped3Tags = precursor[3]

            for pos in mapped5Tags.keys():
                for mappedTag in mapped5Tags[pos]:
                    if(mappedTag[0] in mirnasList):
                        mappedList.append(mappedTag[0])
                        countMirnasFound += 1

            for pos in mapped3Tags.keys():
                for mappedTag in mapped3Tags[pos]:
                    if(mappedTag[0] in mirnasList):
                        mappedList.append(mappedTag[0])
                        countMirnasFound += 1

    print("%s miRNAs found in precursors file" % countMirnasFound)

    distinctMappedList = list(set(mappedList))
    print("%s distinct miRNAs found in precursors file" % \
        len(distinctMappedList))

if __name__ == '__main__':
    if parallel == 1:
        accel = int(multiprocessing.cpu_count()*0.50)
    else:
        pass

    main()
