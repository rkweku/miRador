import os
import re
import subprocess

class Genome:
    """
    Class grouping of genome functions and data structures
    """

    def __init__(self, filename, bowtieBuildPath):
        self.filename = filename

        ### Initialize several output folders if they do not exist yet
        # Create a path for genome if it does not exist already
        if(not os.path.isdir("genome")):
            os.mkdir('genome')
        # Create a path for the inverted repeat if it does not exist already
        if(not os.path.isdir("invertedRepeats")):
            os.mkdir("invertedRepeats")
        # Create a path for the bowtieIndex folder if it does not exist
        if(not os.path.isdir("genome/bowtieIndex")):
            os.mkdir('genome/bowtieIndex')

        self.IRFastaFilename = "invertedRepeats/%s_Inverted_Seqs.fa" %\
            os.path.splitext(os.path.basename(self.filename))[0]
        self.IRAlignmentFilename = "invertedRepeats/%s_Inverted_Aligns.inv" %\
            os.path.splitext(os.path.basename(self.filename))[0]

        # Split the genome file into multiple files - one for each chromosome        
        self.chrFilenamesList, self.chrDict = self.splitGenomeFasta()
        # Create an empty dictionary in IRDictByChr for as many chromosomes
        # exist
        self.IRDictByChr = [{} for i in range(len(self.chrDict))]

        # Build the bowtie index if it does not exist yet
        self.indexFilename = self.buildBowtieIndex(bowtieBuildPath)

    def splitGenomeFasta(self):
        """Split the genome file into multiple files, where each chromosome
           is its own file

        Returns:
            List of filenames where each chromosome was written to and
            a dictionary giving index value for each chromosome as it will
            appear in a list

        """

        # Create an empty list to hold the chromosome filenames
        chrFilenamesList = []
        counter = 1
        emptyCount = 0

        # Initialize a dictionary to be used to link the chromosome
        # name with its positional locations in list data structures like
        # IRDictByChr
        chrDict = {}

        # Read the full file into fastaFile so that we can parse on '>'
        f = open(self.filename, "r")
        fastaFile = f.read()
        f.close()

        # Loop through the FASTA file and split by sequence name in case
        # there are some sequence names with no sequence on the next line
        for line in fastaFile.split(">")[1:]:
            chromoInfo = line.partition('\n')

            # We will call the sequence identifier of the FASTA file
            # the chrName. Technically, this does not have to be a
            # chromosome if the user's file does not contain full
            # chromosomes, but the variable name will at least indicate
            # that
            chrName = chromoInfo[0].split()[0]
            # Add seqID to chrDict for indexing of chromosomes. Use same
            # indexing as the temp sequence name, so subtract 1 as index
            # starts at 0 on lists
            chrDict[chrName] = counter - 1

            # Remove any trailing characters from the sequence
            sequence = chromoInfo[2].rstrip()

            # Create a file for the chromosome and write it to a file
            splitFilename = "genome/%s_%s.fa" % (os.path.splitext(
                os.path.basename(self.filename))[0], counter)
            f_out = open(splitFilename, 'w')
            f_out.write('>%s\n%s\n' % (chrName, sequence))
            f_out.close()

            chrFilenamesList.append(splitFilename)

            counter += 1

        return(chrFilenamesList, chrDict)

    def checkBowtieNeedsUpdate(self, indexFilename):
        """Check if the index filename already exists. Return true if it
        does not and needs to be updated, but return false if doesn't
        need an update.

        Args:
            indexFilename: Name of the index file that would be created
                by bowtie-build for the genome file

        Returns:
            True if no update is needed, false if an update is needed

        """

        pattern = re.compile("%s(.*?).ebwt" % os.path.basename(indexFilename))

        # Check if a bowtie index file exists, and if it does, return false,
        # otherwise, return true so that bowtie-build can be run
        for filepath in os.listdir('genome/bowtieIndex'):
            if(pattern.match(filepath)):
                return(False)

        return(True) 

    def buildBowtieIndex(self, bowtieBuildPath):
        """Code to create a bowtie index for the inverited repeats file.

        Args:
            bowtieBuildPath: Path of bowtie-build
        Returns:
            Path of bowtie index

        """

        # Set the index filename. Remove any file extension and folders
        # from the filename path to ensure the index file is written
        # to the correct folder that is hardcoded here
        filenameStripped = os.path.splitext(self.filename.split('/')[-1])[0]
        indexFilename = "genome/bowtieIndex/%s" % (filenameStripped)

        if(self.checkBowtieNeedsUpdate(indexFilename)):
            print("Building a bowtie index for %s" % (self.filename))
            with open("genome/bowtieIndex/%s_bowtiebuild.log" %\
                    filenameStripped, 'w') as logFile:
                returnCode = subprocess.call([bowtieBuildPath, self.filename,
                    indexFilename], stdout = logFile)

            if(returnCode):
                print("Something went wrong when running bowtie-build. "\
                    "Command was\n%s %s %s" % (bowtieBuildPath, self.filename,
                    indexFilename))
                sys.exit()

            logFile.close()

        return(indexFilename)

    def combineIRTempFiles(self, IRFastaFilenamesList,
            IRAlignmentFilenamesList, runEInvertedFlag):
        """This function combines the temporary einverted files into one
        file for final analysis. However, if the user has opted to not run 
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

                toWriteList = []
                ### *******************
                # I've realized that I can do some filtering here before
                # adding the IRs to the combined file.
                # One such filter is to prevent complex duplex with no
                # secondary stems or large internal loops.
                # I can include a processing step when counter % 5 is 3 
                # to search for more than 5 '-' in a row
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

                        # We cannot have a precursor that has secondary
                        # stems or large loops within. Our criteria
                        # calls for no larger than 6 nucleotides in a row,
                        # so skip these precursors if we see them
                        if("------" not in hairpin3 and "------" not in 
                                hairpin5):

                            # Get the index of the chromosome to add the
                            # inverted repeat to
                            index = self.chrDict[chrName]

                            # Add the inverted repeat to the appropriate
                            # list within IRDictByChr
                            IRName = "precursor-%s" % IRCounter
                            self.IRDictByChr[index][IRName] = (start5, end5,
                                start3, end3, loop, 'w', hairpin5,
                                alignmentIndicators, hairpin3)

                            IRCounter += 1
                            IRName = "precursor-%s" % IRCounter
                            self.IRDictByChr[index][IRName] = (start5, end5,
                                start3, end3, loop, 'c', hairpin5,
                                alignmentIndicators, hairpin3)

                            if(runEInvertedFlag):
                                for entry in toWriteList:
                                    align_out.write(entry)
                            toWriteList = []
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
            for toDelete in garbage:
                if os.path.exists(toDelete):
                    os.remove(toDelete)

        return(IRCounter)


def runEinverted(einvertedPath, chrFilename, match, mismatch, gap,
        threshold, maxRepLen):
    """Fuunction to run einverted for a single chromosome
    
    Args:
        chrFilename: The path to the individual chromosome that will
            be run through einverted
        match: Score to pass to einverted for a match
        mismatch: Penalty score to pass to einverted for a mismatch
        gap: Score to pass to einverted for a gap
        threshold: Minimum total score an inverted repeat must have
            for einverted to record it
        maxRepLen: Maximum length an inverted repeat can have

    Returns:
        The name of the output FASTA file that einverted created, and
        the name of the alignment output file that einverted created.

    """

    outputFastaFilenamesList = []
    outputAlignmentFilenamesLis = []

    # Open FNULL to suppress the output of einverted becuase we do not
    # really need to know it is running for each proc
    FNULL = open(os.devnull, 'w')

    # Names of temporary output files to store results prior to merging
    outputFastaFilename = "invertedRepeats/%s.fa.temp" % os.path.splitext(
        os.path.basename(chrFilename))[0]
    outputAlignmentFilename = "invertedRepeats/%s.alignment.temp" % \
        os.path.splitext(os.path.basename(chrFilename))[0]

    # Call einverted utilizing this current sequence with the user
    # defined arguments from the config file.
    returnCode = subprocess.call([einvertedPath, "-sequence", chrFilename,
        "-gap", str(gap), "-threshold", str(threshold), "-match",
        str(match), "-mismatch", str(mismatch), "-maxrepeat",
        str(maxRepLen), "-outfile", outputAlignmentFilename,
        "-outseq", outputFastaFilename], stdout=FNULL,
        stderr=subprocess.STDOUT)

    # If a return code of anything but 0 is returned, it means there
    # was a problem and it should be investigated. Temp files wiill
    # remain from the run to assist in the debugging process
    if(returnCode != 0):
        print("Something went wrong when running einverted. Command was\n"\
            "%s -sequence %s -gap %s -threshold %s -match %s -mismatch "\
            "%s -maxrepeat %s -outfil %s -outseq %s" % (einvertedPath,
            tempInput, gap, threshold, match, mismatch, maxRepLen,
            outputAlignmentFilename, outputFastaFilename))
        sys.exit()

    # Close FNULL
    FNULL.close()

    return(outputFastaFilename, outputAlignmentFilename)
