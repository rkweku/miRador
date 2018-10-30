import os
import re
import subprocess

class Genome:
    """
    Class grouping of genome functions and data structures
    """

    def __init__(self, filename, bowtieBuildPath):
        self.filename = filename
        self.chrDict = {}
        self.IRFastaFilename = "invertedRepeats/%s_Inverted_Seqs.fa" %\
            os.path.splitext(os.path.basename(self.filename))[0]
        self.IRAlignmentFilename = "invertedRepeats/%s_Inverted_Aligns.inv" %\
            os.path.splitext(os.path.basename(self.filename))[0]
        
        self.genomeSeqList = self.readFasta(self.filename)

        # Create an empty dictionary in IRDictByChr for as many chromosomes
        # exist
        self.IRDictByChr = [{} for i in range(len(self.chrDict))]

        # Build the bowtie index if it does not exist yet
        self.indexFilename = self.buildBowtieIndex(self.filename,
            bowtieBuildPath)

    def readFasta(self, filename):
        """Read a genome FASTA file into memory as a list

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

    def buildBowtieIndex(self, filename, bowtieBuildPath):
        """Code to create a bowtie index for the inverited repeats file.

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
        """Check if the index filename already exists. Return true if it
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

    def runEinverted(self, chrAndSeq, match, mismatch, gap, threshold,
            maxRepLen):
        """Fuunction to run einverted for a sequence
        
        Args:
            chrAndSeq: A tuple where the firstt element of the tuple is 
                the chr (or sequence name) and the second is the sequence
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

        # Separate the name and sequence from the tuple
        name = chrAndSeq[0]
        seq = chrAndSeq[1]

        # Each sequence needs to be in its own fasta file, so we will
        # first create temp input files for each sequence prior to calling
        # einverted
        tempInput = "invertedRepeats/%s.fa" % name
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

        # Delete the temporary FASTA file that was created so that
        # we could run einverted for this chromosome
        if os.path.exists(tempInput):
            os.remove(tempInput)

        # Close FNULL
        FNULL.close()

        return(outputFastaFilename, outputAlignmentFilename)

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

                        # We cannot have a precursor that has secondary
                        # stems or large loops within. Our criteria
                        # calls for no larger than 6 nucleotides in a row,
                        # so skip these precursors if we see them
                        #if("------" in hairpin5):
                        #    break

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
                        #if("------" in hairpin3):
                        #    break

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
            for file in garbage:
                os.remove(file)
