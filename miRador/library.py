import os
import re
import subprocess
import time
import sys

import log

## This function is saved in the event that we ever decide to perform bowtie
## once for all input library sequences. 
#def mergeLibSeqs(LibList):
#    """
#    Create a list of unique tags for the purpose of mapping unique
#    sequences to the genome
#
#    Args:
#        libDictList: A list of dictionaries of each library
#
#    Returns:
#        List of unique sequences
#
#    """
#
#    # Create an empty list for the unique sequences
#    allSeqs = []
#
#    # Loop through each Library object and add all sequences to uniqueSeqs
#    for LibClass in LibList:
#        allSeqs.extend(list(LibClass.libDict.keys()))
#
#    # Create a set of the allSeqs to get a list of unique sequences so
#    # that they may be mapped with bowtie
#    uniqueSeqs = list(set(allSeqs))
#
#    # Open a file to write the unique sequences to and create a counter
#    # to name the sequences
#    f = open("libs/uniqueSeqs.fa", "w")
#    counter = 0
#
#    # Write the unique sequences to a file so that they can be mapped
#    # to the inverted repeats together
#    for sequence in uniqueSeqs:
#        f.write(">s_%s\n%s\n" % (counter, sequence))
#        counter += 1
#
#    return(uniqueSeqs)

class Library:
    """
    Class grouping of library functions and data structures
    """

    def __init__(self, filename, chrDict):
        self.filename = filename
        self.mappedReads = 0

        # Create the name of the temp fasta file. This is needed for
        # the bowtie mapping step
        self.fastaFilename = "miRadorTempFolder/fastaLibs/%s.fa" %\
            os.path.basename(self.filename)
        # Create a folder for the temp fasta files if it doesn't
        # exist yet
        if(not os.path.isdir("miRadorTempFolder/fastaLibs")):
            os.mkdir("miRadorTempFolder/fastaLibs")

        # Identify the type of file format of this library
        self.libType = self.identifyFileType()

        # Identify the format of the input library file
        if(self.libType == "fasta"):
            self.libDict = self.readFasta()
        elif(self.libType == "fastq"):
            self.libDict = self.readFastq()
        elif(self.libType == "tagCount"):
            self.libDict = self.readTagCount()

        # Create the mapped filename for this library
        self.mapFilename = "miRadorTempFolder/bowtieOutput/%s.map" %\
             os.path.basename(self.filename)

        # Using chrDict from the genome file, create a tuple of multiple
        # dictionaries. Each dictionary in the tuple will represent a 
        # chromosome. Within these dictionaries are keys "w" and "c" for
        # each strand of the chromosome. These will then contain empty
        # ditionaries as values which will be populated by positions and
        # the sequences of the reads that map there.
        self.mappedList = [{"w": {}, "c": {}} for i in range(len(chrDict))]

    def identifyFileType(self):
        """Investigate the first few linse of the library filename to
        determine its format so that the proper file parser function
        is used

        Returns:
            A simple string of fasta, fastq, or tagCount

        """

        # Initialize our logger
        logger = log.setupLogger("identifyFileType")

        libType = ""

        with open(self.filename) as f:
            lines = [line for line in f][:4]

            # If the first character of the first line of the file is a >, and
            # both the 2nd and 4th lines are nucleotide sequences, then we 
            # should be able to call this file fasta
            if(lines[0][0] == ">" and re.search("^[ACGT]*$",
                    lines[1].rstrip().upper()) and re.search("^[ACGT]*$",
                    lines[3].rstrip().upper())):
                libType = "fasta"

            # If the second line of the file is a nucleotide sequence, but
            # the 4th is not, then it should be a fastq file
            elif(re.search("^[ACGTN]*$", lines[1].rstrip().upper()) and 
                    not re.search("^[ACGTN]*$", lines[3].rstrip().upper())):
                libType = "fastq"

            # If splitting the first line on a tab results in the first
            # index being just a nucleotide sequence, then the input file
            # is a tag count file
            elif(re.search("^[ACGTN]*$", lines[0].split("\t")[0].upper())):
                libType = "tagCount"

            # If none of the previous tests passed, report an error to
            # the user and kill the run
            else:
                logger.info("The data in %s was not recognized as a fasta, "\
                    "fastq, or tag count file. Please check the file to "\
                    "ensure that it is one of the recognized file types." %\
                    self.filename)
                log.closeLogger(logger)
                sys.exit()

        log.closeLogger(logger)

        return(libType)

    def readFasta(self):
        """Read a library in FASTA format into memory as a dictionary.
        The key will be the sequence, and the value will be the count of
        times in which the tag was seen

        Return:
            Dictionary of the library file

        """

        # Initialize our logger
        logger = log.setupLogger("readFasta")

        # Start timer for function
        funcStart = time.time()

        # Create an empty dictionary to store the full library
        libDict = {}

        # Initialize a counter to simply name the reads in the temp fasta file
        readCount = 1

        # Open the file and loop through line by line to store the library into
        # a dictionary. Each tag will be a key and the number of times seen
        # will be the value
        with open(self.filename) as f, open(self.fastaFilename, "w") as g:
            for count, line in enumerate(f, start=0):
                if(count % 2 == 1):
                    # Store the sequence and count into variables, then add
                    # them to the dictionary. There should not be any
                    # duplicate sequences in this format, however,
                    tag = line.rstrip()

                    # If the tag already exists in the libDict, increment
                    # the counter of the sequence and retain the empty hits 
                    # variable as 0
                    if(tag in libDict):
                        currCount = libDict[tag][0]
                        libDict[tag] = [currCount + 1, 0]
                    else:
                        libDict[tag] = [1, 0]

                        # Write the sequence to a unique reads FASTA file
                        g.write(">s_%s\n%s\n" % (readCount, tag))
                        readCount += 1

        # Stop timer for function
        funcEnd = time.time()

        # Calculate the execution time and print it to the user
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Time to read library %s: %s seconds" % (self.filename,
            execTime))

        log.closeLogger(logger)

        return(libDict)

    def readFastq(self):
        """Read a library in FASTQ format into memory as a dictionary.
        The key will be the sequence, and the value will be the count of
        times in which the tag was seen

        Return:
            Dictionary of the library file

        """

        # Initialize our logger
        logger = log.setupLogger("readFastq")

        # Start timer for function
        funcStart = time.time()

        # Create an empty dictionary to store the full library
        libDict = {}

        # Initialize a counter to simply name the reads in the temp fasta file
        readCount = 1

        # Open the file and loop through line by line to store the library into
        # a dictionary. Each tag will be a key and the number of times seen
        # will be the value
        with open(self.filename) as f, open(self.fastaFilename, "w") as g:
            for count, line in enumerate(f, start=0):
                if(count % 4 == 1):
                    # Store the sequence and count into variables, then add
                    # them to the dictionary. There should not be any
                    # duplicate sequences in this format, however,
                    tag = line.rstrip()

                    # If the tag already exists in the libDict, increment
                    # the counter of the sequence and retain the empty hits 
                    # variable as 0
                    if(tag in libDict):
                        currCount = libDict[tag][0]
                        libDict[tag] = [currCount + 1, 0]
                    else:
                        libDict[tag] = [1, 0]

                        # Write the sequence to a unique reads FASTA file
                        g.write(">s_%s\n%s\n" % (readCount, tag))
                        readCount += 1

        # Stop timer for function
        funcEnd = time.time()

        # Calculate the execution time and print it to the user
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Time to read library %s: %s seconds" % (self.filename,
            execTime))

        log.closeLogger(logger)

        return(libDict)

    def readTagCount(self):
        """Read a library in tag count format into memory. Additionally,
        as each sequence is read, we will write this to a FASTA output
        file for bowtie

        Return:
            Dictionary of library read in

        """

        # Initialize our logger
        logger = log.setupLogger("readTagCount")

        # Start timer for function
        funcStart = time.time()

        # Create an empty dictionary to store the full library
        libDict = {}

        # Initialize a counter to simply name the reads in the temp fasta file
        readCount = 1

        # Open the file and loop through line by line to store the library into
        # a dictionary. Each tag will be a key and the abundance will be the
        # value
        with open(self.filename) as f, open(self.fastaFilename, "w") as g:
            for line in f:
                # Store the sequence and count into variables, then add them
                # to the dictionary. There should not be any duplicate
                # sequences in this format, however, 
                tag = line.split("\t")[0]
                count = int(line.split("\t")[1].strip())
                libDict[tag] = [count, 0]

                # Write the sequence to a unique reads FASTA file
                g.write(">s_%s\n%s\n" % (readCount, tag))
                readCount += 1

        # Stop timer for function
        funcEnd = time.time()

        # Calculate the execution time and print it to the user
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Time to read library %s: %s seconds" % (self.filename,
            execTime))

        log.closeLogger(logger)

        return(libDict)

    def mapper(self, indexFilename, bowtiePath, nthreads):
        """Map small RNAs to the provided index file

        Args:
            indexFilename: Path and name of the index for the genome. 
            bowtiePath: The path of bowtie
            nthreads: The number of threads to use with bowtie
        Returns:
            Filename of mapped data

        """

        # Initialize our logger
        logger = log.setupLogger("mapper")

        # Strip the filename of its folders and create the output map
        # name with that stripped filename in the libs folder
        indexNameStripped = os.path.basename(indexFilename)

        logFilename = "%s_bowtie.log" % os.path.splitext(self.mapFilename)[:-1]

        if(self.libType == "tagCount"):
            logger.info("Mapping small RNAs to the genome files for %s" %\
                (self.fastaFilename))
        else:
            logger.info("Mapping small RNAs to the genome files for %s" %\
                (self.filename))

        with open(logFilename, "w") as logFile:
            # Run bowtie with the following options:
            # -a to report all valid alignments as we want multihits
            # -m 50 to suppress all alignments with more than 50 matches
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
            returnCode = subprocess.call([bowtiePath, indexFilename, "-f",
                self.fastaFilename, "-a", "-m 50", "--best", "--strata",
                "-v 0", "-S", self.mapFilename, -P, nthreads, 
                "--sam-nohead", "--no-unal"], stderr = logFile)

            # If there is a return code, report an error to the user and exit
            if(returnCode):
                logger.error("Something went wrong when running bowtie. "\
                    "Command was\n%s %s -f %s -a -m 50 --best --strata "\
                    "-v 0 -S %s -P %s --sam-nohead --no-unal" %\
                    (bowtiePath, indexFilename, self.fastaFilename, 
                    self.mapFilename, nthreads))

                sys.exit()

        logFile.close()

        log.closeLogger(logger)

        return(logFilename)

    def normalizeReads(self, logFilename):
        """Normalize all of the reads in libDict

        Args:
            logFilename: bowtie log file that will be opened to view the
                total number of alignments

        """

        # Get tne number of genome matched reads from the bowtie log file
        logIn = open(logFilename, "r")
        logFile = logIn.readlines()
        logIn.close()

        toParseLine = logFile[-1]

        self.mappedReads = int(re.search("Reported (.*) alignments",
            toParseLine).group(1))

        for sequence, countsHits in self.libDict.items():
            counts = countsHits[0]
            hits = countsHits[1]

            if(hits):
                # Normalize the counts to RPM
                normalizedAbun = (counts * 1000000) / self.mappedReads
                # Replace the counts in libDict with hits normalized abundnace
                self.libDict[sequence][0] = normalizedAbun / hits

    def createMappedList(self, chrDict):
        """Create a dictionary to hold the mapped sRNAs for simple querying

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

        with open(self.mapFilename, "r") as inFile:
            for line in inFile:
                flag = line.split("\t")[1]
                chrName = line.split("\t")[2]
                position = int(line.split("\t")[3])
                sequence = line.split("\t")[9]

                # If the flag is a star, it means that it did not align
                # and thus must we break from the loop here
                if(chrName == "*"):
                    continue

                # If the 16th bit flag is not set, the strand is w
                elif(not int(flag) & 0b10000):
                    strand = "w"

                # If the 16th bit of the flag is set, the strand is c. We 
                # thus need to make the sequence the reverse complement to
                # ensure we can find the right sequence later
                elif(int(flag) & 0b10000):
                    strand = "c"
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
                    self.mappedList[chrIndex][strand][position].append(
                        sequence)

                # Increment the hits variable for this sequence in libDict
                # to indicate its number of mapped locations found
                self.libDict[sequence][1] += 1
