import os
import subprocess
import time
import sys

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
        self.fastaFilename = "%s.fa" % os.path.splitext(self.filename)[0]

        if(not os.path.isdir("libs")):
            os.mkdir("libs")

        self.libDict, self.sumReads = self.readTagCount()
        self.normalizeCounts()

        # Create the mapped filename for this library
        self.mapFilename = "%s.map" % os.path.splitext(self.filename)[0]

        # Using chrDict from the genome file, create a tuple of multiple
        # dictionaries. Each dictionary in the tuple will represent a 
        # chromosome. Within these dictionaries are keys 'w' and 'c' for
        # each strand of the chromosome. These will then contain empty
        # ditionaries as values which will be populated by positions and
        # the sequences of the reads that map there.
        self.mappedList = [{'w': {}, 'c': {}} for i in range(len(chrDict))]

    def normalizeCounts(self):
        """Normalize the read counts using RPM and replace the counts in
        libDict with the normalized RPM

        """

        # Iterate through each tag sequenced in libDict
        for tag, countHitsList in self.libDict.items():
            count, hits = countHitsList

            # Determine the RPM of each read and replace the count in
            # libDict with its normalized RPM
            rpm = count*1000000/(self.sumReads)

            self.libDict[tag] = [rpm, 0]

    def readTagCount(self):
        """Read a library in tag count format into memory. Because we are
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
        sumReads = 0
        readCount = 1

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
                count = int(line.split('\t')[1].strip())
                libDict[tag] = [count, 0]
                sumReads += count

                g.write(">s_%s\n%s\n" % (readCount, tag))
                readCount += 1

            print("Total entries in library %s: %s" % (self.filename, readCount))
            g.close()

        # Stop timer for function
        funcEnd = time.time()

        # Calculate the execution time and print it to the user
        execTime = round(funcEnd - funcStart, 2)
        print("Time to read library %s: %s seconds" % (self.filename,
            execTime))

        return(libDict, sumReads)

    def mapper(self, indexFilename, bowtiePath, nthreads):
        """Map small RNAs to the provided index file

        Args:
            indexFilename: Path and name of the index for the genome. 
            bowtiePath: The path of bowtie
        Returns:
            Filename of mapped data

        """

        # Strip the filename of its folders and create the output map
        # name with that stripped filename in the libs folder
        indexNameStripped = os.path.basename(indexFilename)

        logFilename = "%s_bowtie.log" % os.path.splitext(self.mapFilename)[:-1]

        # Run bowtie
        print("Mapping small RNAs to the genome files for %s" %\
            (self.fastaFilename))

        with open(logFilename, 'w') as logFile:
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
                "-v 0", "-S", self.mapFilename, "--sam-nohead", "--no-unal"],
                stderr = logFile)

            if(returnCode):
                print("Something went wrong when running bowtie. Command was"\
                    "\n%s %s -f %s -a -m 50 --best --strata -v 0 -S "\
                    "%s -P %s --sam-nohead --no-unal" % (bowtiePath, indexFilename,
                    self.fastaFilename, self.mapFilename, nthreads))
                sys.exit()

        logFile.close()

        return(logFilename)

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
