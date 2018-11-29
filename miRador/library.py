import os
import subprocess
import time
import sys

class Library:
    """
    Class grouping of library functions and data structures
    """

    def __init__(self, filename, chrDict):
        self.filename = filename
        self.fastaFilename = "".join(self.filename.split('.')[:-1]) + '.fa'

        if(not os.path.isdir("libs")):
            os.mkdir("libs")

        self.libDict = self.readTagCount()

        # Create the mapped filename for this library
        self.mapFilename = "libs/%s.map" % ("".join(os.path.splitext(
            os.path.basename(self.filename))[0]))

        # Using chrDict from the genome file, create a tuple of multiple
        # dictionaries. Each dictionary in the tuple will represent a 
        # chromosome. Within these dictionaries are keys 'w' and 'c' for
        # each strand of the chromosome. These will then contain empty
        # ditionaries as values which will be populated by positions and
        # the sequences of the reads that map there.
        self.mappedList = [{'w': {}, 'c': {}} for i in range(len(chrDict))]

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

    def mapper(self, indexFilename, bowtiePath):
        """Map small RNAs to the provided index file

        Args:
            indexFilename: Path and name of the index for the genome. 
                This will be a fragment if fragFasta was run
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
                    "%s --sam-nohead --no-unal" % (bowtiePath, indexFilename,
                    self.fastaFilename, self.mapFilename))
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
