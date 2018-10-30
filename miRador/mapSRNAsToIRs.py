def findOverlappingSeqsAndAbuns(position, subMappedDict, IREnd, libDict,
        IRMappedTags):

    """This function will serve a few purposes. First, we need to track
    the total abundance of sRNAs that map to each arm of an inverted repeat.
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
    """Map small RNAs to the inverted repeats. This will first read in
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

def writeUnfilteredPrecursors(filename, chrDict, IRDictByChr,
        mappedTagsToPrecursors):
    """Write the predicted precursors to a file with the designated file
    name

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

            # Only proceed with this choromosome if there are precursors
            # in this chromosome
            try:
                mappedTagsDict = mappedTagsToPrecursors[chrIndex]
            except:
                continue

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

                # Write the name and the chromosome number
                # followed by the coordinates of the precursor
                for i in range(len(coordinates)):
                    if(i == 0):
                        f.write('%s,' % chrName)
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
