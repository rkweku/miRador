import os
import shutil
import subprocess
import sys

import log

from Bio import SeqIO

def writeFilteredPrecursors(filename, chrDict, IRDictByChr,
        precursorsDictByChr):
    """Write the precursors with their duplex to the output file

    Args:
        filename: The name of the output file
        chrDict: Dictionary of chromosomes and their corresponding
            positions
        IRDictByChr: List of dictionaries with the inverted repeat information
        precursorsDictByChr: Dictionary containing candidate precursors and
            their miRNA:miRNA* duplex for each chromosome

    """

    # Open the output file
    with open(filename, "w") as f, open("%s.fa" % \
            os.path.splitext(filename)[0], "w") as f_fasta:
        # loop through the chromosomes, sorted in numerical order
        for chrName in sorted(chrDict.keys()):
            chrIndex = chrDict[chrName]
            # Only proceed with this choromosome if there are precursors
            # in this chromosome
            try:
                precursorsDict = precursorsDictByChr[chrName]
            except:
                continue

            # Within each chromosome, loop through each precursor
            # and extract information to write to the file
            for precursorName, duplexDict in precursorsDict.items():
                coordinates = IRDictByChr[chrIndex][precursorName]
                duplexCount = 0

                # Write the precursor (its name, coordinates, and alignment)
                # to the precursor file
                f.write("%s,%s," % (precursorName, chrName))
                f.write(">%s\n" % precursorName)
                for i in range(len(coordinates)):
                    if(i == 5 or i == 6 or i == 7 or i == 8):
                        # We want to use U instead of T, but since the strand
                        # and alignment indicators (which are also written
                        # with this line) can't contain T, this should be
                        # safe here
                        f.write("%s\n" % coordinates[i].replace("T", "U"))

                    else:
                        f.write("%s," % coordinates[i]) 

                # Loop through all duplexes that were identified for a
                # precursor and write the stored information to the file
                for mirSeq, duplex in duplexDict.items():
                    arm = duplex[0]
                    starSeq = duplex[1]
                    mirPos = duplex[2]
                    starPos = duplex[3]
                    mirAbun = duplex[4]
                    mirHits = duplex[5]
                    starAbun = duplex[6]
                    starHits = duplex[7]
                    matchCount = duplex[8]
                    mismatchCount = duplex[9]
                    wobbleCount = duplex[10]
                    gapCount = duplex[11]
                    variantMirStrandAbun = duplex[12]
                    variantStarStrandAbun = duplex[13]
                    totalAbunMirStrand = duplex[14]
                    totalAbunStarStrand = duplex[15]
                    loopAbun = duplex[16]
                    variantAbun = variantMirStrandAbun + variantStarStrandAbun
                    totalAbun = totalAbunMirStrand + totalAbunStarStrand +\
                        loopAbun
                    proportion = duplex[17]

                    # In the event that there is more than one duplex for
                    # this precursor, we need a unique identifier for the
                    # miRNA and miRNA*
                    if(len(duplexDict.keys()) > 1):
                        mirName = "candidate-%s_%s" % (precursorName.split(
                            "precursor-")[1], duplexCount)
                        duplexCount += 1
                    else:
                        mirName = "candidate-%s" % precursorName.split(
                            "precursor-")[1]

                    f.write("%s-%s\tSequence: %s\tPosition:%s\t"\
                        "Abundance:%s\n" % (mirName, arm,
                        mirSeq.replace("T", "U"), mirPos,
                        mirAbun))

                    f.write("%s*\tSequence: %s\tPosition:%s\t"\
                        "Abundance:%s\n" % (mirName, starSeq.replace("T", "U"),
                        starPos, starAbun))

                    f.write("Match:%s, Mismatch:%s, Wobble:%s, Gap:%s, "\
                        "1-nt Variant Abundance:%s, Total sRNA Precursor "\
                        "Abundance:%s, Proportion:%s\n\n" % (matchCount,
                        mismatchCount, wobbleCount, gapCount, variantAbun,
                        totalAbun, proportion))

def findSequenceInIR(sequence, IRArm, tagLength):
    """If a sequence cannot be found in a simple search in an inverted
    repeat arm, it is because there are gaps that interrupt the sequence.
    Thus, we must do some work to find the tag, WITH gaps, as well as the
    start and end positions of the tag with gaps to identify potential
    duplexes

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
    localStart = IRArm.replace("-","").find(sequence)

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
            offset].count("-")

        # If the new offset is the same as the
        # old (number of gaps), set the
        # offsetUpdate flag to False to break
        if(offset == oldOffset):
            offsetUpdate = False

    # Initialize non gap base counts
    baseCount = 0

    # Initialize a variable to hold the sequence with gaps
    sequenceWithGaps = ""

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
        if(base != "-"):
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
    """This function gets the alignment between two bases

    Args:
        base1: The first nucleotide to be compared
        base2: The second nucleotide to be compared
    Returns:
        Match, mismatch, or gap

    """

    # Block to check if the alignment is a gap
    if(base1 == "-" or base2 == "-"):
        return("gap")

    # Block to check for simple mismatches (ie anything that cannot
    # possibly be a G-U wobble
    if(base1 == "A" and base2 != "T"):
        return("mismatch")
    elif(base1 == "C" and base2 != "G"):
        return("mismatch")

    # Block to check for a G-C match or G-U wobble
    elif(base1 == "G"):
        if(base2 == "C"):
            return("match")
        elif(base2 == "T"):
            return("wobble")
        else:
            return("mismatch")

    # Block to check for U-A match or U-G wobble
    elif(base1 == "T"):
        if(base2 == "A"):
            return("match")
        elif(base2 == "G"):
            return("wobble")
        else:
            return("mismatch")

    # If no other block was a hit, then we have a match
    else:
        return("match")

def getAlignment(arm5, arm3, alignStart, alignEnd):
    """This function gets the alignment that was agenerated by einverted
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
    """This function will serve two purposes. First, we need to get the
    sum of abundance of all eight 1 nt variants of the candidate sequence.
    Second,this function will return 0 and thus eliminate a candidate
    sequence if it is not the maximum of all 1 nt variants

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

    if(strand == "w"):
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

def filterPrecursors(mappedTagsToPrecursors, IRDict, libDict, overhang):
    """This function will perform the sRNA mapping and abundance filters.
    It will first try to find a miRNA and miRNA* pair by identifying
    tags that map to opposite sides of the precursor. It will also create
    splits of the c and w strand if there are tags that map to both

    Args:
        mappedTagsToPrecursors: Dictionary of tag information mapping
            to the precursor, identified by the precursor name
        IRDict: Dictionary of the inverted repeats in one chromosome
        libDict: The entire library dictionary to be queried for abundances
        overhang: Specific length of overhang that a duplex must have
    Returns:
        Dictionary of all precursors and the miRNA:miRNA* duplexes within
        that pass all filters

    """

    # Initialize our logger
    logger = log.setupLogger("filterPrecursors")

    # Initialize a dictionary to store our final candidaties that pass
    # all filters for this library
    finalCandidates = {}

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
                if(strand == "w"):
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
                # possibilities. Because there can be gaps in the
                # alignment, so we must identify which case (if not
                # both) it is before proceeding
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
                if(oldSequence5 != sequence5.replace("-","")):
                    logger.error("findSequenceInIR messed up for %s. "\
                        "Contact Reza to debug" % oldSequence5)
                    logger.error(precursorName, oldSequence5, sequence5,
                        local5Start, local5End)
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
                        if(strand == "w"):
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
                        if(oldSequence3 != sequence3.replace("-","")):
                            logger.error("findSequenceInIR messed up for %s. "\
                                "Contact Reza to debug" % oldSequence3)
                            logger.error(precursorName, oldSequence3,
                                sequence3, local3Start, local3End)
                            sys.exit()

                        # If there is an overhang on either the sequence,
                        # we have a candidate duplex and will investigate
                        # it further
                        if((strand == "c" and (local3Start - local5Start == 
                                overhang) and (local3End - local5End == 
                                overhang)) or (strand == "w" and (local5End 
                                - local3End == overhang) and (local5Start - 
                                local3Start == overhang))):
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
                            
                            # Only proceed if the alignment meets our filter
                            # specifications
                            if(gapCount + mismatchCount + (wobbleCount * .5) 
                                   <= 5 and gapCount <= 3):
                                # Get the hits information for the 5' and 3'
                                # tags from libDict
                                hits5 = libDict[mapped5Tag[0]][1]
                                hits3 = libDict[mapped3Tag[0]][1]

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
                                proportion = (variant5Abun + variant3Abun) /\
                                    (totalAbun5 + totalAbun3 + loopAbun)

                                # The 5' mapping tag will be kept as a candidate
                                # miRNA if it has at least an abundance of 3 RPM
                                if(tag5Abun >= 3):
                                    duplex = ("5p", mapped3Tag[0],
                                        candidate5Pos, candidate3Pos, tag5Abun,
                                        hits5, tag3Abun, hits3, matchCount, 
                                        mismatchCount, wobbleCount, gapCount,
                                        variant5Abun, variant3Abun, totalAbun5,
                                        totalAbun3, loopAbun, proportion)

                                    # If the sum of the two tags in the
                                    # make up more than 75% of the read
                                    # abundance in the entire precursor,
                                    # add the duplex to the candidates
                                    # dictionary
                                    if(proportion >= .75):
                                        # Add the precursor name as a key to
                                        # finalCandidates if it does not
                                        # yet exist. The value will be a list
                                        # of duplexes found in the precursor,
                                        # but the first element will be the IR
                                        # coordinates
                                        if(precursorName not in
                                                finalCandidates):
                                            finalCandidates[precursorName] = \
                                                {}

                                        finalCandidates[precursorName][\
                                            mapped5Tag[0]] = duplex

                                # The 3' mapping tag will be kept as a
                                # candidate miRNA if it has an abundance
                                # of at least 3 RPM
                                if(tag3Abun >= 3):
                                    duplex = ("3p", mapped5Tag[0],
                                        candidate3Pos, candidate5Pos, tag3Abun,
                                        hits3, tag5Abun, hits5, matchCount,
                                        mismatchCount, wobbleCount, gapCount,
                                        variant3Abun, variant5Abun, totalAbun3,
                                        totalAbun5, loopAbun, proportion)

                                    # If the sum of the two tags in the
                                    # make up more than 75% of the read
                                    # abundance in the entire precursor,
                                    # add the duplex to the candidates
                                    # dictionary
                                    if(proportion >= .75):
                                        # Add the precursor name as a key to
                                        # finalCandidates if it does not
                                        # yet exist. The valu will be a list of 
                                        # duplexes found in the precursor, but
                                        # the first element will be the IR
                                        # coordinates
                                        if(precursorName not in
                                                finalCandidates):
                                            finalCandidates[precursorName] = \
                                                {}

                                        finalCandidates[precursorName][
                                            mapped3Tag[0]] = duplex

    log.closeLogger(logger)

    return(finalCandidates)

def writeCandidates(outputFolder, candidatesByLibDict, filteredPrecursorsDict,
        IRDictByChr, libFilenamesList, chrDict, genomeFilename):
    """Write the candidate miRNAs to the provided filename

    Args:
        outputFolder: Folder for where the data will be written to
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
        genomeFilename: Path to the genome file to access sequences

    """

    outputFilename = "%s/preAnnotatedCandidates.csv" % outputFolder
    fastaFilename = "%s/preAnnotatedCandidates.fa" % outputFolder
    precursorsFilename = "%s/precursors.fa" % outputFolder

    # Open the output files
    with open(outputFilename, "w") as f, open(fastaFilename, "w") as g, \
            open(precursorsFilename, "w") as h:
        # Write the column names
        f.write("miR Name,Chr,Strand,miR Position,miR Sequence,miR Hits,"\
            "miR Length,Star Position,Star Sequence,Star Hits,Star Length,")

        for libName in libFilenamesList:
            libNameNoFolders = os.path.splitext(os.path.basename(libName))[0]
            f.write("miR Abun in {0},Star Abun in {0},1-nt Variants miR Abun "\
                    "in {0},1-nt Variants Star Abun in {0},Total Precursor "\
                    "Abun in {0},Proportion of reads from miR:miR* "\
                    "in {0}".format(libNameNoFolders))

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

            # Loop through each precursor in the chromosome and write
            # it to the file
            for precursorName, candidatesDict in precursorDict.items():
                coordinates = IRDictByChr[chrIndex][precursorName]
                start5 = coordinates[0]
                end3 = coordinates[3]
                strand = coordinates[5]
                multiDuplex = False
                mirCount = 1

                # Loop through every duplex in the candidate dictionary
                # and write them to the file
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
                    mirHits = duplexInfo[5]
                    starHits = duplexInfo[7]

                    # If there is one miRNA candidate in the duplex,
                    # we don't need to create a unique ID other than
                    # the name of the precursor
                    if(len(candidatesDict) == 1):
                        mirName = "candidate-%s-%s" % (precursorName.split(
                            "precursor-")[1], arm)

                    # If there is more than one miRNA candidate in the
                    # duplex, append mirCount to the precursor number
                    # to create a unique identifier for each miRNA
                    # in the duplex
                    else:
                        mirName = "candidate-%s_%s-%s" % (precursorName.split(
                            "precursor-")[1], mirCount, arm)
                        mirCount += 1

                    # Write the candidate name and the sequence to the
                    # fasta file
                    g.write(">%s\n" % mirName)
                    g.write("%s\n" % mirSeq.replace("T", "U"))

                    h.write(">%s\n" % mirName)
                    # Use samtools to identify the sequence from the genome
                    # file 
                    queryPos = "%s:%s-%s" % (chrName, coordinates[0],
                        coordinates[3])
                    proc = subprocess.Popen(["samtools", "faidx",
                        genomeFilename, queryPos], stdout = subprocess.PIPE)

                    # Grab the output parse to just the sequence to be written
                    output = proc.stdout.read()
                    output = output.decode("utf-8")
                    sequence = output.split("\n", 1)[1].replace("\n", "")
                    h.write("%s\n" % sequence)

                    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s," % (mirName,
                        chrName, strand, mirPos, mirSeq.replace("T", "U"),
                        mirHits, len(mirSeq), starPos,
                        starSeq.replace("T", "U"), starHits, len(starSeq)))

                    for libName in libFilenamesList:
                        libNameNoFolders = os.path.splitext(os.path.basename(
                            libName))[0]
                        # If the sequence was predicted in the current
                        # library, fill its columns with the relevant
                        # data
                        if(libNameNoFolders in libList):
                            duplexInfo = filteredPrecursorsDict[\
                                libNameNoFolders][chrName][precursorName][\
                                mirSeq]
                            mirAbun = duplexInfo[4]
                            starAbun = duplexInfo[6]
                            variantMirAbun = duplexInfo[12]
                            variantStarAbun = duplexInfo[13]
                            totalAbun5 = duplexInfo[14]
                            totalAbun3 = duplexInfo[15]
                            loopAbun = duplexInfo[16]
                            variantAbun = variantMirAbun + variantStarAbun
                            totalAbun = totalAbun5 + totalAbun3 + loopAbun
                            proportion = duplexInfo[17]

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
