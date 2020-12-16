#!/usr/local/bin/python3.3

# miRador is a plant miRNA prediction pipeline.
# It functions by first identifying inverted repeats within the genome.
# Then, it maps all small RNAs to the identified inverted repeats to create
# a list of candidate miRNAs.

#### Written by Reza Hammond
#### rkweku@udel.edu 

### miRNA identification criteria: https://doi.org/10.1105/tpc.17.00851

######## IMPORT ################
import configparser
import datetime
import log
import os
import multiprocessing
import time
import shutil
import sys

from itertools import repeat

import annotateCandidates
import filterPrecursors
import genome
import housekeeping
import library
import log
import mapSRNAsToIRs
import setupMiRBase

def createSummary(classificationCountsList, outputFolder, runTime):
    """Create a summary file for this analysis to report how many miRNAs
    of each type was predicted

    Args:
        classificationCountsList: A list of the counts identified in each
            classification category
        outputFolder: Name of the folder where the results will be written to
        runTime: The total runtime of the program

    """

    summaryFilename = "%s/summary.txt" % outputFolder

    with open(summaryFilename, "w") as f_out:
        f_out.write("Total miRNAs predicted: %s\n" % sum(
            classificationCountsList))
        f_out.write("Total known: %s\n" % classificationCountsList[0])
        f_out.write("Total identical to known miRNAs: %s\n" %\
            classificationCountsList[1])
        f_out.write("Total similar to at least one miRNA in this "\
            "organism: %s\n" %classificationCountsList[2])
        f_out.write("Total similar only to miRNAs outside this "\
            "organism: %s\n" % classificationCountsList[3])
        f_out.write("Total novel: %s\n" % classificationCountsList[4])
        f_out.write("Total runtime was %s seconds" % runTime)

def miRador():
    """Parse configuration file and make necessary calls to the various
    helper functions to perform the entire miRNA prediction of the user
    provided input files. This function primarily serves as a wrapper
    to those other functions in other files

    """

    # Initialize our logger
    logger = log.setupLogger("miRador")

    progStart = time.time()

    ######################## Parse Config File ###############################
    configFilename = sys.argv[1]
    config = configparser.ConfigParser()
    config.read(configFilename)

    # Get the preprocessing arguments
    #runPreprocessFlag = config.get("Preprocess", "runPreprocessFlag")

    # Get the genome file name
    genomeFilename = config.get("Genome", "genomeFilename")

    # Get the einverted arguments
    runEInvertedFlag = config.getint("EInverted", "runEInvertedFlag",
        fallback = 1)
    einvertedPresets = (config.get("EInverted", "einvertedPresets",
        fallback = "medium"))

    # If einvertedPresets is set, set the einverted parameters to
    # appropriate levels for prediction of inverted repeats
    if(einvertedPresets.lower() == "medium"):
        match = 3
        mismatch = -4
        gap = 6
        threshold = 45
        maxRepLen = 300

    elif(einvertedPresets.lower() == "low"):
        match = 3
        mismatch = -4
        gap = 6
        threshold = 40
        maxRepLen = 300

    elif(einvertedPresets.lower() == "high"):
        match = 3
        mismatch = -5
        gap = 7
        threshold = 50
        maxRepLen = 300

    # Get the advanced einverted arguments from the config file which
    # will override the presets if these are set
    advancedMatch = config.get("Advanced", "match")
    advancedMismatch = config.get("Advanced", "mismatch")
    advancedGap = config.get("Advanced", "gap")
    advancedThreshold = config.get("Advanced", "threshold")
    advancedMaxRepLen = config.get("Advanced", "maxRepLen")

    # If the advanced settings are set, override whatever has been set in them
    if(advancedMatch):
        match = int(advancedMatch)
    if(advancedMismatch):
        mismatch = int(advancedMismatch)
    if(advancedGap):
        gap = int(advancedGap)
    if(advancedThreshold):
        threshold = int(advancedThreshold)
    if(advancedMaxRepLen):
        maxRepLen = int(advancedMaxRepLen)

    # Get the Libraries arguments and parse the libraries into a list
    # of strings. User can input a list of files or just a directory
    # holding all of the tag count files
    libFilenamesList = []
    libFilenamesString = config.get("Libraries", "libFilenamesList",
        fallback = "")
    libFolder = config.get("Libraries", "libFolder", fallback = "")

    # If individual libraries were given, split the string on commas and
    # store them in libFilenamesList 
    if(libFilenamesString):
        libFilenamesList = libFilenamesString.split(",")

    # If libFolder was specified, loop through the files in the folder
    # and add all files to libFilenamesList
    if(libFolder):
        for file in os.listdir(libFolder):
            libFilenamesList.append("%s/%s" % (libFolder,
                os.path.join(file)))

    # Do a check to confirm the user did not enter the same library
    # multiple times in libFilenamesList. First, we don't want to
    # process a library twice, but we also want to make sure the user
    # also did not intend to place another library in and accidentally
    # just pasted the path to another twice
    numLibs = len(libFilenamesList)
    if(numLibs != len(set(libFilenamesList))):
        logger.error("It appears that a library was input more than once. "
            "Please check your libraries again, remove any duplicate "\
            "entries, and ensure all libraries you intend to process are "
            "present.")
        sys.exit()

    # Grab the information for the BLAST variables
    organism = config.get("miRBase", "organism").lower()
    version = config.get("miRBase", "version", fallback = "CURRENT")

    cleanupFlag = config.getint("General", "cleanupFlag", fallback = 1)
    parallel = config.getint("General", "parallel")
    nthreads = config.get("General", "nthreads")
    blastnPath = os.path.expanduser(config.get("General", "blastnPath"))
    bowtiePath = os.path.expanduser(config.get("General", "bowtiePath"))
    bowtieBuildPath = os.path.expanduser(config.get("General",
        "bowtieBuildPath"))
    einvertedPath = os.path.expanduser(config.get("General", "einvertedPath"))
    makeblastdbPath = os.path.expanduser(config.get("General",
        "makeblastdbPath"))
    perlPath = os.path.expanduser(config.get("General", "perlPath"))
    RNAFoldPath = os.path.expanduser(config.get("General", "RNAFoldPath"))
    RNAPlotPath = os.path.expanduser(config.get("General", "RNAPlotPath"))
    ps2pdfwrPath = os.path.expanduser(config.get("General", "ps2pdfwrPath"))
    outputFolder = config.get("General", "outputFolder", fallback = "")

    # Required overhang between top and bottom strands of miRNA duplex
    # Hardcoded to 2 here, but in such a way that could technically allow
    # modifications
    overhang = 2

    # Perform various housekeeping functions including the checks that all
    # external program dependencies exist, that files being referenced and
    # folders that will be written to exist and are created.
    mirBaseDict, outputFolder = housekeeping.housekeeping(genomeFilename,
        libFilenamesString, libFolder, libFilenamesList, bowtiePath,
        bowtieBuildPath, einvertedPath, blastnPath, makeblastdbPath,
        perlPath, RNAFoldPath, RNAPlotPath, ps2pdfwrPath, outputFolder,
        organism, version)

    # Set the number of cores, if parallel is on
    if(parallel):
        nproc = int(round(int(multiprocessing.cpu_count()*.5),1))

    # Create genome object
    GenomeClass = genome.Genome(genomeFilename, bowtieBuildPath)

    ##########################################################################

    ############### Find inverted repeats in genome file #####################

    ##########################################################################
    # Run EInverted if the flag is set
    if(runEInvertedFlag):
        # Create an empty list for both the inverted repeat FASTA files
        # and alignment files
        IRFastaFilenamesList = []
        IRAlignmentFilenamesList = []

        # If parallel is set, run einverted using the parallel version
        if(parallel):
            logger.info("Running einverted in parallel")

            if(len(GenomeClass.chrFilenamesList) < nproc):
                pool = multiprocessing.Pool(len(GenomeClass.chrFilenamesList))
            else:
                pool = multiprocessing.Pool(nproc)

            res = pool.starmap_async(genome.runEinverted,
                zip(repeat(einvertedPath), GenomeClass.chrFilenamesList,
                repeat(match), repeat(mismatch), repeat(gap),
                repeat(threshold), repeat(maxRepLen)))

            results = res.get()

            pool.close()

            # Loop through the results and add the inverted repeat filenames
            # to their respective lists
            for result in results:
                IRFastaFilenamesList.append(result[0])
                IRAlignmentFilenamesList.append(result[1])

        else:
            logger.info("Running einverted sequentially")

            # Loop through each chromosome and run einverted on each, one at 
            # a time
            for chrFilename in GenomeClass.chrFilenamesList:
                IRName, IRSeq = genome.runEinverted(einvertedPath, 
                    chrFilename, match, mismatch, gap, threshold, maxRepLen)

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
        IRAlignmentFilenamesList, runEInvertedFlag)

    #########################################################################

    ######################## Map small RNAs to genome #######################
 
    #########################################################################

    filteredPrecursorsDict = {}
    candidatesByLibDict = {}

    # Populate candidatesByLibDict chromosomes with empty dictionaries
    for chrName, chrIndex in GenomeClass.chrDict.items():
        candidatesByLibDict[chrName] = {}

    # Initialize libCounter to help inform users how far along the run is
    libCounter = 1

    for libraryFilename in libFilenamesList:
        libNameNoFolders = os.path.splitext(os.path.basename(
            libraryFilename))[0]

        logger.info("Beginning to process %s, library %s of %s." % (
            libraryFilename, libCounter, len(libFilenamesList)))
        Lib = library.Library(libraryFilename, GenomeClass.chrDict)

        filteredPrecursorsDict[libNameNoFolders] = {}

        for chrName in sorted(GenomeClass.chrDict.keys()):
            filteredPrecursorsDict[libNameNoFolders][chrName] = {}

        # Map small RNAs to the genome
        logger.info("Running bowtie on %s" % Lib.filename)
        funcStart = time.time()

        logFilename = Lib.mapper(GenomeClass.indexFilename, bowtiePath,
            nthreads)

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Runtime of bowtie for %s: %s seconds" % \
            (Lib.mapFilename, execTime))

        logger.info("Creating the mapped list for %s" % Lib.filename)
        funcStart = time.time()

        # Create a dictionary with the sequence of all tags that
        # map to a position on every chromosome
        Lib.createMappedList(GenomeClass.chrDict)

        # Normalize the reads in libDict
        Lib.normalizeReads(logFilename)

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Time to create the mappedList: %s seconds" % (execTime))
    
        #######################################################################

        ################# Map small RNAs to inverted repeats ##################

        #######################################################################

        logger.info("Mapping sRNAs to the inverted repeats")

        funcStart = time.time()

        mappedTagsToPrecursors = []

        # Parallelization of this module has been removed as the overhead of
        # transferring mappedList to each proc is quite significant while
        # its runtime on one proc is extremely quick
#        if(parallel):
#            # Run mapSRNAsToIRs in parallel
#            pool = multiprocessing.Pool(nproc)
#
#            res = pool.starmap_async(mapSRNAsToIRs.mapSRNAsToIRs,
#                zip(GenomeClass.IRDictByChr, Lib.mappedList,
#                repeat(Lib.libDict)))
#
#            mappedTagsToPrecursors = res.get()
#
#            pool.close()

        # Map the sRNAs for this library to the inverted repeats predicted
        # for this genome.
        # Note that the format of this where we run one chromosome at a
        # time is a holdover from the parallelization effort that was made
        # for this function. While this for loop can now be moved into the
        # function, I am keeping it outside quite simply because the tabbing
        # within this function became quite deep and I'd rather avoid going
        # another level deeper
        for i in range(len(GenomeClass.chrDict)):
            mappedTagsToPrecursors.append(mapSRNAsToIRs.mapSRNAsToIRs(
                GenomeClass.IRDictByChr[i], Lib.mappedList[i], Lib.libDict))

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Time to map sRNAs to inverted repeats: %s seconds" \
            % (execTime))

        logger.info("Writing precursors to a file")

        # Create a file for all precursors to be written to that have at least
        # one sRNA that maps to both strands
        unfilteredFilename = "%s/libs/%s_all_precursors.txt" % (outputFolder,
            libNameNoFolders)

        mapSRNAsToIRs.writeUnfilteredPrecursors(unfilteredFilename,
            GenomeClass.chrDict, GenomeClass.IRDictByChr,
            mappedTagsToPrecursors)

        #######################################################################

        ################### Filter precursor candidates #######################

        #######################################################################

        logger.info("Filtering candidate precursors")

        funcStart = time.time()

        # Parallelization of this module has been removed as the overhead of
        # transferring mappedList to each proc is quite significant while
        # its runtime on one proc is extremely quick
#        if(parallel):
#            pool = multiprocessing.Pool(nproc)
#
#            res = pool.starmap_async(filterPrecursors.filterPrecursors,
#                zip(mappedTagsToPrecursors, GenomeClass.IRDictByChr,
#                repeat(overhang)))
#
#            results = res.get()
#
#            pool.close()
#
#            for chrName in sorted(GenomeClass.chrDict.keys()):
#                chrIndex = GenomeClass.chrDict[chrName]
#
#                # Get the index of chrDict[chrName]
#                filteredPrecursorsDict[libNameNoFolders][chrName] = \
#                    results[chrIndex][1]

        # Filter the precursors, one chromosome at a time
        # Note that the format of this where we run one chromosome at a
        # time is a holdover from the parallelization effort that was made
        # for this function. While this for loop can now be moved into the
        # function, I am keeping it outside quite simply because the tabbing
        # within this function became quite deep and I'd rather avoid going
        # another level deeper
        for chrName in sorted(GenomeClass.chrDict.keys()):
            # Get the index of each chromosome that will be processed
            # sequentially
            chrIndex = GenomeClass.chrDict[chrName]
            precursorList = mappedTagsToPrecursors[chrIndex]
            IRDict = GenomeClass.IRDictByChr[chrIndex]

            filteredPrecursorsDict[libNameNoFolders][chrName] = \
                filterPrecursors.filterPrecursors(precursorList, IRDict,
                Lib.libDict, overhang)

        funcEnd = time.time()
        execTime = round(funcEnd - funcStart, 2)
        logger.info("Time to filter inverted repeats: %s seconds" % \
            (execTime))

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

    
        # Create a file for all precursors that have been identified as having
        # a valid miRNA:miRNA* duplex to be written to
        filteredFilename = "%s/libs/%s_candidate_precursors.txt" % (
            outputFolder, libNameNoFolders)

        funcStart = time.time()
        filterPrecursors.writeFilteredPrecursors(filteredFilename,
            GenomeClass.chrDict, GenomeClass.IRDictByChr,
            filteredPrecursorsDict[libNameNoFolders])

        # Increment the library counter
        libCounter += 1

    filterPrecursors.writeCandidates(outputFolder, candidatesByLibDict,
        filteredPrecursorsDict, GenomeClass.IRDictByChr, libFilenamesList,
        GenomeClass.chrDict)

    ##########################################################################

    ################### Annotate candidate miRNAs ############################

    ##########################################################################
    logger.info("Annotating candidate miRNAs")

    funcStart = time.time()

    subjectSequencesFilename = "miRBase/miRBasePlantMirnas.fa"
    queryMirnasFilename = "%s/preAnnotatedCandidates.fa" % outputFolder
    dbFilename = "miRBase/miRBasePlantMirnas.db"

    # Create a list of candidate miRNAs and mirBase miRNAs
    querySeqsList = annotateCandidates.createListForAlign(queryMirnasFilename)
    subjectSeqsList = annotateCandidates.createListForAlign(
        subjectSequencesFilename)

    similarityDict = {}

    if(parallel):
        logger.info("Running sequence alignment in parallel")

        pool = multiprocessing.Pool(nproc)

        res = pool.starmap_async(annotateCandidates.pairwiseAlignmentParallel,
            zip(querySeqsList, repeat(subjectSeqsList), repeat(organism)))

        results = res.get()

        pool.close()

        for result in results:
            similarityDict.update(result)

    else:
        similarityDict = annotateCandidates.pairwiseAlignment(querySeqsList,
            subjectSeqsList, organism)

    # Properly annotate the candidate miRNAs with the data in similarityDict
    classificationCountsList = annotateCandidates.annotateCandidates(
        outputFolder, similarityDict, organism, mirBaseDict,
        GenomeClass.IRDictByChr, numLibs, GenomeClass.chrDict,
        GenomeClass.chrFilenamesList, perlPath, RNAFoldPath, RNAPlotPath,
        ps2pdfwrPath)

    # Delete the single chromosome files used by einverted and the
    # draw functions to clean up temp file
    for chrFilename in GenomeClass.chrFilenamesList:
        os.remove(chrFilename)

    funcEnd = time.time()
    execTime = round(funcEnd - funcStart, 2)
    logger.info("Time to annotate candidate miRNAs: %s seconds" % (execTime))

    progEnd = time.time()
    execTime = round(progEnd - progStart, 2)

    logger.info("Total runtime was %s seconds" % execTime)

    log.closeLogger(logger)

    # Write a summary file with details of the analysis
    createSummary(classificationCountsList, outputFolder, execTime)

    # If cleanup flag is set, remove the temp folder containing bowtie
    # output data as well as potential temp fasta files created when 
    # a user input tag count files
#    if(cleanupFlag):
#        shutil.rmtree("miRadorTempFolder")
