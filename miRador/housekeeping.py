import logging
import os
import shutil
import sys
import time

import log
import setupMiRBase

def housekeeping(genomeFilename, libFilenamesString, libFolder,
        libFilenamesList, bowtiePath, bowtieBuildPath, einvertedPath,
        perlPath, outputFolder, organism, version):
    """Perform various housekeeping functions including the checks that all
    external program dependencies exist, that files being referenced and
    folders that will be written to exist and are created. Additionally,
    it will also call setupMirBase.py to create download the current
    version of miRBase if needed to prepare for the annotation of our
    candidate miRNAs

    Args:
        genomeFilename: The path of the genome file
        libFilenamesString: The raw text of the library files that the
            user would have supplied in the config file
        libFolder: The folder of the library files if the user chose
            to supply that instead of individual library paths
        libFilenamesList: The list of library paths that have already
            been parsed either from libFilenamesString or supplied by
            finding files that end with chopped.txt from libFolder
        bowtiePath: The path of bowtie on the system
        einvertedPath: The path of einverted on the system
        perlPath: The path of perl on the system
        outputFolder: The config entry for the output folder. Can be blank
        organism: The three letter identifier of the organism being studied
        version: The version of miRBase to be queried

    """

    # Initialize our logger
    logger = log.setupLogger("housekeeping")

    # Make sure the genome file exists as defined
    if(not os.path.isfile(genomeFilename)):
        logger.error("%s could not be found! Please check that the "\
            "file path was input correctly" % genomeFilename)
        sys.exit()

    # Do not allow execution if both libFilenamesString and libFolder
    # are defined to anything other than empty strings
    if(libFilenamesString and libFolder):
        logger.error("You specified both libFolder and libNamesList, but "\
            "only one can exist. Delete one and try running again")
        sys.exit()

    # Loop through all libraries in libFilenamesList and confirm that they
    # exist before running
    for libName in libFilenamesList:
        if(not os.path.isfile(libName)):
            logger.error("%s could not be found! Please check that the "\
                "file path was input correctly" % libName)
            sys.exit()
    if(len(libFilenamesList) == 1):
        logger.warning("Only one library was provided. While miRador "\
            "can run with this, miRador will not\noutput any miRNAs that are "\
            "predicted outside of any known families as we require\n"\
            "identification in multiple libraries for novel annotation.\nIf "\
            "this organism does not exist yet in miRBase, then no miRNAs "\
            "will be predicted.\nPausing execution for 20 seconds if you "\
            "want to stop this run and add libraries. (Use ctrl+c to stop)\n")
        time.sleep(20)

    if(not shutil.which(bowtiePath)):
        logger.error("bowtie could not be found at the provided path: %s\n"\
            "Correct before running again" % bowtiePath)
        sys.exit()

    if(not shutil.which(bowtieBuildPath)):
        logger.error("bowtie-build could not be found at the provided path: "\
            "%s\nCorrect before running again" % bowtieBuildPath)
        sys.exit()

    if(not shutil.which(einvertedPath)):
        logger.error("einverted could not be found at the provided path: "\
            "%s\nCorrect before running again" % einvertedPath)
        sys.exit()

    if(not shutil.which(perlPath)):
        logger.error("perl could not be found at the provided path: %s\n"\
            "Correct before running again" % perlPath)
        sys.exit()

    ### Create the necessary folders if they don't already exist
    # Create a folderfor genome if it does not exist already
    if not os.path.isdir("genome"):
        os.mkdir('genome')
    # Create a folder for the inverted repeat if it does not exist already
    if(not os.path.isdir("invertedRepeats")):
        os.mkdir("invertedRepeats")
    # Create a folder for the miRBase folder if it does not exist yet
    if(not os.path.isdir("miRBase")):
        os.mkdir("miRBase")

    # If the user has filled the outputFolder option, check to see if it
    # has results from an older run and then delete them
    if(outputFolder):
        # Confirm that the output folder's name is not the same as
        # libFolder. This will ensure nothing of importance is
        # accidentally deleted 
        if(outputFolder == libFolder):
            logger.error("outputFolder and libFolder cannot be the same "\
                "folder. Please rename outputFolder and run again")
            sys.exit()

        # Delete the libs folder if it exists already
        if(os.path.isdir("%s/libs" % outputFolder)):
            shutil.rmtree("%s/libs" % outputFolder)

        # Delete the images folder if it exists already
        if(os.path.isdir("%s/images" % outputFolder)):
            shutil.rmtree("%s/images" % outputFolder)

        # Delete the various output files if they exist already
        if(os.path.isfile("%s/%s_blastResults.txt" % (outputFolder,
                outputFolder))):
            os.remove("%s/%s_blastResults.txt" % (outputFolder, outputFolder))
        if(os.path.isfile("%s/%s_finalAnnotatedCandidates.csv" % \
                (outputFolder, outputFolder))):
            os.remove("%s/%s_finalAnnotatedCandidates.csv" % (outputFolder,
                outputFolder))
        if(os.path.isfile("%s/%s_finalAnnotatedCandidates.fa" % \
                (outputFolder, outputFolder))):
            os.remove("%s/%s_finalAnnotatedCandidates.fa" % (outputFolder,
                outputFolder))
        if(os.path.isfile("%s/%s_preAnnotatedCandidates.csv" % \
                (outputFolder, outputFolder))):
            os.remove("%s/%s_preAnnotatedCandidates.csv" % (outputFolder,
                outputFolder))
        if(os.path.isfile("%s/%s_preAnnotatedCandidates.fa" % \
                (outputFolder, outputFolder))):
            os.remove("%s/%s_preAnnotatedCandidates.fa" % (outputFolder,
                outputFolder))

    ###########################################################################

    # Create a path for an output folder if it does not exist already
    # (Almost certainly shoul dnot as it would require the same run second)
    else:
        outputFolder = datetime.datetime.now().strftime(
            'output_%Y-%m-%d_%H-%M-%S')
    if not os.path.isdir(outputFolder):
        os.mkdir(outputFolder)
    if not os.path.isdir("%s/libs" % outputFolder):
        os.mkdir("%s/libs" % outputFolder)
    if not os.path.isdir("%s/images" % outputFolder):
        os.mkdir("%s/images" % outputFolder)

    mirBaseDict = setupMiRBase.setupMiRBase(organism, version)

    return(mirBaseDict)
