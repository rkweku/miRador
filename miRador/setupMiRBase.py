import ftplib
import gzip
import os
import shutil
import sys

import annotateCandidates
import log

def checkNeedUpdate(version):
    """Check if the version file from miRBase exists in our miRBase folder.

    Args:
        version: Version of miRBase to be queried. Generally should
            be "CURRENT"

    Returns:
        True if we need to update our miRBase folder and false if no update
        is required

    """

    # Initialize our logger
    logger = log.setupLogger("checkNeedUpdate")

    ftp = ftplib.FTP("mirbase.org")
    ftp.login("anonymous", "")
    try:
        ftp.cwd("pub/mirbase/%s/" % version)
    except:
        logger.error("Input version does not appear to exist in miRBase. " \
            "Check version number in ini file and try again")
        sys.exit()

    # Pull the list of files from miRBase
    filenamesList = ftp.nlst()

    # Loop through the files in this version of miRBase to find the
    # version file
    for filename in filenamesList:
        # First, check if the miRBase plant file even exists yet. If not,
        # don't bother going any further. We need to update
        if(os.path.exists("miRBase/miRBasePlantMirnas.fa")):
            if(filename.startswith("0_THIS_IS_RELEASE")):
                versionFile = filename

                # If the verson file exists exactly as it was found in miRBase,
                # return false as we do not need to update our miRBase files
                if(os.path.exists("miRBase/%s" % filename)):
                    ftp.close()
                    return(False)

    ftp.close()

    log.closeLogger(logger)

    return(True)

def downloadOrganismsAndMirnas(version):
    """Get the organisms.txt.gz and mature.fa.gz files from the current
    version of miRBase, unzip the files, and save them to the  miRBase 
    folder

    Args:
        version: Version of miRBase to be downloaded. Generally should
            be "CURRENT"

    """

    # Initialize our logger
    logger = log.setupLogger("downloadOrganismsAndMirnas")

    ftp = ftplib.FTP("mirbase.org")
    ftp.login("anonymous", "")
    ftp.cwd("pub/mirbase/%s/" % version)

    # Before we actually get the organism file, we will actually need to
    # downlaod the version file to our directory to prevent updates in
    # successive runs
    filenamesList = ftp.nlst()
    for filename in filenamesList:
        if(filename.startswith("0_THIS_IS_RELEASE")):
            ftp.retrbinary("RETR %s" % filename, open("miRBase/%s" % \
                filename, "wb").write) 

    # Try to download the organisms file for this organism.
    try:
        ftp.retrbinary("RETR organisms.txt.gz", open(
            "miRBase/organisms.txt.gz", "wb").write)

    # If there is no file 
    except ftplib.all_errors as e:
        logger.error("Error while downloading organisms.txt.gz...\n%s\nTry "\
            "changing the version of miRBase you are trying to download" % e)
        sys.exit()

    # Use gzip to unzip the file and save it as organisms.tsv
    with gzip.open("miRBase/organisms.txt.gz", "rb") as f_gz:
        with open("miRBase/organisms.tsv", "wb") as f_unzip:
            shutil.copyfileobj(f_gz, f_unzip)

    # Remove organisms.txt.gz as we have already unzipped it
    os.remove("miRBase/organisms.txt.gz")

    # Try to download the mature miRNA file for this organism
    try:
        ftp.retrbinary("RETR mature.fa.gz", open(
            "miRBase/mature.fa.gz", "wb").write)

    # If there is no file 
    except ftplib.all_errors as e:
        logger.error("Error while downloading mature.fa.gz...\n%s\nTry "\
            "changing the version of miRBase you are trying to download" % e)
        sys.exit()

    # Use gzip to unzip the file and save it as mature.fa
    with gzip.open("miRBase/mature.fa.gz", "rb") as f_gz:
        with open("miRBase/mature.fa", "wb") as f_unzip:
            shutil.copyfileobj(f_gz, f_unzip)

    # Remove organisms.txt.gz as we have already unzipped it
    os.remove("miRBase/mature.fa.gz")

    ftp.close()

    log.closeLogger(logger)

def findPlantSpeciesFromOrganisms():
    """Create a list of plant species using the data available in
    organisms.txt

    Returns:
        List of tuples indicating the three letter identifier and 
        genus + species of all plant organisms in miRBase.

    """

    plantList = []

    # Use the readFile function in annotateCandidates to read the organism
    # file, which is a tsv file
    organismsFile = annotateCandidates.readFile("miRBase/organisms.tsv", "\t")

    # Read through each line of the organismsFile to find all plant species,
    # but skip the heder
    for organism in organismsFile[1:]:
        threeLetterIdentifier = organism[0]
        tree = organism[3]

        if("viridiplantae" in tree.lower()):
            plantList.append(threeLetterIdentifier)

    return(plantList)

def downloadPlantSpecies(version, plantList):
    """Download all plant species from the current version of miRBase

    Args:
        version: Version of miRBase to be downloaded. Generally should
            be "CURRENT"
        plantList: Three letter identifiers of all plant species in miRBase

    """

    # Connect to the miRBase FTP page
    ftp = ftplib.FTP("mirbase.org")
    ftp.login("anonymous", "")
    ftp.cwd("pub/mirbase/%s/genomes/" % version)

    # Loop through the organisms to download their gff3 files from miRBase
    for organism in plantList:
        threeLetterIdentifier = organism

        filename = "%s.gff3" % threeLetterIdentifier

        # Try to download the gff3 file for this organism.
        try:
            with open("miRBase/%s" % filename, 'wb') as fhandle:
                ftp.retrbinary("RETR " + filename, fhandle.write)

        # If there is an error when opening an ftp file, it is because the
        # GFF file does not exist for this organism. This seems to be
        # caused by an organism not having a sequenced genome 
        except ftplib.all_errors as e:
            os.remove("miRBase/%s" % filename)

    ftp.close()

def parsePrecursorGFF(filename):
    """This function parses a GFF3 file from any miRBase's genome folder.
    These GFF files contain the coordinates of all primary miRNA
    transcripts, the mature miRNA from that precursor, and in some
    cases, the miRNA* from that transcript.

    Args:
        filename: Name of the GFF file to parse

    Returns:
        Dictionary with the miRNA name as a key and a tuple as its value.
        This tuple includes the coordinates of the miRNA with the chromosome
        name, strand, start, and end positions of the miRNA

    """

    mirBaseDict = {}

    # Read the GFF file into memory so that we can begin to parse it
    gffFile = annotateCandidates.readFile(filename, "\t")

    for line in gffFile:
        # Skip line that being with # as these are just comments
        if(line[0].startswith("#")):
            continue

        # There 6 fields that we need to track from each entry
        chrID = line[0]
        type = line[2]
        start = line[3]
        end = line[4]
        strand = line[6]
        attributes = line[8]

        # We are only using the miRNAs as we need to annotate them
        # based upon their position, so only add sequences to
        # mirBaseDict if it is in a line with type miRNA
        if(type == "miRNA"):
            # We will split the name on 'Name', but it should be abundantly
            # clear here that this will NOT work on non GFF3 files in miRBaes
            # as those appear to use ID rather than Name, and they also put
            # the names in quotes. If miRBase also chooses to change the
            # field, this could result in a failure
            name = attributes.split("Name=")[1].split(";")[0]

            # Add the miRNA and its coordinates to mirBaseDict
            if(name in mirBaseDict):
                mirBaseDict[name].append((chrID, strand, start, end))
            else:
                mirBaseDict[name] = [(chrID, strand, start, end)]

    return(mirBaseDict)

def createFastaFile(plantList):
    """Loop through mature.fa and add all plant miRNAs to a file called
    miRBasePlantMirnas.fa

    Args:
        plantList: A list of tuples containing the three letter identifiers
            of all plants in miRBase

    """

    f = open("miRBase/mature.fa", "r")
    wholeFile = f.readlines()
    f.close()

    g = open("miRBase/miRBasePlantMirnas.fa", "w")

    # Loop through the file in steps of 2
    for i in range(0,len(wholeFile),2):
        # Strip the > and newline character from seqID
        seqID = wholeFile[i][1:].rstrip()
        # Strip the newline character from the sequence
        sequence = wholeFile[i+1].rstrip()

        if(seqID.split("-")[0] in plantList):
            g.write(">%s\n%s" % (seqID, sequence))
            if(i + 1 != len(wholeFile) - 1):
                g.write("\n")

    g.close() 

def setupMiRBase(organism, version):
    """Download all plant miRNAs from the provided version of miRBase to
    generate the subject files for BLAST when attentmping to annotate
    our candidate miRNAs

    Args:
        organism: The three letter identifier of the organism being studied
        version: Version of miRBase to be downloaded. Generally should 
            be "CURRENT"

    """

    # Initialize our logger
    logger = log.setupLogger("setupMirBase")

    updateStatus = checkNeedUpdate(version)
    gffFilename = "miRBase/%s.gff3" % organism
    mirBaseDict = {}

    # Check if the miRBase files need to be updated
    if(updateStatus):
        logger.info("Downloading the relevant miRBase files")
        # If there are still files in the miRBase folder but it needs to
        # be updated, then clear the contents prior to populating it again
        if(os.listdir("miRBase")):
            shutil.rmtree("miRBase")
            os.mkdir("miRBase")

        # Download the organism file and mature miRNA file
        downloadOrganismsAndMirnas(version)
        # Find all the plant species that are in miRBase and then download
        # their GFF files IF they exist
        plantList = findPlantSpeciesFromOrganisms()
        downloadPlantSpecies(version, plantList)
        createFastaFile(plantList)

    # We can only perform identity searches, with positional information, if
    # the GFF file exists. So, first check if the GFF file actually exists,
    # then read the file into a dictionary with parsePrecursorGFF
    if(os.path.isfile(gffFilename)):
        mirBaseDict = parsePrecursorGFF(gffFilename)

    log.closeLogger(logger)

    return(mirBaseDict)
