import fnmatch
import ftplib
import gzip
import os
import shutil
import sys
import urllib.request

import annotateCandidates

def checkNeedUpdate(version):
    """Check if the version file from miRBase exists in our miRBase folder.

    Args:
        version: Version of miRBase to be queried. Generally should
            be "CURRENT"

    Returns:
        True if we need to update our miRBase folder and false if no update
        is required

    """

    ftp = ftplib.FTP("mirbase.org")
    ftp.login("anonymous", "")
    ftp.cwd("pub/mirbase/%s/" % version)

    # Pull the list of files from miRBase
    filenamesList = ftp.nlst()

    # Loop through the files in this version of miRBase to find the
    # version file
    for filename in filenamesList:
        if(filename.startswith("0_THIS_IS_RELEASE")):
            versionFile = filename

        # If the verson file exists exactly as it was found in miRBase,
        # return false as we do not need to update our miRBase files
        if(os.path.exists("miRBase/%s" % filename)):
            ftp.close()
            return(False)

    ftp.close()
    return(True)

def downloadOrganismsAndMirnas(version):
    """Get the organisms.txt.gz and mature.fa.gz files from the current
    version of miRBase, unzip the files, and save them to the  miRBase 
    folder

    Args:
        version: Version of miRBase to be downloaded. Generally should
            be "CURRENT"

    """

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
        print("Error while downloading organisms.txt.gz...\n%s" % e)
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
            "miRBase/miRBaseMirnas.fa.gz", "wb").write)

    # If there is no file 
    except ftplib.all_errors as e:
        print("Error while downloading mature.fa.gz...\n%s" % e)
        sys.exit()

    # Use gzip to unzip the file and save it as miRBaseMirnas.fa
    with gzip.open("miRBase/miRBaseMirnas.fa.gz", "rb") as f_gz:
        with open("miRBase/miRBaseMirnas.fa", "wb") as f_unzip:
            shutil.copyfileobj(f_gz, f_unzip)

    # Remove organisms.txt.gz as we have already unzipped it
    os.remove("miRBase/miRBaseMirnas.fa.gz")

    ftp.close()

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
        name = organism[2]
        tree = organism[3]

        if("viridiplantae" in tree.lower()):
            plantList.append((threeLetterIdentifier, name))

    return(plantList)

def downloadPlantSpecies(version, plantList):
    """Download all plant species from the current version of miRBase

    Args:
        version: Version of miRBase to be downloaded. Generally should
            be "CURRENT"
        plantList: A list of tuples of the three letter identifier and
            organism name of all plant species in miRBase

    """

    # Connect to the miRBase FTP page
    ftp = ftplib.FTP("mirbase.org")
    ftp.login("anonymous", "")
    ftp.cwd("pub/mirbase/%s/genomes/" % version)

    # Loop through the organisms to download their gff3 files from miRBase
    for organism in plantList:
        threeLetterIdentifier = organism[0]
        name = organism[1]

        # Try to download the gff3 file for this organism.
        try:
            ftp.retrbinary("RETR %s.gff3" % threeLetterIdentifier, open(
                "miRBase/%s.gff3" % threeLetterIdentifier, "wb").write)

        # If there is an error when opening an ftp file, it is because the
        # GFF file does not exist for this organism. This seems to be
        # caused by an organism not having a sequenced genome 
        except ftplib.all_errors as e:
            os.remove("miRBase/%s.gff3" % threeLetterIdentifier)

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
            name = attributes.split("Name=")[1].split(";")[0]

            # Check to ensure a miRNA is not annotated twice in the GFF file.
            # Report an error to the user if there is as there should not be
            if(name in mirBaseDict):
                print("Error. Found more than one miRNA with the same name.\n"\
                    "%s was found twice in %s" % (name, filename))
                sys.exit()

            # Save the miRNA in mirBaseDict with the name as the key and its
            # coordinates as the value
            mirBaseDict[name] = (chrID, strand, start, end)

    return(mirBaseDict)

def setupMiRBase(organism, version):
    """Download all plant miRNAs from the provided version of miRBase to
    generate the subject files for BLAST when attentmping to annotate
    our candidate miRNAs

    Args:
        organism: The three letter identifier of the organism being studied
        version: Version of miRBase to the downloaded. Generally should 
            be "CURRENT"

    """

    updateStatus = checkNeedUpdate(version)
    gffFilename = "miRBase/%s.gff3" % organism
    mirBaseDict = {}

    # Check if the miRBase files need to be updated
    if(updateStatus):
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

    # We can only perform identity searches, with positional information, if
    # the GFF file exists. So, first check if the GFF file actually exists,
    # then read the file into a dictionary with parsePrecursorGFF
    if(os.path.isfile(gffFilename)):
        mirBaseDict = parsePrecursorGFF(gffFilename)

    return(mirBaseDict)
