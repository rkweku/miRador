## miRador Plant miRNA identification tool
[![License](https://img.shields.io/badge/license-Apache%202-4EB1BA.svg)](https://www.apache.org/licenses/LICENSE-2.0.html

## Description
miRador is a miRNA prediction tool developed to be the first of two compontents in a miRNA prediction and validation pipeline. While miRNAs can be predicted quite accurately with miRador alone, utilizing PARE data to predict and validate targets of novel miRNAs with sPARTA is the ultimate use of this package. Through runs of both programs, users will be able to provide evidence of cleavage at predicted targets of novel miRNAs.

## Dependencies
There are several dependencies of miRador, all of which are checked prior to running. However, in order to meet all dependencies, users should download the following:

### Standalone packages
`bowtie`: http://bowtie-bio.sourceforge.net/index.shtml
`perl`: https://www.perl.org/get.html
`python3`: https://www.python.org/downloads/
`ViennaRNA`: https://www.tbi.univie.ac.at/RNA/#download

### Python3 packages
`biopython`: https://biopython.org/wiki/Download

## Configuration
Running miRador requires modifyinga configuration file, initially named `miRador.ini`

### Options
