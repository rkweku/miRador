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
|Option                   |Default |Description                                                                                   |
|-------------------------|--------|----------------------------------------------------------------------------------------------|
|genomeFilename           |        |The path and name of the genome file                                                          |
|runEInvertedFlag         |1       |Flag to be set if you wish to run EInverted                                                   |
|gap                      |6       |EInverted score for gaps                                                                      |
|match                    |3       |EInverted score for matches                                                                   |
|mismatch                 |-4      |EInverted penalty score for mismatches                                                        |
|threshold                |40      |Einverted scoring threshold for identifying inverted repeats                                  |
|maxRepLen                |300     |Maximum length that an inverted repeat can be                                                 |
|libFilenamesList         |        |List of library file names and their path for each                                            |
|libFolder                |        |The name of the folder holding all of the chopped.txt files                                   |
|organism                 |        |First letter of genus and first 2 letters of species                                          |
|parallel                 |        |Flag to utilize parallelization                                                               |
|nthrads                  |        |Number of threads to utilize when running bowtie                                              |
|bowtiePath               |        |Path of bowtie                                                                                |
|bowtieBuildPath          |        |Path of bowtie-build                                                                          |
|einvertedPath            |        |Path of einverted                                                                             |
|perlPath                 |        |Path of perl                                                                                  |
|outputFolder             |        |Name of specific folder to write data to. *If folder exists, data within will be overwritten* |

## Examples
When all options in the configuration file are set, running miRador is quite simple. From the miRador base directory, type:
```
python3 miRador miRador.ini
```
