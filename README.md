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
|libFilenamesList         |        |List of library file names and their path for each. Should not be set if libFolder is set     |
|libFolder                |        |The name of the folder holding all of the chopped.txt files. Should not be set if libFilenamesList is used|
|organism                 |        |First letter of genus and first 2 letters of species                                          |
|version                  |CURRENT |Version of miRBase to use for annotation                                                      |
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
## Output
miRador writes its output to a folder provided by you, the user, or to a folder with the data and time as a means of providing a unique folder name. The contents of the folder are:
- finalAnnotatedCandidates.csv: A CSV file containing the annotated miRNAs that miRador was able to identify
- preAnnotatedCandidates.csv: A CSV fil containing the miRNAs that miRAdor was able to identify, but this file does not contain any annotation using the known miRNAs from miRBase. Additionally, no filtration is done for novel miRNAs as miRNA annotation has not been performed at the time of this file's creation.
- preAnnotatedCandidates.fa: Simply a FASTA file containing all of the miRNAs in preAnnotatedCandidates. Generally should not be used for any purposes outside of miRador itself
- images: This folder contains RNAfold predicted images of the precursor miRNAs and the miRNA:miRNA* duplex on those precursors. Note that the images created here may not actually match the einverted predicted inverted repeat as these algorithms are different. However, the images created here can still generally be quite useful when making a decision on a predicted miRNA's status.
- libs: While somewhat useful, this folder is unlikely to be used on a run-to-run basis as its results do not consider annotation of the candidate miRNAs. This folder will contain two files for each library provided for prediction. The first is libName_all_precursors.txt which shows every einverted predicited inverted repeat, its coordinates, the alignment between the two strands of the inverted repeat, and every small RNA that maps to the inverted repeat, including their mapping starting positions and their hits normalized abundances. The second file, libName_candidate_precursors.txt, is a much more filtered version of the first file. In this file, every precursor miRNA and miRNA:miRNA* duplex for that precursor that was identified in this library after all filters (except library replication) is detailed. Here we can find the precursor miRNA's coordinates, the alignment between the two strands of the precursor miRNA, the candidate miRNA and its star sequence as well as their positions and hits normalized abundances, and the specific alignment between the miRNA and the miRNA*, according to einverted.
