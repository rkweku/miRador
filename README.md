# miRador - Plant miRNA identification tool
[![License](https://img.shields.io/badge/license-Apache%202-4EB1BA.svg)](https://www.apache.org/licenses/LICENSE-2.0.html

## Description
miRador is a miRNA prediction tool developed to be the first of two compontents in a miRNA prediction and validation pipeline. While miRNAs can be predicted quite accurately with miRador alone, utilizing PARE data to predict and validate targets of novel miRNAs with sPARTA is the ultimate use of this package. Through runs of both programs, users will be able to provide evidence of cleavage at predicted targets of novel miRNAs.

## Installation
miRador comes with a conda environment file which can be utilized to ensure all dependencies are satisfied, with versions that miRador was developed for, prior to running. Install either anaconda or miniconda, if you do not have it installed already on your system, following either of the links below:
Anaconda: https://www.anaconda.com/products/distribution
miniconda: https://docs.conda.io/en/latest/miniconda.html

Once conda is installed, a miRador environment can be created with the following command from within the miRador directory:
`conda env create -n mirador -f environment.yml`

When this command is complete, you are ready to run miRador. All of the dependencies in the following section should now be satisfied.

___

In the event that there are issues installing via the provided .yml file, you may also try to create your own conda environment that will be miRador ready. First, create a conda environment:  
`conda create --name mirador`  

Next, activate the envioronment
`conda activate mirador`

Set the conda channel priority to flexible as some environments will fail without this change. We will change it back to strict at the end.
`conda config --set channel_priority flexible`

Then, install the following packages:  
`conda install python=3.7.13`  
`conda install -c conda-forge ghostscript=9.54.0 perl=5.32.1 pypdf2=2.11.1 typing_extensions=4.5.0`  
`conda install -c bioconda biopython=1.78 blast=2.13.0 bowtie=1.3.1 emboss=6.6.0 samtools=1.6 perl-io-string=1.08 viennarna=2.5.1`  

Finally, reset the channel priority back to strict
`conda config --set channel_priority strict`
___
If you opt to not utilize conda, the dependencies for miRador can be downloaded separately and their executable paths can be initialized in miRador.ini

### Dependencies
There are several dependencies of miRador, all of which are checked prior to running. If you opt to not utilize conda, you must download and install the following packages

#### Standalone packages
`blast`: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

`bowtie`: http://bowtie-bio.sourceforge.net/index.shtml

`einverted`: http://emboss.sourceforge.net/download/

`ps2pdfwr`: https://www.ghostscript.com/download.html

`perl`: https://www.perl.org/get.html

`python3`: https://www.python.org/downloads/

`samtools`: htslib.org/download/

`ViennaRNA`: https://www.tbi.univie.ac.at/RNA/#download

#### Python3 packages
`PyPDF2`: https://pypi.org/project/PyPDF2/

Upon python installation, Python packages may be installed via pip. If your system does not have pip installed already, or you do not have pip for your version of python, you can follow the installation instructions here: https://pip.pypa.io/en/stable/installing/.

If you do not have sudo privileges on your system, the --user option can be used when installing packages with pip. This will add the package to your local directory python directory and allow you to install any package without the need of sudo. To do this, open your terminal and simply type: `pip3 install --user PackageName`

#### Perl Modules
`IO::String`: https://metacpan.org/pod/IO::String

This can be installed via CPAN. See instructions here:
https://docs.huihoo.com/livejournal/server/lj.install.perl_setup.modules.html

## Configuration
Running miRador requires modifying a configuration file, initially named `miRador.ini`

### Options
|Option                   |Default |Description                                                                                   |
|-------------------------|--------|----------------------------------------------------------------------------------------------|
|genomeFilename           |        |The path and name of the genome file                                                          |
|runEInvertedFlag         |1       |Flag to be set if you wish to run EInverted                                                   |
|einvertedPresets         |        |Presets for einverted parameters to be set. low, medium, or high                              |
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
|RNAFoldPath              |        |Path of RNAFold                                                                               |
|RNAPlotPath              |        |Path of RNAPlot                                                                               |
|blastnPath               |        |Path of blastn                                                                                |
|makeblastdbPath          |        |Path of makeblastdb                                                                           |
|ps2pdfwrPath             |        |Path of ps2pdfwr                                                                              |
|outputFolder             |        |Name of specific folder to write data to. *If folder exists, data within will be overwritten* |

### Bowtie
Bowtie has been configured to be run with the following arguments. While these can be changed, there is no simple option to do so within the ini file and thus changes would need to be done within library.py. With that said, the specific options are as follows:
- -a to report all valid alignments as we want multihits to the genome. We expect few multi-matches to the genome
- --best and --strata ensures only the best alignments are reported and so that less optimum but passable alignments do not appear
- -v 0 Allow no mismatch
- --sam-nohead removes the header from the SAM file. This is useful because we have to merge the fragment alignments for parallel runs
- --no-unal suppresses sequences with no alignemnt. This helps to keep the map file manageable and filter out these sequences earlier for efficiency

### Einverted
Einverted is utilized to predict an initial set of inverted repeats from the genome FASTA file.
There are three presets which can be selected by the user, low, medium or high, which assigns pre-defined preset scores for matches and mismatches (utilize negative integer), penalty score for gaps (utilize positive integer), and the overall scoring threshold for reporting an inverted repeats. 
low: match = 3, mismatch = -4, gap = 6, threshold = 40
medium: match = 3, mismatch = -4, gap = 6, threshold = 45
high: match = 3, mismatch = -5, gap = 7, threshold = 50

A `maxRepLen` variable is editable as well, though we recommend this remain at its default value of 300.

While we generally recommend running miRador with a medium threshold, you are able to edit the individual scores and penalties yourself with the respective variables in the configuation file. If a preset is used, you may simply override the default values of any, or all scores, by placing a value yourself. For instance, you can use the medium presets, but set `mismatch = -5` to run the pipeline with the medium scores, but the mismatch score will be overridden by your provided score.

## Examples
When all options in the configuration file are set, running miRador is quite simple. From the miRador base directory, type:
```
python3 miRador miRador.ini
```

The configuation file that is included in this repository include the first chromosome of the Arabidopsis TAIR10 genome and two sRNA sequencing libraries which can be used to test that miRador will run properly. These files do need to be unzipped before running. This can be done quite simply with the following command: `gunzip -r test`

## Output
miRador writes its output to a folder provided by you, the user, or to a folder with the data and time as a means of providing a unique folder name. The contents of the folder are:
- finalAnnotatedCandidates.csv: A CSV file containing the annotated miRNAs that miRador was able to identify
- preAnnotatedCandidates.csv: A CSV fil containing the miRNAs that miRAdor was able to identify, but this file does not contain any annotation using the known miRNAs from miRBase. Additionally, no filtration is done for novel miRNAs as miRNA annotation has not been performed at the time of this file's creation.
- preAnnotatedCandidates.fa: Simply a FASTA file containing all of the miRNAs in preAnnotatedCandidates. Generally should not be used for any purposes outside of miRador itself
- images: This folder contains RNAfold predicted images of the precursor miRNAs and the miRNA:miRNA* duplex on those precursors. Note that the images created here may not actually match the einverted predicted inverted repeat as these algorithms are different. However, the images created here can still generally be quite useful when making a decision on a predicted miRNA's status.
- libs: While somewhat useful, this folder is unlikely to be used on a run-to-run basis as its results do not consider annotation of the candidate miRNAs. This folder will contain two files for each library provided for prediction. The first is libName_all_precursors.txt which shows every einverted predicited inverted repeat, its coordinates, the alignment between the two strands of the inverted repeat, and every small RNA that maps to the inverted repeat, including their mapping starting positions and their hits normalized abundances. The second file, libName_candidate_precursors.txt, is a much more filtered version of the first file. In this file, every precursor miRNA and miRNA:miRNA* duplex for that precursor that was identified in this library after all filters (except library replication) is detailed. Here we can find the precursor miRNA's coordinates, the alignment between the two strands of the precursor miRNA, the candidate miRNA and its star sequence as well as their positions and hits normalized abundances, and the specific alignment between the miRNA and the miRNA*, according to einverted.
