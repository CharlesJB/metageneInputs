# Dependencies

Make sure those dependencies are installed **before** you try to install the
`metageneInputs` package:

```r
require(devtools)
devtools::install_github("CharlesJB/FantomEnhancers.hg19", ref = "v0.2.0")
devtools::install_github("CharlesJB/FantomTSS.hg19", ref = "v0.2.0")
devtools::install_github("CharlesJB/metagene",
                         ref = "05b1e7594a72901b8ccd37eead171f56017693bc")
install_packages("dplyr")
install_packages("tidyr")
```

All other dependencies should be installed automatically during `metageneInputs`
installation.

# Introduction

The goal of this package is to show how to reproduce the input files that will
be used for the Plos Computationnal Biology paper:

* Groups of regions
* Scripts to download bam files (Makefile)
* Designs

# Usage

First, you need to install the package. This could take ~10-30 minutes
depending on your machine since the whole analysis will be done during the
package installation (i.e.: none of the results of the analysis are saved on
the github repository).

```r
require(devtools)
# TODO: Add ref when available
devtools::install_github("CharlesJB/metageneInputs", build_vignettes = TRUE)
```

To avoid accidently overwritting files with similar names to those produced by
this package, the results are saved in a temporary directory. The name of the
directory changes each time the package is built. The name of the directory
used for the analysis is saved in the vignettes. To view the vignettes, from a
R session:
```r
browseVignettes("metageneInputs")
```

# Input files descriptions

## `splitted_regions.RData`

This object contains 2 `GRangesList`s: one for enhancers and one for promoters.
Enhancers are regions are obtained with the `FantomEnhancers.hg19` package and
promoters are obtained with the `metagene` package.

Each element of the `GRangesList` correpond to an expression group in GM12878
cell lines (based on Famtom5 datasets):

* (0: not expressed)
* (25: First quartile)
* (50: Second quartile)
* (75: Third quartile)
* (100: Fourth quartile)

We calculate the tresholds of expression of each group (in
transcript per million (TPM)) using the non-zero values from every Fantom5
transcription start sites (TSS; using the `FantomTSS.hg19` package) for every
available cell lines except GM12878 and splitting them in quartiles. Regions
from GM12878 cell lines are then splitted in each groups based on their TPM
value.

## `Makefile`

We need to download ENCODE datasets from the GM12878 cell line (everything but
histones). Using the `ENCODExplorer` package, we extract the URLs for every
files and save them in a Makefile. We also saved the URL for the bam files to
use for the validation in the same Makefile.

To download the files, simply copy the Makefile in the directory where you want 
to download the files and type in a shell:
```shell
make
```

This will download all the bam files required to reproduce the metagene analysis
(~300G).

## `designs.RData`

This objects contains 2 elements:
* `designs_encode`: A list with the designs for every ENCODE experiments to be
analyzed.
* `design_validation`: The design for the validation.

## `tpm_grangeslist`

This file contains a `GRangesList` named `TPM` that consists of 2 `GRanges`, one
for the promoters and one for the enhancers. Each `GRanges` contains the mean
TPM value for each regions.
