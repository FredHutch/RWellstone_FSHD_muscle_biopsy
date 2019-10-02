# RWellstone_FSHD_muscle_biopsy

This repository supports the transparent and reproducible analysis for the manuscript "Longitudinal measures of RNA expression and disease activity in FSHD muscle biopsies".  It has two major branches: the master branch is a starndard R package containing RData (`../data/*.rda`) files and scripts (`../inst/scripts/*.R`) used for the manuscript; the _gh-pages_ branch hosts a [gitbook](https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy)  explaining the details of the analysis and R codes. 

As a standard R experimental data package, RWEllstone_FSHD_muscle_biopsy main branch consists with `\data`, and `\inst` folders and `DESCRIPTION` file. The `\data` folder contains mainly the RNA-seq count data wrapped within a `DESeqDataSet` instance, `sanitized.dds.rda`. The histopathology scores and MRI signals are stored in `mri_pathology.rda` and as column data embadded in `sanitized.dds` instance.

The major Bioinformatics analysis is performed by using Bioconductor _DESeq_ and _goseq_ packages, and visualization by the R _ggplot2_packages.

## gitbook
The gitbook is hosted on the _gh-pages_ branch and can be accessed [here](https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy). This book consists with chapters explaining different parts of analysis for the manuscripts. Most of the major codes are listed here. The goad of this book is to support the transparency and reporducibily of the research involved in computational analysis. 

## Software requirement
- R-3.5.1 or above
- Bioconductor 3.7. See the DESCRIPTION file for depend, import and suggest packages.

## Dataset
