# RWellstone_FSHD_muscle_biopsy

This repository was made to support the transparency and reproducibility of our computational analysis for the manuscript [Longitudinal measures of RNA expression and disease activity in FSHD muscle biopsies](https://doi.org/10.1093/hmg/ddaa031), published in Jan. 2020 by _Human Molecular Genetic_ (HMG).  It has two major branches - main and _gh-pages_. The main branch is a standard R package containing RData (`../data/*.rda`) files, scripts (`../inst/scripts/*.R`) used for the manuscript, and R markdown files (`../inst/gitbook/*.Rmd`) that makes the gitbook, the narratives and codes that back up our discussion on our the analysis. The  _gh-pages_ branch hosts the [gitbook](https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy). 

As a standard R experimental data package, __RWEllstone_FSHD_muscle_biopsy__ master branch consists with `/data`, and `/inst` folders and `DESCRIPTION` file. The `/data` folder contains many instances built upon the and _R_ and _Bioconductor_ platform, especially the DESeq2's RNA-seq count container, `DESeqDataSet`, as listed below.  The histopathology scores and MRI signals are stored in `mri_pathology.rda`.

The major Bioinformatics analysis is performed by using Bioconductor _DESeq_ and _goseq_ packages, and visualization by the R _ggplot2_ packages.

## Package structure
```
- /data     
  |- sanitized.dds: a DESeqDataSet instance of the biopsy RNA-seq gene expression matrix (containing first and the follow-up visits) and MRI characteristics and 
     histopathology socres.  
  |- cluster_df: a data.frame instance of k-means clustering results 
  |- FSHD_markers: a DataFrame instance of 53 DUX4-targeted FSHD robust biomarkers from Yao et al. (2014) paper
  |- marker_path_scores: sample marker scores (DUX4/inflamm/extracellular matrix/cell cycle) 
                         based on gene expression alog with MRI characteristics and 
                         histopathology scores
  |- mri_pathology: a data.frame instance of MRI characteristics and histopathology scores
  |- pax7_targets: a data.frame instance of PAX7-targeted gene expression matrix
  |- year2.dds: a DESeqDataSet instance of the follow-up visit gene expression 
                matrix
  
- /inst
  |- /gitbook: Rmd files that make this gitbook
  |- /extdata: supplemental tables
  |- /scripts: shell and original R scripts
```  

## gitbook
The gitbook is hosted on the _gh-pages_ branch and can be accessed [here](https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy). This book consists with chapters explaining different parts of analysis for the manuscripts. Most of the major codes are listed here. The goad of this book is to support the transparency and reproducibility of the research involved in computational analysis. 

## Software requirement
- R-3.5.1
- Bioconductor 3.7. See the DESCRIPTION file for depend, import and suggest packages.
