# AMA

This is an automated workflow tool for downloading and processing bioinformatic data. Tis workflow was developed by Agata and Max in collaboration with the bioinformatics group BioLab at the University of Melbourne under the supervision of Heroen Verbruggen.


## Introduction

The workflow was designed to be fast and memory-optimized. The download works with an input list of SRA accessions from the NIH database. 10% of this data is downloaded in order to validate the individual SRAs according to the parameters perc_identity, qcov_hsp_perc. Subsequently, promising data sets are entirely downloaded and analyzed. The pipeline was designed to be modular, allowing various bioinformatic tools to be used to preprocess and subsequently analyze the data. By default, Fastp is used in the workflow for preprocessing and then Bast for analyzing and creating a results database. However, any other precrossesing tool and analysis tool can also be used, only the Snakemake workflow has to be adapted. 


## Running the Code

This project is powered by [Nix](https://nixos.org) making it incredibly easy to reproduce. Simply run `nix build` after installing *nix* and cloning the repository, will create the Document in the 'result' folder. If you are unwilling (or unable) to install `nix` running `Rscript src/script.R` should yield the same result if the dependencies are installed correctly

### Dependencies

#### Conda packages

* ggplot2
* tidyverse
* tikzDevice
* gplots
* FactoMineR
* ggfortify
* cluster
* readxl
* factoextra
* pheatmap
