# AMA

This is an automated workflow tool for downloading and processing bioinformatic data. Tis workflow was developed by Agata and Max in collaboration with the bioinformatics group BioLab at the University of Melbourne under the supervision of Heroen Verbruggen.


## Introduction

The workflow was designed to be fast and memory-optimized. The download works with an input list of SRA accessions from the NIH database. 10% of this data is downloaded in order to validate the individual SRAs according to the parameters _perc_identity_, _qcov_hsp_perc_. Subsequently, promising data sets are entirely downloaded and analyzed. The pipeline was designed to be modular, allowing various bioinformatic tools to be used to preprocess and subsequently analyze the data. By default, **Fastp** is used in the workflow for preprocessing and then **Blast** for analyzing and creating a results database. However, any other precrossesing tool and analysis tool can also be used, only the snakemake workflow has to be adapted. 


## Requirements

- **Conda**

1. Please checkout the [Conda Documentation](https://github.com/conda/conda-docs).

  2. To execute all tasks in one single conda environment the `AMA.yaml` contains all required packages and the corresponding channels
   
     - If you want to update your current environment manually you should add the following **conda packages**:
       
       - python
       - snakemake
       - sra-tools
       - fastp
       - blast
       - curl
       - wget
       - entrez-direct

  3. Otherwise navigate to the location of the pulled AMA.yaml file and execute `conda env create -f AMA.yaml`


  4. To activate the created conda environmentconda and getting started execute `conda activate AMA`




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
