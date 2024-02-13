# AMA

This is an **snakemake** and **python** automated workflow tool for downloading and processing bioinformatic data. This workflow was developed by Agata and Max in collaboration with the bioinformatics group BioLab at the University of Melbourne under the supervision of Heroen Verbruggen.


## Introduction

The workflow was designed to be fast and memory-optimized. The download works with an input list of SRA accessions from the NIH database. A fraction of this data specified in parameter _download_perc_ is downloaded in order to validate the individual SRAs according to the parameters _perc_identity_, _qcov_hsp_perc_. Subsequently, promising data sets are entirely downloaded and analyzed. The pipeline was designed to be modular, allowing various bioinformatic tools to be used to preprocess and subsequently analyze the data. By default, **Fastp** is used in the workflow for preprocessing and then **Blast** for analyzing and creating a results database. However, any other precrossesing tool and analysis tool can also be used, only the snakemake workflow has to be adapted. 


## Installation

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


## Requirements

Before the snakemake workflow can be started, the configfiles of the various modular steps must be customized with **individual output paths** and **desired processing parameters**.

- **Config_download**

    - `csv_file:` The path to the csv file containing the results (SRA accessions) generated out of the search string from the NIH database search `/genbank_store/biodiversity_project/data/SRA_list_example.csv`

    - `output_directory:` The path including a created folder in which the results shall be saved, for example `/genbank_store/biodiversity_project/result1`

    - `download_perc:` The percentage of how much of the datasets should be downloaded for each individual SRR accession we would suggest `10` for the validation download and `100` for the full download later on


- **Config_alignment**

    - `csv_file`: Same like above `/genbank_store/biodiversity_project/data/SRA_list_example.csv`

    - `output_directory:` Same like above `/genbank_store/biodiversity_project/result1`

    - `query_fasta:` The path to the reference sequence with which the blast search compares the sequences from the NIH database `/genbank_store/biodiversity_project/data/example_reference.fa`

    - `perc_identity:` Treshold, how much percentage of the sequence should matching `95`

    - `qcov_hsp_perc:` Treshold, how much percentage of the query fasta is covered `90`



## Example usage

- **getting started**
- **downloading**
- **processing**



## Output Structure


### Dependencies
