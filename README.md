# AMA

This is an **snakemake** and **python** automated workflow tool for downloading and processing bioinformatic data. This workflow was developed by Agata and Max in collaboration with the bioinformatics group BioLab at the University of Melbourne under the supervision of Heroen Verbruggen.


## Introduction

The workflow was designed to be fast and memory-optimized. The download works with an input list of SRA accessions from the NIH database. A fraction of this data specified in parameter _download_perc_ is downloaded in order to validate the individual SRAs according to the parameters _perc_identity_, _qcov_hsp_perc_. Subsequently, promising data sets are entirely downloaded and analyzed. The pipeline was designed to be modular, allowing various bioinformatic tools to be used to preprocess and subsequently analyze the data. By default, **Fastp** is used in the workflow for preprocessing and then **Blast** for analyzing and creating a results database. However, any other precrossesing tool and analysis tool can also be used, only the snakemake workflow has to be adapted. 


## Installation

- **Conda**

   - Please checkout the [Conda Documentation](https://github.com/conda/conda-docs).

   - To execute all tasks in one single conda environment the `AMA.yaml` contains all required packages and the corresponding channels
   

     > Note: If you want to update your current environment manually you should add the following **conda packages**:
       >
       > ```bash   
       > - python
       > - snakemake
       > - sra-tools
       > - fastp
       > - blast
       > - curl
       > - wget
       > - entrez-direct
       > ```


   - Otherwise navigate to the location of the pulled AMA.yaml file and execute `conda env create -f AMA.yaml`


   - To activate the created conda environmentconda and getting started execute `conda activate AMA`


## Requirements

- Set up all required data

   - you need the input `csv_file` containing all of the retrieved SRA accessions found with the selected search string within the NIH database search.

   - you need the input `query_fasta` containing the reference sequence which the processing algorithm (in our case **Blast**) needs to compare with the individual sequences from the SRA accessions and calculate their match.



> Note: Before the snakemake workflow can be started, the config files of the various modular steps must be customized with **individual output paths** and **desired processing parameters**



- **partial_download_config**

   - `csv_file:` The path to the csv file containing the results (SRA accessions) generated with the selected search string within the NIH database search `/genbank_store/biodiversity_project/data/SRA_list_example.csv`

   - `output_directory:` The path in which the results shall be saved it is important that the input contains a previously created folder as destination, for example `/genbank_store/biodiversity_project/result1`

   - `download_perc:` The percentage of how much of the datasets should be downloaded for each individual SRR accession we would suggest `10` for the validation download and `100` for the full download later on


- **alignment_config**

   - `csv_file`: Same like above `/genbank_store/biodiversity_project/data/SRA_list_example.csv`

   - `output_directory:` Same like above `/genbank_store/biodiversity_project/result1`

   - `query_fasta:` The path to the reference sequence against Blast compares the individual sequences and calculates their match `/genbank_store/biodiversity_project/data/example_reference.fa`

   - `perc_identity:` Treshold, how much percentage of the sequence should matching `95`

   - `qcov_hsp_perc:` Treshold, how much percentage of the query fasta is covered `90`



## Example usage

> [!IMPORTANT]
> - Activate the previously created conda environment 
> - Make sure the csv_file and query_fasta is provided
> - Adjust all parameters and paths in the config files Helpful advice for doing things better or more easily


> [!TIP]
> - keep in mind that the following code is just an explaination you have to adjust all paths and inputs according to your operation system


- **Partial Download**

   The first step and modular part of the workflow will be the partial download running the skript to validate the datasets fast and memory-optimized with the snakemake syntax:

     >
     > ```bash
     > snakemake --snakefile /genbank_store/biodiversity_project/download/smk_download_new_final.smk --configfile /genbank_store/biodiversity_project/download/partial_download_config.yaml --cores 16
     > ```


- **Processing**

   The next step is the preprocessing and alignment with FastP and Blast there are plenty of paramenters to adjust preprocessing trimming and quality control of the input sequences which can change the output significant if you want to adjust some variables have a look at [FastP] (https://github.com/OpenGene/fastp) and [Blast] (https://blast.ncbi.nlm.nih.gov/Blast.cgi) and change it in the `smk_alignment2_final.smk`:

     >
     > ```bash
     > snakemake --snakefile /genbank_store/biodiversity_project/workflow/smk_alignment2_final.smk --configfile /genbank_store/biodiversity_project/workflow/alignment_config.yaml --cores 16
     > ```


- **Full Download**



## Output Structure


### Dependencies
