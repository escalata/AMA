# AMA

This is an **snakemake** and **python** automated workflow tool for downloading and processing bioinformatic data. This workflow was developed by Agata and Max in collaboration with the bioinformatics lab [Verbruggen Laboratory](https://hverbruggen.github.io//) at the University of Melbourne under the supervision of Heroen Verbruggen.


## Introduction

The workflow was designed to be fast and memory-optimized. The download works with an input list of SRA accessions from the NIH database. A fraction of this data specified in parameter _download_perc_ is downloaded in order to validate the individual SRAs according to the parameters _perc_identity_, _qcov_hsp_perc_. Subsequently, promising data sets are entirely downloaded and analyzed. The pipeline was designed to be modular, allowing various bioinformatic tools to be used to preprocess and subsequently analyze the data. By default, **Fastp** is used in the workflow for preprocessing and then **Blast** for analyzing and creating a results database. However, any other precrossesing tool and analysis tool can also be used, only the snakemake workflow has to be adapted. 


## Installation

- **Conda**

   - Please check out the [Conda Documentation](https://github.com/conda/conda-docs).

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
       > - bc
       > - entrez-direct       
       > ```


   - Otherwise navigate to the location of the pulled AMA.yaml file and execute `conda env create -f AMA.yaml`


   - To activate the created conda environmentconda and getting started execute `conda activate AMA`


## Requirements

- Set up all required data

   - you need the input `csv_file` containing all of the retrieved SRA accessions found with the selected search string within the NIH database search.

   - you need the input `query_fasta` containing the reference sequence which the processing algorithm (in our case **Blast**) needs to compare with the individual sequences from the SRA accessions and calculate their match.



> Note: Before the snakemake workflow can be started, the config files of the various modular steps must be customized with **individual output paths** and **desired processing parameters**



- **config_download**

   - `csv_file:` The path to the csv file containing the results (SRA accessions) generated with the selected search string within the NIH database search `/genbank_store/biodiversity_project/data/SRA_list_example.csv`

   - `output_directory:` The path in which the results shall be saved it is important that the input contains a previously created folder as destination, for example `/genbank_store/biodiversity_project/result1`

   - `download_perc:` The percentage of how much of the datasets should be downloaded for each individual SRR accession we would suggest `10` for the validation download and `100` for the full download later on


- **config_alignment**

   - `csv_file`: It is important to reference the csv file `downloaded_sra_ids.csv` created by the download script here, as this is an updated version of the original csv file in which the failed downloads have already been deleted. `/genbank_store/biodiversity_project/result1/downloaded_sra_ids.csv`

   - `output_directory:` Same like in the config_download `/genbank_store/biodiversity_project/result1`

   - `query_fasta:` The path to the reference sequence against Blast compares the individual sequences and calculates their match `/genbank_store/biodiversity_project/data/example_reference.fa`

   - `perc_identity:` Threshold, how much percentage of the sequence should matching `90`

   - `qcov_hsp_perc:` Threshold, how much percentage of the query fasta is covered `90`



## Example usage

> [!IMPORTANT]
> - Activate the previously created conda environment 
> - Make sure the csv_file and query_fasta is provided
> - Adjust all parameters and paths in the config files Helpful advice for doing things better or more easily


> [!TIP]
> - Keep in mind that the following code is just an explaination you have to adjust all paths and inputs according to your operation system


- **Partial Download**

   The first step and modular part of the workflow will be the partial download running the snakefile `download.smk` skript with `download_perc:` parameter 10-20 to validate the datasets fast and memory-optimized with the snakemake bash syntax:

     >
     > ```bash
     > snakemake --snakefile /genbank_store/biodiversity_project/download/download.smk --configfile /genbank_store/biodiversity_project/download/config_download.yaml --cores 16
     > ```


- **Processing**

   The next step is the preprocessing and alignment with FastP and Blast there are plenty of paramenters to adjust preprocessing trimming and quality control of the input sequences which can change the output significant if you want to adjust some variables have a look at [FastP](https://github.com/OpenGene/fastp) and [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and change it in the snakefile `alignment.smk` script then run it with the snakemake bash syntax:

     >
     > ```bash
     > snakemake --snakefile /genbank_store/biodiversity_project/workflow/alignment.smk --configfile /genbank_store/biodiversity_project/workflow/config_alignment.yaml --cores 16
     > ```


- **Delete**

   If you want to rerun the processing on the same downloads again. For example try again with different parameter. You can do it by simply executing the provided python script `delete.py`. It will delete all previous generated processing workflow results by the alignment script but dont harm the downloads. Afterwards just rerun the snakemake alignment script again. The python script is located in the alignment directory, just run it with the bash syntax:

     >
     > ```bash
     > python delete.py
     > ```


- **Full Download**

   If your data looks promising the next part of the modular workflow will be the full download running the snakefile `download.smk` skript with `download_perc:` parameter 100 to get access to the whole dataset. Also keep in mind to change the `output_directory:`of the config_download to a new created empty folder to don`t mix up with the validation download. Run the script with the snakemake bash syntax:

     >
     > ```bash
     > snakemake --snakefile /genbank_store/biodiversity_project/download/download.smk --configfile /genbank_store/biodiversity_project/download/config_download.yaml --cores 16
     > ```


- **Primer Detection


## Output structure

The pipeline will create folders per SRA run accessions and generate results using the run accession as the prefix (**XXX**). The download and analysis of data will be stored in the given superfolder output directory `results1`. Each modular part of the analysis workflow is stored in their corresponding folders within the single run accessions. This makes the data clean seperated and easy to change a single part of the workflow. The results of the analysis workflow among all SRA run accessions are processed and stored and visualized in the `alignment_results.txt`.

 
### The complete outline of directory structure picture:

![Output structure and complete file hirarchie](https://github.com/escalata/AMA/blob/main/Picture_complete_hirarchieV2.png)


## Workflow

### The complete bioinformatic workflow picture:

<img src="https://github.com/escalata/AMA/blob/main/Picture_workflowV2.png" width="600" height="1360">


### Dependencies

- [Verbruggen Laboratory](https://hverbruggen.github.io//)

- [ARA](https://github.com/maurya-anand/ARA/tree/main)

- [FastP](https://github.com/OpenGene/fastp)

- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

- [Conda](https://github.com/conda/conda-docs)
