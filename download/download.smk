import os 
import shutil
from datetime import datetime

configfile: "config_download.yaml"


# Define global variables from configuration file
CSV_FILE = config["csv_file"]
OUTPUT_DIRECTORY = config["output_directory"]
PERC = config["download_perc"]
# Globale variable for starttime
start_time = datetime.now()


# Function to extract SRA IDs from the CSV file
def extract_sra_ids(csv_file):
    sra_ids = []
    with open(csv_file, 'r') as file:
        # Iterate through each line in the file to extract SRA IDs
        for line in file:
            sra_id = line.strip()
            # Add SRA ID to the list if it's not empty and doesn't start with "acc"
            if sra_id and not sra_id.startswith("acc"):
                sra_ids.append(sra_id)
    return sra_ids


# Rule to define the final output files of the workflow
rule all:
    input:
        expand("{output_directory}/download_summary.txt", output_directory=OUTPUT_DIRECTORY)


# Define the handlers
onstart:
    print("Starting SRA download workflow.")

onsuccess:
    shell("""
        echo 'Download Summary' > {config[output_directory]}/download_summary.txt
        echo '---------------------------------' >> {config[output_directory]}/download_summary.txt
        echo "Successful downloads: $(cat {config[output_directory]}/success_sra_ids.txt | wc -l)" >> {config[output_directory]}/download_summary.txt
        echo "Failed downloads: $(cat {config[output_directory]}/failed_sra_ids.txt | wc -l)" >> {config[output_directory]}/download_summary.txt
        cat {config[output_directory]}/download_summary.txt
    """)
    # Calculate the runtime
    end_time = datetime.now()
    duration = end_time - start_time
    hours, remainder = divmod(duration.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    # Print the runtime
    print("---------------------------------")
    print(f"Total runtime: {hours} hours, {minutes} minutes, {seconds} seconds")
    
onerror:
    print("Error encountered during SRA download workflow.")


# Rule to create directories for each ID and initialize metadata files
rule create_directories_and_initialize_metadata:
    input:
        csv_file = CSV_FILE
    output:
    params:
        output_directory = OUTPUT_DIRECTORY
    shell:
        """
        output_file_label=$(basename {input.csv_file} .csv)
        # output_file_label=$(basename {params.output_file} .fa)
        sra_id={wildcards.sra_id}
        metadata_dir={params.output_directory}/$sra_id/metadata
        mkdir -p $metadata_dir
        metadata_file=$metadata_dir/$output_file_label.metadata.txt
        echo "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash" > $metadata_file
        # touch {output.metadata_initialized}
        """


# Rule to process each SRA ID: download metadata and partial datasets
rule process_sra_id:
    input:
        csv_file = CSV_FILE
    output:
        fastq = "{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq",
        metadata_file = "{output_directory}/{sra_id}/metadata/{sra_id}.run.metadata.txt"
    params:
        output_directory = OUTPUT_DIRECTORY,
        csv_file = CSV_FILE,
        perc = PERC,  # Percentage of data to download
        retries = 3,  # Number of retry attempts
        wait_time = 5  # Wait time between retries in seconds  
    shell:
        """
        # Initialize Variables
        sra_id={wildcards.sra_id}
        metadata_dir={params.output_directory}/$sra_id/metadata
        individual_metadata_file=$metadata_dir/$sra_id.run.metadata.txt
        output_dir={params.output_directory}/$sra_id/raw_fastq


        # Create Directories
        mkdir -p $metadata_dir
        mkdir -p $output_dir


        # Dummy Structur function for unsuccessful downloads
        create_dummy_structure() {{
            touch $metadata_dir/$sra_id.run.metadata.txt
            touch $output_dir/$sra_id"_1.fastq"            
        }}


        # Download Metadata
        retries=0
        success=0
        while [ $retries -lt {params.retries} ]; do
            metadata=$(efetch -db sra -id "$sra_id" -format runinfo | tail -n +2)
            if [ $? -eq 0 ] && [ -n "$metadata" ]; then
                echo "$metadata" >> $individual_metadata_file
                echo "Metadata successfully downloaded for: $sra_id"
                success=1
                break
            else
                echo "Attempt $(($retries + 1)) of {params.retries} failed to download metadata for: $sra_id"
                retries=$(($retries + 1))
                sleep {params.wait_time}
            fi
        done

        if [ $success -eq 0 ]; then
            echo "Failed to download metadata for: $sra_id after maximum attempts"
            create_dummy_structure
            touch {params.output_directory}/$sra_id/unsuccessful.txt
            exit 0
        fi


        # Calculate the percentage of downloaded datasets
        max_spots=$(echo "$metadata" | cut -d ',' -f4)
        calculated_spots=$(echo "scale=0; $max_spots * {params.perc} / 100" | bc)


        # Download raw_fastq
        retries=0
        success=0
        if [ "$calculated_spots" -gt 0 ]; then
            while [ $retries -lt {params.retries} ] && [ $success -eq 0 ]; do
                if fastq-dump "$sra_id" -X "$calculated_spots" --split-files --outdir $output_dir; then
                    echo "{params.perc}% of dataset ($max_spots spots) successfully downloaded for: $sra_id"
                    success=1
                    touch {params.output_directory}/$sra_id/successful.txt
                    break
                else
                    echo "Attempt $(($retries + 1)) of {params.retries} failed to download {params.perc}% of dataset for: $sra_id"
                    retries=$(($retries + 1))
                    sleep {params.wait_time}
                fi
            done
            if [ $success -eq 0 ]; then
                echo "Failed to download {params.perc}% of dataset for: $sra_id after maximum attempts"
                create_dummy_structure
                touch {params.output_directory}/$sra_id/unsuccessful.txt
            fi
        else
            echo "Failed to download could not determine the total spots for: $sra_id"
            create_dummy_structure
            touch {params.output_directory}/$sra_id/unsuccessful.txt
        fi
        """

# Rule to wait until all downloads are completed like a checkpoint because errors occured cause of parallelisation
rule mark_all_downloads_complete:
    input:
        expand("{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq", output_directory=OUTPUT_DIRECTORY, sra_id=extract_sra_ids(CSV_FILE)),
        expand("{output_directory}/{sra_id}/metadata/{sra_id}.run.metadata.txt", output_directory=OUTPUT_DIRECTORY, sra_id=extract_sra_ids(CSV_FILE))
    output:
        temp("{output_directory}/complete_download.txt")
    shell:
        """
        echo "All downloads are complete." > {output}
        """


# Rule to compile lists of successful and unsuccessful downloads based on processed_download.txt markers
rule compile_download_lists:
    input:
        complete_marker = "{output_directory}/complete_download.txt",
    output:
        success_sra_list = "{output_directory}/success_sra_ids.txt",
        failed_sra_list = "{output_directory}/failed_sra_ids.txt"
    run:
        success_ids = []
        failed_ids = []
        for sra_id in extract_sra_ids(CSV_FILE):
            sra_dir_path = os.path.join(OUTPUT_DIRECTORY, sra_id)
            success_path = os.path.join(sra_dir_path, "successful.txt")
            failure_path = os.path.join(sra_dir_path, "unsuccessful.txt")
            if os.path.isfile(success_path):
                success_ids.append(sra_id)
            elif os.path.isfile(failure_path):
                failed_ids.append(sra_id)
            else:
                print(f"Neither success nor failure marker found for {sra_id}.")
        
        with open(output.success_sra_list, 'w') as f_success:
            for id in success_ids:
                f_success.write(f"{id}\n")
        
        with open(output.failed_sra_list, 'w') as f_failed:
            for id in failed_ids:
                f_failed.write(f"{id}\n")


# Rule to summarize downloads
rule summarize_downloads:
    input:
        success_sra_list = "{output_directory}/success_sra_ids.txt",
        failed_sra_list = "{output_directory}/failed_sra_ids.txt",
        # Add the dependency on the finished downloads and metadata files
        fastq_files = expand("{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq", output_directory=OUTPUT_DIRECTORY, sra_id=extract_sra_ids(CSV_FILE)),
        metadata_files = expand("{output_directory}/{sra_id}/metadata/{sra_id}.run.metadata.txt", output_directory=OUTPUT_DIRECTORY, sra_id=extract_sra_ids(CSV_FILE))
    output:
        summary = "{output_directory}/download_summary.txt"
    run:
        # Delete successful.txt out of downloads
        sra_ids = extract_sra_ids(CSV_FILE)
        for sra_id in sra_ids:
            sra_dir_path = os.path.join(OUTPUT_DIRECTORY, sra_id)
            success_file_path = os.path.join(sra_dir_path, "successful.txt")
            if os.path.exists(success_file_path):
                os.remove(success_file_path)
            

            # Cleanup directories with unsuccessful downloads
            unsuccessful_file_path = os.path.join(sra_dir_path, "unsuccessful.txt")
            if os.path.exists(unsuccessful_file_path):
                shutil.rmtree(sra_dir_path)
                print(f"Deleted {sra_dir_path} due to unsuccessful download.")


        # Read all successful and failed downloads and write it to download_summary.txt
        with open(input.success_sra_list) as f:
            success_count = sum(1 for line in f)
        with open(input.failed_sra_list) as f:
            failed_count = sum(1 for line in f)
        with open(output.summary, "w") as f:
            f.write("Download Summary\n")
            f.write("---------------------------------\n")
            f.write(f"Successful downloads: {success_count}\n")
            f.write(f"Failed downloads: {failed_count}\n")