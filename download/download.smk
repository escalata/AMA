import os 
configfile: "config.yaml"

# Define global variables from configuration file
CSV_FILE = config["csv_file"]
OUTPUT_DIRECTORY = config["output_directory"]
PERC = config["download_perc"]


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
        # Generate paths for .fastq files based on SRA IDs extracted from the CSV file
        expand("{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq", output_directory=OUTPUT_DIRECTORY, sra_id=extract_sra_ids(CSV_FILE)),
        expand("{output_directory}/{sra_id}/raw_fastq/{sra_id}_2.fastq", output_directory=OUTPUT_DIRECTORY, sra_id=extract_sra_ids(CSV_FILE)),
        # Path to the file listing all project IDs
        "{output_directory}/project_ids.txt".format(output_directory=OUTPUT_DIRECTORY)


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
        fastq1 = "{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq",
        fastq2 = "{output_directory}/{sra_id}/raw_fastq/{sra_id}_2.fastq",
        metadata_file = "{output_directory}/{sra_id}/metadata/{sra_id}.run.metadata.txt"
    params:
        output_directory = OUTPUT_DIRECTORY,
        csv_file = CSV_FILE,
        perc = PERC, # Percentage of data to download
        retries = 5, # Number of retry attempts
        wait_time = 5 # Wait time between retries in seconds
    shell:
        """
        sra_id={wildcards.sra_id}
        metadata_dir={params.output_directory}/$sra_id/metadata
        individual_metadata_file=$metadata_dir/$sra_id.run.metadata.txt
        mkdir -p $metadata_dir

        output_dir={params.output_directory}/$sra_id/raw_fastq
        mkdir -p $output_dir

        retries=0
        max_retries={params.retries}
        while [ $retries -lt $max_retries ]; do
            metadata=$(efetch -db sra -id "$sra_id" -format runinfo | tail -n +2)
            if [ $? -eq 0 ] && [ -n "$metadata" ]; then
                echo "$metadata" >> $individual_metadata_file
                echo "Metadata successfully downloaded for: $sra_id"
                break
            else
                echo "Attempt $(($retries + 1)) failed to download metadata for: $sra_id"
                retries=$(($retries + 1))
                sleep {params.wait_time}
            fi
        done

        if [ $retries -eq $max_retries ]; then
            echo "Failed to download metadata for: $sra_id after $max_retries attempts"
            exit 1
        fi

        retries=0
        total_spots=$(echo "$metadata" | cut -d ',' -f4)
        # Berechnet die Anzahl der herunterzuladenden Spots basierend auf dem gewünschten Prozentsatz mit Unterstützung für Fließkommazahlen
        max_spots=$(echo "scale=0; $total_spots * {params.perc} / 100" | bc)

        if [ "$max_spots" -gt 0 ]; then
            while [ $retries -lt $max_retries ]; do
                fastq-dump "$sra_id" -X "$max_spots" --split-files --outdir $output_dir
                if [ $? -eq 0 ]; then
                    echo "{params.perc}% of dataset successfully downloaded for: $sra_id"
                    break
                else
                    echo "Attempt $(($retries + 1)) failed to download {params.perc}% of dataset for: $sra_id"
                    retries=$(($retries + 1))
                    sleep {params.wait_time}
                fi
            done
        else
            echo "Could not determine the total spots for: $sra_id"
            exit 1
        fi

        if [ $retries -eq $max_retries ]; then
            echo "Failed to download {params.perc}% of dataset for: $sra_id after $max_retries attempts"
            exit 1
        fi
        """

# Rule to parse metadata files and extract project IDs
rule parse_metadata:
    input:
        # Input metadata files for all SRA IDs
        expand("{output_directory}/{sra_id}/metadata/{sra_id}.run.metadata.txt", 
               output_directory=OUTPUT_DIRECTORY, 
               sra_id=extract_sra_ids(CSV_FILE))
    output:
        # Output file listing all project IDs
        project_ids_file = "{output_directory}/project_ids.txt"
    run:
        # Extract project IDs from metadata files and write them to the output file
        project_ids = {}
        for metadata_file in input:
            # Skip processing for the header 'acc'
            if 'acc.run.metadata.txt' in metadata_file:
                continue

            with open(metadata_file, 'r') as file:
                for line in file:
                    parts = line.strip().split(',')
                    if len(parts) > 21:  # Ensure there are enough fields
                        sra_id = os.path.basename(metadata_file).split('.')[0]
                        project_id = parts[21]  # ProjectID field (22th field, python index 21)
                        project_ids[sra_id] = project_id

        # Write project IDs to the output file
        with open(output.project_ids_file, 'w') as file:
            for sra_id, project_id in project_ids.items():
                file.write(f"{sra_id}: {project_id}\n")


# Rule to signify the end of the workflow
rule finish_workflow:
    input:
        # Ensure all .fastq files have been downloaded
        expand("{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq", 
               output_directory=OUTPUT_DIRECTORY, 
               sra_id=extract_sra_ids(CSV_FILE)),
        expand("{output_directory}/{sra_id}/raw_fastq/{sra_id}_2.fastq", 
               output_directory=OUTPUT_DIRECTORY, 
               sra_id=extract_sra_ids(CSV_FILE))
    output:
    shell:
        "echo 'Partial download and metadata collection complete.'"
