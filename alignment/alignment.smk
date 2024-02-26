import os
from Bio import SeqIO
from datetime import datetime

configfile: "config_alignment.yaml"


# Define global variables from configuration file
csv_file = config["csv_file"]
output_directory = config["output_directory"]
query_fasta = config["query_fasta"]
# Globale variable for starttime
start_time = datetime.now()


# Function to extract SRA IDs from a given CSV file
def extract_sra_ids(csv_file):
    sra_ids = []
    with open(csv_file, 'r') as file:
        for line in file:
            sra_id = line.strip()
            # Add SRA ID to the list if it's not empty and doesn't start with "acc"
            if sra_id and not sra_id.startswith("acc"):
                sra_ids.append(sra_id)
    return sra_ids


# Calculate the average sequence length in a FASTA file
def calculate_average_length(fasta_file):
    total_length = 0
    num_sequences = 0
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        total_length += len(seq_record.seq)
        num_sequences += 1
    average_length = int(total_length / num_sequences if num_sequences else 0)
    return average_length


# Fragment sequences in a FASTA file into smaller pieces
def split_sequences(input_fasta, output_fasta, window_size):
    step_size = window_size // 2
    with open(output_fasta, 'w') as out_f:
        for seq_record in SeqIO.parse(input_fasta, 'fasta'):
            seq_length = len(seq_record.seq)
            for start in range(0, seq_length, step_size):
                end = min(start + window_size, seq_length)
                fragment = seq_record.seq[start:end]
                out_f.write(f'>{seq_record.id}_{start}-{end}\n{fragment}\n')
                if end == seq_length:
                    break


# Rule to define the final output files of the workflow
rule all:
    input:
        alignment_results = os.path.join(config["output_directory"], "alignment_results.txt"),
        #results_sra_ids_csv = os.path.join(config["output_directory"], "publication_info/results_sra_ids.csv"),
        matching_project_ids = os.path.join(config["output_directory"], "publication_info/matching_project_ids.txt"),
        full_project_ids = os.path.join(config["output_directory"], "publication_info/full_project_ids.txt"),
        publications = os.path.join(config["output_directory"], "publication_info/matching_publications.txt")


# Define the handlers
onstart:
    print("Starting SRA processing workflow.")
onsuccess:
   # Read the number of successfully downloaded SRA IDs
    with open(os.path.join(config["output_directory"], "success_sra_ids.txt"), 'r') as file:
        success_sra_ids = [line.strip() for line in file.readlines()]
    num_success_sra_ids = len(success_sra_ids)

    # Read the number of entries in the results_sra_ids.csv
    with open(os.path.join(config["output_directory"], "results_sra_ids.csv"), 'r') as file:
        # Skip the header
        next(file, None)
        sra_ids_from_results = [line.strip() for line in file]
    num_sra_ids_from_results = len(sra_ids_from_results)

    # Read the BLAST Results
    with open(os.path.join(config["output_directory"], "alignment_results.txt"), 'r') as file:
        blast_results = file.readlines()
    num_blast_results = len(blast_results)

    # Calculate the runtime
    end_time = datetime.now()
    duration = end_time - start_time
    hours, remainder = divmod(duration.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    # Print the results
    print("Processing Summary")
    print("---------------------------------")
    print(f"Number of successful SRA ID downloads: {num_success_sra_ids}")
    print(f"Number of SRA IDs with results: {num_sra_ids_from_results}")
    print(f"Number of results from the BLAST search with perc_identity= {config['blast_search']['perc_identity']} and qcov_hsp_perc= {config['blast_search']['qcov_hsp_perc']}: {num_blast_results}")
    print("---------------------------------")
    print(f"Total runtime: {hours} hours, {minutes} minutes, {seconds} seconds")    
onerror:
    print("Error encountered during SRA processing workflow.")


# Rule to preprocess raw FASTQ files using fastp
rule preprocess_fastq:
    input:
        raw_fastq1 = "{output_directory}/{sra_id}/raw_fastq/{sra_id}_1.fastq"
    output:
        clean_fastq1 = "{output_directory}/{sra_id}/fastp/{sra_id}_clean_1.fastq",
        clean_fastq2 = "{output_directory}/{sra_id}/fastp/{sra_id}_clean_2.fastq"
    shell:
        """
        # Create the directory for the output, if it does not already exist.
        mkdir -p {wildcards.output_directory}/{wildcards.sra_id}/fastp

        # Check if raw_fastq2 exists and execute fastp with one input or two inputs
        if [ -f "{wildcards.output_directory}/{wildcards.sra_id}/raw_fastq/{wildcards.sra_id}_2.fastq" ]; then
           
            # Fastp with two inputs
            fastp -i {wildcards.output_directory}/{wildcards.sra_id}/raw_fastq/{wildcards.sra_id}_1.fastq -o {output.clean_fastq1} -I {wildcards.output_directory}/{wildcards.sra_id}/raw_fastq/{wildcards.sra_id}_2.fastq -O {output.clean_fastq2} -h {wildcards.output_directory}/{wildcards.sra_id}/fastp/fastp.html -j {wildcards.output_directory}/{wildcards.sra_id}/fastp/fastp.json
        
        else
            # Fastp with one input
            fastp -i {wildcards.output_directory}/{wildcards.sra_id}/raw_fastq/{wildcards.sra_id}_1.fastq -o {output.clean_fastq1} -h {wildcards.output_directory}/{wildcards.sra_id}/fastp/fastp.html -j {wildcards.output_directory}/{wildcards.sra_id}/fastp/fastp.json

            # Create a empty dummy fastq2 for further processes
            touch {output.clean_fastq2}
        fi
        """


#Rule to convert preprocessed clean FASTQ files to FASTA format
rule convert_fastq_to_fasta:
    input:
        clean_fastq1 = "{output_directory}/{sra_id}/fastp/{sra_id}_clean_1.fastq",
        clean_fastq2 = "{output_directory}/{sra_id}/fastp/{sra_id}_clean_2.fastq"
    output:
        fasta1 = "{output_directory}/{sra_id}/blast/{sra_id}_1.fasta",
        fasta2 = "{output_directory}/{sra_id}/blast/{sra_id}_2.fasta"
    shell:
        """
        # Create the directory for the output, if it does not already exist.
        mkdir -p {output_directory}/{wildcards.sra_id}/blast

        # Converts FASTQ files into FASTA format
        seqtk seq -a {input.clean_fastq1} > {output.fasta1}
        seqtk seq -a {input.clean_fastq2} > {output.fasta2}

        # Check whether the clean_fastq files are empty and delete them
        if [ ! -s {input.clean_fastq1} ]; then
            rm {input.clean_fastq1}
        fi

         if [ ! -s {input.clean_fastq2} ]; then
            rm {input.clean_fastq2}
        fi
        """


# Rule to merge individual FASTA files into a single file
rule merge_fasta_files:
    input:
        fasta1 = "{output_directory}/{sra_id}/blast/{sra_id}_1.fasta",
        fasta2 = "{output_directory}/{sra_id}/blast/{sra_id}_2.fasta"
    output:
        merged_fasta = "{output_directory}/{sra_id}/blast/{sra_id}_merged.fasta"
    shell:
        """
        # Merges the individual FASTA files to a single file
        cat {input.fasta1} {input.fasta2} > {output.merged_fasta}

        # Check whether the fasta files are empty and delete them
        if [ ! -s {input.fasta1} ]; then
            rm {input.fasta1}
        fi
        
        if [ ! -s {input.fasta2} ]; then
            rm {input.fasta2}
        fi
        """


# Rule to merge all individual merged FASTA files into one comprehensive file
rule merge_all_fasta:
    input:
        expand("{output_directory}/{sra_id}/blast/{sra_id}_merged.fasta", output_directory=output_directory, sra_id=extract_sra_ids(csv_file))
    output:
        all_merged_fasta = "{output_directory}/blast_db/all_merged.fasta"
    shell:
        """
        mkdir -p {output_directory}/blast_db
        cat {input} > {output.all_merged_fasta}
        """


# Rule to create a BLAST database from the comprehensive merged FASTA file
rule create_combined_blast_db:
    input:
        all_merged_fasta = "{output_directory}/blast_db/all_merged.fasta"
    output:
        db_nhr = "{output_directory}/blast_db/combined_blast_db.nhr",
        db_nin = "{output_directory}/blast_db/combined_blast_db.nin",
        db_nsq = "{output_directory}/blast_db/combined_blast_db.nsq",
    shell:
        """
        makeblastdb -in {input.all_merged_fasta} -dbtype nucl -out {output_directory}/blast_db/combined_blast_db
        """


# Rule to calculate the average sequence length from the comprehensive FASTA file
rule calculate_average_length_from_db_fasta:
    input:
        db_fasta="{output_directory}/blast_db/all_merged.fasta"
    output:
        window_size_file="{output_directory}/blast_db/window_size.txt"
    run:
        window_size = calculate_average_length(input.db_fasta)
        with open(output.window_size_file, "w") as f:
            f.write(str(window_size))


# Rule to fragment the query sequences based on the average sequence length
rule fragment_query_sequences:
    input:
        query_fasta=config["query_fasta"],
        window_size_file="{output_directory}/blast_db/window_size.txt"
    output:
        fragmented_fasta="{output_directory}/blast_db/query_fragmented.fasta"
    run:
        with open(input.window_size_file) as f:
            window_size = int(f.read().strip())
        split_sequences(input.query_fasta, output.fragmented_fasta, window_size)


# Rule to perform a BLAST search of fragmented query sequences against the combined database
rule blastn_search:
    input:
        db_nhr = "{output_directory}/blast_db/combined_blast_db.nhr",
        db_nin = "{output_directory}/blast_db/combined_blast_db.nin",
        db_nsq = "{output_directory}/blast_db/combined_blast_db.nsq",
        query = "{output_directory}/blast_db/query_fragmented.fasta"
    output:
        alignment_results = "{output_directory}/alignment_results.txt"
    params:
        perc_identity = config['blast_search']['perc_identity'],
        qcov_hsp_perc = config['blast_search']['qcov_hsp_perc']
    shell:
        """
        blastn -query {input.query} -db {output_directory}/blast_db/combined_blast_db \
        -out {output.alignment_results} -outfmt 6 \
        -perc_identity {params.perc_identity} -qcov_hsp_perc {params.qcov_hsp_perc}
        """


# Parse metadata files extract project IDs and write it into publication_info folder
rule parse_metadata:
    input:
        # Input metadata files for all SRA IDs
        expand("{output_directory}/{sra_id}/metadata/{sra_id}.run.metadata.txt", 
               output_directory=output_directory,
               sra_id=extract_sra_ids(csv_file))
    output:
        # Output file listing all project IDs
        project_ids_file = "{output_directory}/publication_info/full_project_ids.txt"
    run:
        # Create the publication_info directory if it does not exist
        if not os.path.exists(os.path.dirname(output.project_ids_file)):
            os.makedirs(os.path.dirname(output.project_ids_file))

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


# Extract project IDs matching the BLAST search results
rule extract_matching_project_ids:
    input:
        alignment_results = "{output_directory}/alignment_results.txt"
    output:
        matching_project_ids = "{output_directory}/publication_info/matching_project_ids.txt"
    shell:
        """
        awk -F"\t" '{{print $2}}' {input.alignment_results} | sed 's/\..*//' | sort -u > {output_directory}/results_sra_ids.csv
        grep -F -f {output_directory}/results_sra_ids.csv {output_directory}/publication_info/full_project_ids.txt > {output.matching_project_ids}
        """


# Fetch publication IDs for projects based on matching project IDs
rule fetch_publication_ids_for_projects:
    input:
        matching_project_ids = os.path.join(config["output_directory"], "publication_info/matching_project_ids.txt")
    output:
        publications = os.path.join(config["output_directory"], "publication_info/matching_publications.txt")
    shell:
        """
        set -x
        while IFS=',' read -r sra_id project_id; do
            echo "Fetching publication ID for project: $project_id"
            efetch -db bioproject -id "$project_id" -format xml > temp.xml
            pub_id=$(cat temp.xml | xtract -pattern ProjectDescr -element Publication@id)
            if [ ! -z "$pub_id" ]; then
                echo "$project_id, $pub_id" >> {output.publications}
            else
                echo "$project_id, No publication ID found" >> {output.publications}
            fi
            rm temp.xml
        done < {input.matching_project_ids}
        """
