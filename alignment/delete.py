import yaml
import os
import shutil


# Function that performs reset operations when required
def reset_operations(output_directory):
    print("Resetting operations and outputs...")
    # Delete subfolders (fastp, blast) in each SRA_ID directory
    if os.path.exists(output_directory):
        for sra_id in os.listdir(output_directory):
            sra_id_path = os.path.join(output_directory, sra_id)
            if os.path.isdir(sra_id_path):
                fastp_path = os.path.join(sra_id_path, "fastp")
                blast_path = os.path.join(sra_id_path, "blast")
                for path in [fastp_path, blast_path]:
                    if os.path.exists(path):
                        shutil.rmtree(path)
                        print(f"Deleted: {path}")
    
    # Delete result files (alignment_results.txt, results_sra_ids.csv) in the output directory
    for file_name in ["alignment_results.txt", "results_sra_ids.csv"]:
        file_path = os.path.join(output_directory, file_name)
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"Deleted: {file_path}")
    
    # Delete the publication_info folder
    publication_info_path = os.path.join(output_directory, "publication_info")
    if os.path.exists(publication_info_path):
        shutil.rmtree(publication_info_path)
        print(f"Deleted: {publication_info_path}")

    # Delete the blast_db folder
    blast_db_path = os.path.join(output_directory, "blast_db")
    if os.path.exists(blast_db_path):
        shutil.rmtree(blast_db_path)
        print(f"Deleted: {blast_db_path}")
    print("---------------------------------")
    print("Reset completed.")


# Load configuration from the config_alignment.yaml file
def load_config(filename):
    with open(filename, 'r') as f:
        return yaml.safe_load(f)

# Start
if __name__ == "__main__":
    config = load_config("config_alignment.yaml")
    output_directory = config.get("output_directory", "default_output_directory")
    reset_operations(output_directory)