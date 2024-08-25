import os
import multiprocessing
import pandas as pd

def download_sra_run(srr_id, output_dir, check_dir):
    # Check if the file already exists
    fastq_file_path = os.path.join(check_dir, f"{srr_id}_1.fastq")
    
    if not os.path.exists(fastq_file_path):
        # If the file doesn't exist, download it
        os.system(f"fasterq-dump -O {output_dir} {srr_id}")
        print(f"{srr_id} downloaded.")
    else:
        print(f"{srr_id} already exists. Skipping download.")

def download_srr_files(output_dir, phenodata_file, check_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the phenodata file
    phenodata_df = pd.read_csv(phenodata_file)

    # Get the list of SRR IDs from phenodata_df
    srr_ids = phenodata_df['Run'].tolist()

    # Set the number of processes to use (adjust according to your CPU cores)
    num_processes = multiprocessing.cpu_count()

    # Create a multiprocessing pool
    pool = multiprocessing.Pool(processes=num_processes)

    # Map the download function to the SRR IDs using the multiprocessing pool
    pool.starmap(download_sra_run, [(srr_id, output_dir, check_dir) for srr_id in srr_ids])

    # Close the pool
    pool.close()
    pool.join()


output_dir = '/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP349864_man3'
if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
check_dir = '/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP349864_trimmed_man3'
if not os.path.exists(check_dir):
        os.makedirs(check_dir)

phenodata_file = '/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP349864_metadata.csv'

download_srr_files(output_dir, phenodata_file, check_dir)
