
import os
import pandas as pd
import multiprocessing
import shutil

def process_srr_id(srr_id, input_dir, reference_genome_dir, output_dir):
    # Define input and output paths
    input_fastq1 = os.path.join(input_dir, f"{srr_id}_1.fastq")
    input_fastq2 = os.path.join(input_dir, f"{srr_id}_2.fastq")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_bam = f"{output_dir}/{srr_id}/{srr_id}.transcript.bam"

    # Run rsem-calculate-expression command if output files don't exist
    if not os.path.exists(output_bam):
        
        cmd = f"rsem-calculate-expression --paired-end -p 8 --output-genome-bam --star --append-names --estimate-rspd --keep-intermediate-files --time {input_fastq1} {input_fastq2} {reference_genome_dir} {output_dir}/{srr_id}"
        os.system(cmd)
        print(f"Processed SRR ID: {srr_id}")
        
        # Move specific files to new directory
        output_dir_srr = os.path.join(output_dir, srr_id)
        os.makedirs(output_dir_srr, exist_ok=True)

        # Move genes.results and transcript.bam
        shutil.move(f"{output_dir}/{srr_id}.genes.results", output_dir_srr)
        shutil.move(f"{output_dir}/{srr_id}.transcript.bam", output_dir_srr)

        # Move files from temp directory
        temp_dir = os.path.join(output_dir, f"{srr_id}.temp")
        shutil.move(f"{temp_dir}/{srr_id}_un_1.fq", output_dir_srr)
        shutil.move(f"{temp_dir}/{srr_id}_un_2.fq", output_dir_srr)
        shutil.move(f"{temp_dir}/{srr_id}Log.out", output_dir_srr)

        # Removing all intermediate files
        os.remove(os.path.join(output_dir, f"{srr_id}.genome.bam"))
        os.remove(os.path.join(output_dir, f"{srr_id}.isoforms.results"))
        os.remove(os.path.join(output_dir, f"{srr_id}.log"))
        shutil.rmtree(os.path.join(output_dir, f"{srr_id}.stat"))
        
        # Check if all files have been moved out of temp directory
        temp_files = [f"{srr_id}_un_1.fq", f"{srr_id}_un_2.fq", f"{srr_id}Log.out"]
        if all(os.path.exists(os.path.join(output_dir_srr, file)) for file in temp_files):
            # Remove temp directory
            shutil.rmtree(temp_dir)


def process_srr_ids(csv_file, input_dir, reference_genome_dir, output_dir):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Extract SRR IDs from the 'Run' column
    srr_ids = df['Run'].tolist()
    
    num_processes = 4
    
    # Process SRR IDs parallelly
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(process_srr_id, [(srr_id, input_dir, reference_genome_dir, output_dir) for srr_id in srr_ids])
    pool.close()
    pool.join()

    # Process each SRR ID sequentially
    #for srr_id in srr_ids:
        #process_srr_id(srr_id, input_dir, reference_genome_dir, output_dir)
        
        
# Define paths
csv_file = '/path/to/SRP_metadata.csv'
input_dir = '/path/to/trimmed_SRP_directory/'
reference_genome_dir = '/path/to/reference_genome_directory/star/UCSChg38_human_idx' # Star genome index for UCSChg38 genome index
output_dir = '/path/to/RSEM_host_directory/'

# Process SRR IDs
process_srr_ids(csv_file, input_dir, reference_genome_dir, output_dir)
