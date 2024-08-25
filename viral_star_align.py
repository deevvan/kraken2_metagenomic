import os
import pandas as pd
from multiprocessing import Pool

def run_star(args):
    srr_id, input_dir, output_dir, genome_dir = args
    # Create output directory if it doesn't exist
    output_srr_dir = os.path.join(output_dir, srr_id)
    os.makedirs(output_srr_dir, exist_ok=True)
    
    
    # Check if BAM file exists
    bai_file = f"{output_dir}/{srr_id}/{srr_id}Aligned.sortedByCoord.out.bam.bai"
    if os.path.exists(bai_file):
        print(f"Viral BAM file already exists for {srr_id}. Skipping STAR alignment.")
        return
    
            
    # Define input fastq files
    input_fastq1 = os.path.join(input_dir, srr_id, f"{srr_id}_un_1.fq")
    input_fastq2 = os.path.join(input_dir, srr_id, f"{srr_id}_un_2.fq")

    # Define output directory and file prefix for STAR command
    output_prefix = os.path.join(output_srr_dir, srr_id)

    # Construct STAR command
    star_cmd = f"STAR --genomeDir {genome_dir} --readFilesIn {input_fastq1} {input_fastq2} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix {output_prefix} --limitBAMsortRAM 7999715745"
    # Run STAR command
    os.system(star_cmd)
    

if __name__ == "__main__":
    # Define input and output directories
    input_dir = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/expression"
    # input_dir = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/trial"
    output_dir = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/expression_viral"
    
    genome_dir = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/reference_genome/star_cov"

    # Read the SRR IDs from the metadata CSV file
    metadata_csv = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP398274_metadata.csv"
    # metadata_csv = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/trial/SRP398274_trial.csv"
    df = pd.read_csv(metadata_csv)
    srr_ids = df["Run"].tolist()

    # Create argument list for multiprocessing pool
    args_list = [(srr_id, input_dir, output_dir, genome_dir) for srr_id in srr_ids]

    # Run STAR command using multiprocessing pool
    with Pool(processes=2) as pool:
        pool.map(run_star, args_list)