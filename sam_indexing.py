import os
import pandas as pd
from multiprocessing import Pool

def run_samidx(args):
    srr_id, input_dir, output_dir = args
    
    # Check if BAM file exists and if exists remove BAI file
    bam_file = os.path.join(input_dir, srr_id, f"{srr_id}Aligned.sortedByCoord.out.bam")
    
    bai_file = f"{bam_file}.bai"
    if os.path.exists(bai_file):
        os.remove(bai_file)
    
    # Generate bam index file needed for sam idxstats
    sam_bai = f"samtools index {bam_file}"
    os.system(sam_bai)
        
    # Generate idxstats using samtools
    idxstats_cmd = f"samtools idxstats {bam_file} > {output_dir}/{srr_id}_idxstats.txt"
    os.system(idxstats_cmd)

if __name__ == "__main__":
    # Define input and output directories
    input_dir = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/expression_viral"
    output_dir = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/viralCopy_files"
    
    # Get list of subdirectories (SRR IDs)
    srr_ids = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

    # Create argument list for multiprocessing pool
    args_list = [(srr_id, input_dir, output_dir) for srr_id in srr_ids]

    # Run STAR command using multiprocessing pool
    with Pool(processes=8) as pool:
        pool.map(run_samidx, args_list)
