import os
import pandas as pd
from multiprocessing import Pool

# Main directories
main_dir = '/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/bowtie2_unhost/unhost_mouse'

unhost_gtf_dir = '/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP414264_mouse_gtf'
os.makedirs(unhost_gtf_dir, exist_ok=True)

# Convert input SAM files to BAM files
def sam_to_bam(srr_id):
    input_sam = os.path.join(main_dir, f"{srr_id}_mouse.sam")
    output_bam = os.path.join(main_dir, f"{srr_id}.bam")
    
    # Check if the output BAM file already exists
    if os.path.exists(output_bam):
        print(f"{output_bam} already exists")
        return
    
    os.system(f"samtools sort -@ 8 {input_sam} -o {output_bam}")

def run_stringtie_covid(srr_id):
    input_bam = os.path.join(main_dir, f"{srr_id}.bam")
    output_file = os.path.join(unhost_gtf_dir, f"{srr_id}_mouse.gtf")
    #annotation_file = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/reference_genome/UCSChg38_human.gtf"
    annotation_file = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/reference_genome/GRCm39_mouse.gff"
    
    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"{output_file} already exists")
        return

    # Run StringTie #### IMP Note: -e is needed to avoid novel transcripts ####
    os.system(f"stringtie -l {srr_id} {input_bam} -o {output_file} -G {annotation_file} -e -p 10")

def process_srr_id(srr_id):
    # Run sam_to_bam and then run_stringtie_covid
    #sam_to_bam(srr_id)
    run_stringtie_covid(srr_id)

if __name__ == "__main__":
    metadata_csv = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP414264_metadata.csv"

    # Read the metadata CSV and get the list of SRR IDs
    df = pd.read_csv(metadata_csv)
    srr_ids = df['Run'].tolist()
    
    # Calculate the number of processes based on the number of available CPU cores divided by 2
    num_processes = 8
    
    # Process each SRR ID by first converting SAM to BAM and then running StringTie
    with Pool(processes=num_processes) as pool:
        pool.map(process_srr_id, srr_ids)
