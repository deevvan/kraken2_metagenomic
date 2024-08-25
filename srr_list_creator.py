import os

# Path to the directory containing the GTF files
gtf_directory = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/SRP414264_mouse_gtf"

# List all GTF files in the directory
gtf_files = [f for f in os.listdir(gtf_directory) if f.endswith('.gtf')]

# Write the sample list to a text file
with open("sample_list.txt", "w") as f:
    for gtf_file in gtf_files:
        srr_id = gtf_file.split('_')[0]
        
        f.write(f"{srr_id} {os.path.join(gtf_directory, gtf_file)}\n")
