import os
import csv
from multiprocessing import Pool

def average_read_length(fastq_file):
    total_length = 0
    read_count = 0

    with open(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                # Skip header lines
                continue
            total_length += len(line.strip())
            read_count += 1

    if read_count == 0:
        return 0  # Return 0 if no reads found
    else:
        return total_length / read_count

def process_fastq_file(args):
    directory, fastq_file = args
    srr_id = fastq_file.split('_')[0]
    if '1.fastq' in fastq_file:
        paired_end = 'front'
    else:
        paired_end = 'reverse'
    avg_length = average_read_length(os.path.join(directory, fastq_file))
    return [srr_id, paired_end, avg_length]

def main():
    # Directory containing FASTQ files
    directory = "/mmfs1/scratch/deewan.bajracharya/SRR_trimmed_Man3"
    #directory = "/mmfs1/scratch/deewan.bajracharya/SRR_trial"
    # Output CSV file
    output_csv = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/RSEM/mRNA_reads_count.csv"
    
    # List all files in the directory ending with .fastq
    fastq_files = [f for f in os.listdir(directory) if f.endswith('.fastq')]

    # Create a list of arguments for each file
    args_list = [(directory, f) for f in fastq_files]

    # Create a multiprocessing pool
    with Pool(processes=6) as pool:
        # Calculate average read length for each FASTQ file
        results = pool.map(process_fastq_file, args_list)

    # Write header to CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Filename', 'Paired_End', 'Average_Read_Length'])

        # Write data to CSV file
        for result in results:
            writer.writerow(result)

    print(f"Average read lengths have been written to: {output_csv}")

if __name__ == "__main__":
    main()
