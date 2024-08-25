# Using genes.results files created by rsem-calculate-expression --star function, read_count (w/o normalization for DESeq2) matrix is created here:

import os
import pandas as pd

# Define the input directory path
input_dir = "/path/to/RSEM_host_directory/"

# Initialize an empty dictionary to store data for each SRR ID
data = {}

# Iterate over subdirectories in the input directory
for subdir in os.listdir(input_dir):
    subdir_path = os.path.join(input_dir, subdir)
    if os.path.isdir(subdir_path):
        # Get the file path for the genes.results file
        genes_results_file = os.path.join(subdir_path, f"{subdir}.genes.results")
        if os.path.isfile(genes_results_file):
            # Read the genes.results file into a DataFrame
            df = pd.read_csv(genes_results_file, sep='\t')
            # Extract gene_id and expected_count columns
            df = df[['gene_id', 'expected_count']]
            # Rename the columns to include the SRR ID
            df.columns = ['gene_id', subdir]
            # Set the gene_id column as the index
            df.set_index('gene_id', inplace=True)
            # Store the DataFrame in the dictionary
            data[subdir] = df

# Concatenate all DataFrames into one DataFrame
count_table = pd.concat(data.values(), axis=1, sort=False)

# Write the count table to a CSV file
count_table.to_csv("/path/to/RSEM_host_directory/count_matrix_host.csv")
