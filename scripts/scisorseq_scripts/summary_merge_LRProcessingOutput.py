import os
import pandas as pd

# Function to extract values from a summary file
def extract_values(file_path):
    values = {}
    with open(file_path, 'r') as file:
        next(file)  # Skip the first line (header)
        for line in file:
            key, value = line.strip().split("    ") # 4 spaces
            values[key] = int(value)
    return values


#LRProcessing_pass_reject_summary

# List all summary files in the current directory and its subdirectories
summary_files = [os.path.join(root, file) for root, dirs, files in os.walk(os.getcwd()) for file in files if file.startswith('LRProcessing_') and file.endswith('_summary')]
print(summary_files)
summary_files.sort() # Sort by filenames 

data = pd.DataFrame()

# Iterate through each summary file
for file_path in summary_files:
    values = extract_values(file_path)

    # Extract the filename from the full path
    file_name = file_path.split("/")[-1]
    print(file_name)

    data[file_name] = pd.Series(values)

# Define mapping for the first column to new column names
row_mapping = {
    "mapping.bestperRead.bam": "mapped_barcoded",
    "mapping.bestperRead.noRiboRNA.bam": "mapped_barcoded",
    "mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz": "spliced_mapped_barcoded",
    "newIsoforms_vs_Anno_ignoreAnno/exonStretches.gz": "spliced_mapped_barcoded"
}

# Rename index using row_mapping
data['rename'] = data.index.map(row_mapping.get)

data.to_csv('merged_LRProcessing_summary_table.csv', index_label='File')
