import os
import pandas as pd

# Function to extract values from a summary file
def extract_values(file_path):
    values = {}
    with open(file_path, 'r') as file:
        for line in file:
            if len(line.strip().split("\t")) == 1:
                key = "-".join(line.strip().split(" ")[:3])
                value = line.strip().split(" ")[-1]
            else: 
                key, value = line.strip().split("\t") 
            values[key] = int(value)
    return values


#LRProcessing_pass_reject_summary

# List all summary files in the current directory and its subdirectories
summary_files = [os.path.join(root, file) for root, dirs, files in os.walk(os.getcwd()) for file in files if file.endswith('UMI_filteredED2.stats.tsv')]
print(summary_files)
summary_files.sort() # Sort by filenames 

data = pd.DataFrame()

# Iterate through each summary file
for file_path in summary_files:
    values = extract_values(file_path)

    # Extract the filename from the full path
    file_name = file_path.split("/")[-1].split(".")[0]
    print(file_name)

    data[file_name] = pd.Series(values)

data.to_csv('merged_stats_summary_table.csv', index_label='File')
