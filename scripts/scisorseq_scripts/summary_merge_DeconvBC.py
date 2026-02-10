import os
import pandas as pd

# Function to extract values from a summary file
def extract_values(file_path):
    values = {}
    with open(file_path, 'r') as file:
        for line in file:
            key, value = line.strip().split(': ')
            values[key] = float(value)
    return values

# List all summary files in the current directory and its subdirectories
summary_files = [os.path.join(root, file) for root, dirs, files in os.walk(os.getcwd()) for file in files if file.startswith('DeconvBC_') and file.endswith('_summary')]
print(summary_files)

summary_files.sort() # Sort by filenames 

data = pd.DataFrame()

# Iterate through each summary file
for file_path in summary_files:
    values = extract_values(file_path)

    # Extract only the file name from the full path
    # file_name = "_".join(os.path.basename(file_path).split("_")[6:8])
    file_name = file_path.split("/")[-1]
    print(file_name)
    # Add the values to the DataFrame
    data[file_name] = pd.Series(values)

data.to_csv('merged_barcode_summary_table.csv', index_label='File')
