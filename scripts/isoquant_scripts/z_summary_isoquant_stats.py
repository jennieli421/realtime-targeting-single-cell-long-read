import os
import sys
import pandas as pd

def process_isoquant_stats_file(file_path):
    """Read a summary file and process it to extract relevant data."""
    df = pd.read_csv(file_path, sep='\t', header=0)
    print(df.head(3))
    # Extract the sample name
    print(df.columns.values)
    sample_name = df.columns.values[2]
    metric_col = df.iloc[:, 1]
    value_col = df.iloc[:, 2]
    # Create a dictionary with metrics as keys and values as values
    data = dict(zip(metric_col, value_col))
    
    return sample_name, data


def consolidate_data(input_dir):
    """Consolidate data from all summary files in the given directory."""
    all_data = {}
    samples = []
    # List all summary files in the current directory and its subdirectories
    files_to_summarize = [os.path.join(input_dir, file) for input_dir, dirs, files in os.walk(os.getcwd()) for file in files if file.startswith("isoquant_stats_summary")]
    # print(files_to_summarize)
    files_to_summarize.sort() # Sort by filenames 
    
    for file_path in files_to_summarize:
        print(file_path)
        sample_name, data = process_isoquant_stats_file(file_path)
        samples.append(sample_name)
        all_data[sample_name] = data

    df = pd.DataFrame(all_data)
    df.index.name = 'Metric'
    return df

def main(input_dir):
    output_file = 'summary_isoquant_stats.tsv'

    # Consolidate the data and save to a TSV file
    df = consolidate_data(input_dir)
    df.to_csv(output_file, sep='\t')

    print(f"Combined summary table created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python combine_summaries.py <input_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    main(input_dir)

