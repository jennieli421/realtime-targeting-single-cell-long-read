import pandas as pd
import glob

# Define the file pattern and read all matching files into a list of DataFrames
file_pattern = "*_Ontarget_readsPerTimebin.txt"
all_files = glob.glob(file_pattern)

# Read all files into dataframes and store them in a list
dfs = []
for filename in all_files:
    df = pd.read_csv(filename, sep="\t", header=None)
    sample = df.iloc[0, -1]
    df = df.iloc[:, :-1]  # drop last column 
    df.columns = ["time_bin", "hour_range", sample]
    dfs.append(df)

merged_df = dfs[0]
for df in dfs[1:]:
    merged_df = pd.merge(merged_df, df, on=["time_bin", "hour_range"], how='inner')

order = [
		"time_bin", "hour_range", 
		"accept_bcFilter_splice_allinfo", "accept_bcFilter_splice", "accept_noBCFilter_splice", 
		"control_bcFilter_splice_allinfo","control_bcFilter_splice", "control_noBCFilter_splice"
]

merged_df = merged_df.reindex(columns=order) # Reorder columns

merged_df.to_csv("summary_allReadTypes_ontarget_readsPerTimebin", sep="\t", index=False, header=True)