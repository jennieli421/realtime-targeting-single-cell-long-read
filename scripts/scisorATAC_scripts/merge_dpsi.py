import os
import pandas as pd

# Directory paths
output_folder = '/home/xil4009/scratch_tilgnerlab/celltype_dpsi_P56M8VIS_P56M9VIS/'
m8_folder = '/home/xil4009/store_tilgnerlab/SC_P56_M8_VIS/scisorATAC'
m9_folder = '/home/xil4009/store_tilgnerlab/SC_P56_SC_M9_VIS/scisorATAC'

# List folders in m8 and m9 directories
m8_folders = [entry.name for entry in os.scandir(m8_folder) if entry.is_dir()]
m9_folders = [entry.name for entry in os.scandir(m9_folder) if entry.is_dir()]

# Ensure that only common folders are considered
all_comparisons = set(m8_folders).union(m9_folders)

# Iterate over common folders
for folder in all_comparisons:
    folder_path_m8 = os.path.join(m8_folder, folder)
    folder_path_m9 = os.path.join(m9_folder, folder)

    # Define file paths for m8 and m9
    file_path_m8 = os.path.join(folder_path_m8, 'cases_vs_controls.counts.dPSI')
    file_path_m9 = os.path.join(folder_path_m9, 'cases_vs_controls.counts.dPSI')

    # Check if either file is empty
    if os.path.getsize(file_path_m8) == 0 or os.path.getsize(file_path_m9) == 0:
        print(f"Skipping {folder} because one of the files is empty.")
        continue

    # Read files into pandas dataframes
    df_m8 = pd.read_csv(file_path_m8, sep='\t', header=None, index_col=0)
    df_m9 = pd.read_csv(file_path_m9, sep='\t', header=None, index_col=0)

    # Join dataframes on index (first column), filling missing values with "NA"
    merged_df = df_m8.join(df_m9, how='outer', lsuffix='_m8', rsuffix='_m9').fillna('NA')

    # Add header
    merged_df.columns = ['m8_case_in', 'm8_case_ex', 'm8_clt_in', 'm8_clt_ex', 'm8_case_total', 'm8_clt_total', 'm8_dPSI',
                        'm9_case_in', 'm9_case_ex', 'm9_clt_in', 'm9_clt_ex', 'm9_case_total', 'm9_clt_total', 'm9_dPSI']

    # Save merged dataframe to a new file in the output directory
    output_file_path = os.path.join(output_folder, f"{folder}_counts_dPSI_merged.tsv")
    merged_df.to_csv(output_file_path, sep='\t', header=True)

    print(f"Merged {folder}")

print("All files merged successfully!")







# import os
# import pandas as pd

# # Directory paths
# output_folder = '/home/xil4009/scratch_tilgnerlab/P56M8VIS_P56M9VIS_layer_dpsi/'
# folder1 = f'{output_folder}/P56M8VIS_dpsi/'
# folder2 = f'{output_folder}/P56M9VIS_dpsi/'

# # cols: exon, m8_case_in, m8_case_ex, m8_clt_in, m8_clt_ex, m8_case_total, m8_clt_total, m8_dPSI, m9_case_in, m9_case_ex, m9_clt_in, m9_clt_ex, m9_case_total, m9_clt_total, m9_dPSI 

# # List files in folder1
# files = os.listdir(folder1)

# # Iterate over each file in folder1
# for file in files:
#     file_path1 = os.path.join(folder1, file)
#     file_path2 = os.path.join(folder2, file)
#     # output_path = os.path.join(output_folder, file)
    
#     # Check if either file is empty
#     if os.path.getsize(file_path1) == 0 or os.path.getsize(file_path2) == 0:
#         print(f"Skipping {file} because one of the files is empty.")
#         continue

#     # Read files into pandas dataframes
#     df1 = pd.read_csv(file_path1, sep='\t', header=None, index_col=0)
#     df2 = pd.read_csv(file_path2, sep='\t', header=None, index_col=0)
    
#     # Join dataframes on index (first column), filling missing values with "NA"
#     merged_df = df1.join(df2, how='outer', lsuffix='_1', rsuffix='_2').fillna('NA')
    
#     # Add header
#     merged_df.columns = ['m8_case_in', 'm8_case_ex', 'm8_clt_in', 'm8_clt_ex', 'm8_case_total', 'm8_clt_total', 'm8_dPSI',
#                          'm9_case_in', 'm9_case_ex', 'm9_clt_in', 'm9_clt_ex', 'm9_case_total', 'm9_clt_total', 'm9_dPSI']
    
#     merged_df.to_csv(f'{file}_merged', sep='\t', header=True)

#     print(f"Merged {file}")

# print("All files merged successfully!")
