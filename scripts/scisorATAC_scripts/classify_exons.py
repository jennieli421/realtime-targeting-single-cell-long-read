import os
import pandas as pd

input_folder = '/home/xil4009/scratch_tilgnerlab/celltype_dpsi_P56M8VIS_P56M9VIS/'
M8 = f'{input_folder}/M8_dpsiAbove05_caseAbove10_ctlAbove10/'
M9 = f'{input_folder}/M9_dpsiAbove05_caseAbove10_ctlAbove10/'

# List files in input folder
files = [file for file in os.listdir(input_folder) if file.endswith('_dPSI_merged.tsv')]

# Iterate over each file in the input folder
for file in files:
    input_path = os.path.join(input_folder, file)
    df = pd.read_csv(input_path, sep='\t', index_col=0)
    
    M8_df = df[(abs(df['m8_dPSI']) >= 0.5) & (df['m8_case_total']>10) & (df['m8_clt_total']>10) ]
    M9_df = df[(abs(df['m9_dPSI']) >= 0.5) & (df['m9_case_total']>10) & (df['m9_clt_total']>10)]

    # Write extracted dataframe to output file
    M8_df.to_csv(f'{M8}/{"_".join(file.split("_")[:3])}_M8_dpsiAbove05_caseAbove10_ctlAbove10', sep='\t', header=False)
    M9_df.to_csv(f'{M9}/{"_".join(file.split("_")[:3])}_M9_dpsiAbove05_caseAbove10_ctlAbove10', sep='\t', header=False)

    print(f"Extracted rows from {file}")

print("Extraction completed successfully!")







####

# # Directory paths
# input_folder = '/home/xil4009/scratch_tilgnerlab/celltype_dpsi_P56M8VIS_P56M9VIS/'
# both_dpsi_above_05 = f'{input_folder}/both_dpsi_above_05/'

# # List files in input folder
# files = [file for file in os.listdir(input_folder) if file.endswith('_dPSI_merged.tsv')]

# # Iterate over each file in the input folder
# for file in files:
#     input_path = os.path.join(input_folder, file)
    
#     # Read file into pandas dataframe
#     df = pd.read_csv(input_path, sep='\t', index_col=0)
    
#     # Extract rows where both "m8_dPSI" and "m9_dPSI" satisfy criteria
#     extracted_df = df[(abs(df['m8_dPSI']) >= 0.5) & (abs(df['m9_dPSI']) >= 0.5)]
    
#     # Write extracted dataframe to output file
#     extracted_df.to_csv(f'{both_dpsi_above_05}/{file}_both_dpsi_above05', sep='\t', header=False)

#     print(f"Extracted rows from {file}")

# print("Extraction completed successfully!")
