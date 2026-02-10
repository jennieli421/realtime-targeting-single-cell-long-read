import os
import sys
import argparse
import gzip
import pandas as pd
from Bio import SeqIO
from traceback import print_exc
import pickle
import matplotlib.pyplot as plt
from joypy import joyplot


"""Reads the TSV file and returns a dictionary with read_id as keys and timestamp as values."""
def read2timestamp(tsv_file, pickle_file):

    if os.path.exists(pickle_file):
        print(f"Loading timestamp dictionary from {pickle_file}")
        with open(pickle_file, 'rb') as f:
            timestamp_dict = pickle.load(f)
    else:
        print(f"Generating timestamp dictionary to {pickle_file}")
        # df = pd.read_csv(tsv_file, sep='\t', dtype={'timestamp': 'float64'})
        df = pd.read_csv(tsv_file, sep='\t')
        print(df)
         # Check for any non-numeric values in the timestamp column
        if df['timestamp'].dtype != 'float64':
            print("incorrect dtype")
            df['timestamp'] = pd.to_numeric(df['timestamp'], errors='coerce')
            print(df)
        
        # Drop rows where the timestamp is NaN after conversion
        df = df.dropna(subset=['timestamp'])

        df['timestamp'] = df['timestamp'] - df['timestamp'].min()
        timestamp_dict = dict(zip(df['read_id'], df['timestamp']))

        with open(pickle_file, 'wb') as f:
            pickle.dump((timestamp_dict), f)

    return timestamp_dict

"""Processes the FASTQ file and writes the output CSV with read_id, timestamp, and read length."""
def timestamp(read_handler, output_fname, timestamp_dict):
    output_data = []

    for r in read_handler:
        read_id = str(r.id)
        seq = str(r.seq)
        timestamp = timestamp_dict.get(read_id, "NA")  # Use "NA" if timestamp is not found
        if timestamp == "NA": print("timestamp not found")
        seq_length = len(seq)
        output_data.append([read_id, timestamp, seq_length])

    output_df = pd.DataFrame(output_data, columns=["read_id", "timestamp", "seq_len"])
    print(f'Saving output to {output_fname}')
    output_df.to_csv(output_fname, index=False, sep='\t')

    # output_df['5h_time_bin'] = (output_df['timestamp'] // 18000).astype(int)  # 18000 seconds = 5 hours
    # output_df.sort_values(by='timestamp', inplace=True)  # Sort by timestamp in ascending order
    # print(output_df)
    # output_df.to_csv(f'{output_fname}_5h_timebin', index=False, sep='\t')


"""Divide readID into x-hour bins. """
def hours_time_bin(input_file, prefix, hour=5):

    df = pd.read_csv(input_file, sep='\t')
    print(df)
    # df = df.dropna(subset=['timestamp'])
    
    # Check if 'time_bin' column exists and remove it if it does
    if 'time_bin' in df.columns:
        df.drop(columns=['time_bin'], inplace=True)
    
    hour = int(hour) # ensure `hour` is an int
    bin_name = f'{hour}h_time_bin'

    df[bin_name] = (df['timestamp'] // (hour*3600)).astype(int)  # 18000 seconds = 5 hours
    df.sort_values(by='timestamp', inplace=True)  # Sort by timestamp in ascending order
    
    # Convert time_bin to actual hour numbers
    df['hour_range'] = df[bin_name].apply(lambda x: f"{x*hour}-{(x+1)*hour}")
    hour_range_order = [f"{i*hour}-{(i+1)*hour}" for i in range(df[bin_name].max() + 1)]
    # print(hour_range_order)
    df['hour_range'] = pd.Categorical(df['hour_range'], categories=hour_range_order, ordered=True)
    print(df)

    df.to_csv(f"{prefix}_readID_timestamp_readLength_{hour}hTimeBin", index=False, sep='\t') 


# # Convert time_bin to actual hour numbers
# df['hour_range'] = df['time_bin'].apply(lambda x: f"{x*5}-{(x+1)*5}")
# hour_range_order = [f"{i*5}-{(i+1)*5}" for i in range(df['time_bin'].max() + 1)]


from ChromaPalette.chroma_palette import *

"""Analyzes read lengths over time, Make ridge plot. """
def plot_ridge_read_length(input_file, experiment, sample, plotDir, x_range=[0, 4000]):
    df = pd.read_csv(input_file, sep='\t')

    hour_range_order = df['hour_range'].unique()
    print(hour_range_order)
    df['hour_range'] = pd.Categorical(df['hour_range'], categories=hour_range_order, ordered=True)

    # print(df.shape[0])
    # df = df[df["seq_len"]<15135] # filter for M2_lowEXP data 
    # print(df.shape[0])
    # print(df)
    # print("filtered!")
    

    # Set Font
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.axisbelow'] = True

    # Set Color 
    # nc = len(hour_range_order) # number of colors 
    colors = color_palette(name="WheatFields",N=14) # SoftSerenity CalmBalance(r)  SimplePastel  WheatFields Pond BlueGradient(r)
    #colors.reverse()

    # Plot ridge 
    plt.figure(figsize=(12, 10))
    fig, axes = joyplot( df, by="hour_range", column="seq_len", range_style='own',
                    overlap=1, x_range=x_range, color=colors,
                    linecolor='k', linewidth=1.2, #colormap= plt.cm.Pastel1, # plt.cm.viridis, fade=False, 
    )
    # Add grid lines  
    for ax in axes:
        ax.xaxis.grid(True, color='#EEEEEE', linestyle='-', linewidth=1.4)
    
    fig.suptitle(f"Read length over time ({sample})", fontsize=14)
    plt.title(f"{experiment}", fontsize=10)  # Subheader
    plt.xlabel("Read length (bp)")
    fig.text(0.02, 0.5, "Time range (hours)", va='center', rotation='vertical')
    plt.subplots_adjust(left=0.15, right=0.95, top=0.88, bottom=0.1)

    output_prefix = f"{plotDir}/{experiment}_{sample}" 
    print(output_prefix)
    plt.savefig(f"{output_prefix}.png", dpi=1200)
    # output_prefix = f"{plotDir}/pdf/{experiment}_{sample}" 
    plt.savefig(f"{output_prefix}.pdf")
    


def parse_args(sys_argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", help="Mode to run the script in", required=True, type=str)
    parser.add_argument("--fastq", help="fastq/fastq.gz", required=False, type=str)
    parser.add_argument("--tsv", help="readfish.tsv", required=False, type=str)
    parser.add_argument("--sample", help="output prefix", required=True, type=str)
    # parser.add_argument("--pickle", "-p", help="pickle file for timestamp dictionary", required=False, type=str)
    parser.add_argument("--read_timestamp", help="input file for read_length_overtime mode", required=False, type=str)
    parser.add_argument("--timestamp_files", nargs='+', help="input files", required=False, type=str)

    # optional args for `plot_ridge_read_length`
    parser.add_argument("--exp", help="experiment name", required=False, type=str)
    parser.add_argument("--x_range", help="x-axis range for ridge plot", required=False, type=str)
    # optional args for `hours_time_bin`
    parser.add_argument("--hour", help="hour range for time bin", required=False, type=int)
    args = parser.parse_args(sys_argv)
    return args


def main(sys_argv):
    args = parse_args(sys_argv)

    if args.mode == "read_timestamp": 
        fastq_file = args.fastq
        tsv_file = args.tsv

        if not os.path.isfile(fastq_file) or not os.path.isfile(tsv_file):
            print(f"Error: '{fastq_file}' or '{tsv_file}' does not exist.")
            sys.exit(1)

        output_fname = args.sample + ".readID_timestamp_readLength.tsv"

        pickle_file = os.path.splitext(tsv_file)[0] + ".pkl" # load pkl if exist, else create one; consistant name with readfish tsv

        print(f"Processing {fastq_file} for timestamping readID")

        # Read the TSV file and create a dictionary
        timestamp_dict = read2timestamp(tsv_file, pickle_file)

        # Open the FASTQ file and process it
        handle = fastq_file
        if fastq_file.endswith('.gz'):
            handle = gzip.open(fastq_file, "rt")
        read_handler = SeqIO.parse(handle, "fastq")

        timestamp(read_handler, output_fname, timestamp_dict)

        print(f"====== read_timestamp finished for {fastq_file} ======")

    elif args.mode == "hours_time_bin":
        input_file = args.read_timestamp
        output_prefix = args.sample
        hour = args.hour

        # os.makedirs("read_id_by_time_bin", exist_ok=True) # directory to store divided readIDs
        print(f"Processing {input_file} for assigning time bin")
        hours_time_bin(input_file, output_prefix, hour)
        print(f"====== hours_time_bin finished for {input_file} ======")

    elif args.mode == "plot_ridge_read_length":
        plotRoot = "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/all_plots/"
        plotDir = os.path.join(plotRoot, "read_length_over_time")
        if not os.path.exists(plotDir): os.makedirs(plotDir)
        # if not os.path.exists(f"{plotDir}/pdf"): os.makedirs(f"{plotDir}/pdf") # pdf in sepa

        input_file = args.read_timestamp

        # output_prefix = f"{plotDir}/{args.exp}_{args.sample}" 
        # print(output_prefix)

        plot_ridge_read_length(input_file, args.exp, args.sample, plotDir)
        # print(f"====== plot_ridge_read_length finished for {input_file} ======")

    else:
        print("Invalid mode")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



""" Combine ridge plots for three samples """
"""
def plot_ridge_3samples(input_files, experiment, plotDir, x_range=[0, 5000]):
    # Check that we have three files, experiments, and samples
    N = len(input_files)
    print(f"{N} files receved.")

    # Set up the figure for subplots
    fig, axes = plt.subplots(nrows=N, ncols=1, figsize=(12, 6*N))

    for i, (input_file, ax) in enumerate(zip(input_files, axes)):
        print(i, input_file)
        df = pd.read_csv(input_file, sep='\t')

        hour_range_order = df['hour_range'].unique()
        df['hour_range'] = pd.Categorical(df['hour_range'], categories=hour_range_order, ordered=True)

        # Set Font
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['axes.axisbelow'] = True

        # Set Color 
        colors = color_palette(name="WheatFields", N=13)
        
        # Plot ridge 
        # ax = axes[i]
        fig, axes = joyplot(
            data=df, by="hour_range", column="seq_len", range_style='own',
            overlap=1, x_range=x_range, color=colors,
            linecolor='k', linewidth=1.2, 
        )

        # Add grid lines  
        # ax.xaxis.grid(True, color='#EEEEEE', linestyle='-', linewidth=1.4)
        
        # Titles and labels
        # ax.set_title("title", fontsize=12)
        # ax.set_xlabel("Read length (bp)" if i == 2 else "")
        # ax.set_ylabel("Time range (hours)")
        
    # Set the main title for the figure
    fig.suptitle(f"Read length over time {experiment}", fontsize=16)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.3)
    
    output_prefix = f"{plotDir}/{experiment}_all" 
    print(output_prefix)
    plt.savefig(f"{output_prefix}.png", dpi=300)
    plt.savefig(f"{output_prefix}.pdf")
"""

# def read_id_by_time_bin(input_file, output_prefix):
#     df = pd.read_csv(input_file, sep='\t')

#     # Ensure 'time_bin' column is of integer type
#     df['time_bin'] = df['time_bin'].astype(int)
    
#     # Get unique time bins
#     time_bins = df['time_bin'].unique()
    
#     for bin in time_bins:
#         bin_df = df[df['time_bin'] == bin]
#         bin_df = bin_df.sort_values(by='read_id') # sort IDs 
#         output_fname = f"{output_prefix}_timebin_{bin}"
        
#         # Save the read IDs for this time bin to a file
#         bin_df[['read_id', 'time_bin']].to_csv(f'read_id_by_time_bin/{output_fname}', sep='\t', index=False, header=False)
#         print(f"Saved read IDs for time bin {bin} to {output_fname}")

# def read_id_by_time_bin(input_file, output_prefix):
#     df = pd.read_csv(input_file, sep='\t')

#     # Ensure 'time_bin' column is of integer type
#     df['time_bin'] = df['time_bin'].astype(int)
    
#     # Sort the dataframe by 'time_bin' and then by 'read_id'
#     df = df.sort_values(by=['read_id', 'time_bin'])
    
#     # Define output file name
#     output_fname = f"{output_prefix}_all_timebins.tsv"
    
#     # Save the sorted dataframe with 'read_id' and 'time_bin' columns to a single file
#     df[['read_id', 'time_bin']].to_csv(output_fname, sep='\t', index=False)
#     print(f"Saved all read IDs and time bins to {output_fname}")



    ########################### 
    # # Ensure 'seq_len' column contains only numeric values or NA
    # if not pd.api.types.is_numeric_dtype(df['seq_len']):
    #     print("Error: 'seq_len' column contains non-numeric values.")
    #     return
    
    # df['time_bin'] = (df['timestamp'] // 18000).astype(int)  # 18000 seconds = 5 hours
    # baseline_length = df[df['time_bin'] == 0]['seq_len'].mean()
    # threshold = 2 * baseline_length
    
    # result = df.groupby('time_bin').apply(
    #     lambda x: pd.Series({
    #         'n_read': len(x),
    #         'n_read_length_abnormal': (x['seq_len'] > threshold).sum(),
    #         'pct_read_length_abnormal': (x['seq_len'] > threshold).mean() * 100
    #     })
    # ).reset_index()
    
    # output_fname = output_prefix + ".read_length_overtime.tsv"
    # result.to_csv(output_fname, index=False, sep='\t')

    ###########################
    # print(df['seq_len'])
    # print(df)
    
    # # Plotting read length over time as a box plot
    # plt.figure(figsize=(12, 6))
    # df.boxplot(column='seq_len', by='time_bin', grid=False)
    # plt.title('Read Length Over Time')
    # plt.suptitle('')  # Suppress the automatic title to avoid overlap
    # plt.xlabel('Time Bin (5 hours per bin)')
    # plt.ylabel('Read Length')
    # plt.xticks(rotation=45)
    # plt.tight_layout()
    # plt.savefig(output_fname.replace('.tsv', '_boxplot.png'))
    # plt.savefig(f'{prefix}_readLength_overtime.png')
    # plt.show()