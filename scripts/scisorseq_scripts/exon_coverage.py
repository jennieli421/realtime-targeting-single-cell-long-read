import os
import sys
import argparse
import gzip
import pandas as pd
import numpy as np
from traceback import print_exc


""" exons per gene, reads per exon"""
def allinfo_exonblock_to_annotation_exonblock(input_file, out_prefix, targets):
    # bedtools intersect output 
    df = pd.read_csv(input_file, sep="\t", header=None, names=["allinfo_chr","allinfo_start","allinfo_stop","allinfo_readID","allinfo_geneID","allinfo_strand","anno_chr","anno_start","anno_stop","anno_geneID","anno_type","anno_strand", "overlap"])
    if df.shape[1] != 13:
        raise ValueError(f"DataFrame must have exactly 13 columns, but it has {df.shape[1]}.")
    # make sure version number is removed 
    df['allinfo_geneID'] = df['allinfo_geneID'].str.split('.').str[0]

    df['allinfo_readID'].nunique() # 441,247
    df.shape # (1,012,636, 16)

    df['allinfo_block_size'] = df['allinfo_stop'] - df['allinfo_start']
    df['anno_block_size'] = df['anno_stop'] - df['anno_start']
    df['overlap_frac_anno'] = df['overlap'] / df['anno_block_size']

    keys = ["allinfo_readID","allinfo_geneID","anno_geneID","anno_type","overlap","allinfo_block_size",'anno_block_size', 'overlap_frac_anno']
    # df[keys]

    # Filter rows where allinfo_geneID is equal to anno_geneID 
    df1 = df.query('allinfo_geneID == anno_geneID')
    print(df1[keys])
    print("df1 shape:", df1.shape) # (1,011,657, 16)

    # Filter the rows where overlap is equal to the maximum overlap for that group
    max_overlap = df1.groupby('allinfo_readID')['overlap'].transform(max)
    df2 = df1[df1['overlap'] == max_overlap]
    print(df2[keys])
    print("df2 shape:", df2.shape)  # (802,169, 16)

    # Filterthe rows with maximum overlap_frac_annotation for that group 
    max_frac = df2.groupby('allinfo_readID')['overlap_frac_anno'].transform(max)
    df3 = df2[df2['overlap_frac_anno'] == max_frac]
    print(df3[keys])
    print("df3 shape:", df3.shape)  # (626,909, 16)

    # unique readid_exonblock different exonID 
    # df3[df3.duplicated(subset=['allinfo_readID'], keep=False)][['allinfo_start', 'allinfo_stop', 'allinfo_readID', 'anno_geneID', 'anno_type', 'anno_start', 'anno_stop']]
    # drop exonID column, dedup 
    # df4 = df3.drop(columns=['anno_exonID']).drop_duplicates()
    # print("df4 shape:", df4.shape)  

    # extract readID 
    df4 = df3.copy()
    df4['readID'] = df4.apply(lambda row: row['allinfo_readID'].split('_')[0], axis=1)
    # # Strip anything after "." in the anno_geneID column
    # df4['anno_geneID'] = df4['anno_geneID'].apply(lambda x: x.split('.')[0])
    # construct a unique exon identifier 
    df4['anno_exonblock'] = df4.apply(lambda row: f"{row['anno_chr']}_{row['anno_start']}_{row['anno_stop']}_{row['anno_geneID']}_{row['anno_strand']}", axis=1)
    df4.to_csv(f"{out_prefix}_allinfo_exonblock_to_annotation_exonblock.tsv", header=True, index=False, sep='\t')

    # df_ontarget with rows containing anno_geneID in the provided gene_ids list
    df_ontarget = df4[df4['anno_geneID'].isin(targets)]
    print("df_ontarget shape:", df_ontarget.shape)
    df_ontarget.to_csv(f"{out_prefix}_Ontarget_allinfo_exonblock_to_annotation_exonblock.tsv", header=True, index=False, sep='\t')

    anno_exonblock_count = df_ontarget.groupby('anno_geneID')['anno_exonblock'].nunique().reset_index()
    anno_exonblock_count.columns = ['anno_geneID', f'{out_prefix}_unique_anno_exonblock_count']
    # print(anno_exonblock_count.head())
    anno_exonblock_count.to_csv(f"{out_prefix}_exonPerOntargetGene.tsv", header=True, index=False, sep='\t')

    readID_count = df_ontarget.groupby('anno_exonblock')['readID'].nunique().reset_index()
    readID_count.columns = ['anno_exonblock', f'{out_prefix}_unique_readID_count']
    readID_count.to_csv(f"{out_prefix}_readPerOntargetExon.tsv", header=True, index=False, sep='\t')



""" get mapped reads, readID_geneID (not used) """
def mapped_to_annotation_geneID(input_file, out_prefix, targets):
    # bedtools intersect output 
    df = pd.read_csv(input_file, sep="\t", header=None, names=["bed_chr","bed_start","bed_stop","bed_readID","bed_score","bed_strand","anno_chr","anno_start","anno_stop","anno_geneID","anno_gene_type","anno_strand", "overlap"])
    if df.shape[1] != 13:
        raise ValueError(f"DataFrame must have exactly 13 columns, but it has {df.shape[1]}.")

    print("df shape:", df.shape)  # (1040970, 13)
    print("unique read ids: ", df['bed_readID'].nunique()) # 334,928

    # filter out non-mapped rows 
    df = df[df['overlap'] != 0 ]
    print("df shape (mapped): ", df.shape)  # (511012, 13)
    print("unique read ids: ", df['bed_readID'].nunique()) # 158608
    
    # KEEP rows of protein_coding genes
    # df = df[df['anno_gene_type'] == "protein_coding" ]
    # print("unique read ids: ", df['bed_readID'].nunique()) # 158608
    # print("df shape pc only: ", df.shape)  # (511012, 13)

    df['bed_block_size'] = df['bed_stop'] - df['bed_start']
    df['anno_block_size'] = df['anno_stop'] - df['anno_start']
    df['overlap_frac_anno'] = df['overlap'] / df['anno_block_size']

    keys = ["bed_readID", "anno_geneID", "anno_gene_type", "overlap",'overlap_frac_anno']
    print(df[keys])

    # Filter the rows where overlap is equal to the maximum overlap for that group
    max_overlap = df.groupby('bed_readID')['overlap'].transform(max)
    df1 = df[df['overlap'] == max_overlap]
    print("df1 shape (max overlap):", df1.shape)  # (161852, 16)
    print("unique read ids: ", df1['bed_readID'].nunique()) # 158608
    # print(df1[keys])

    # Filterthe rows with maximum overlap_frac_annotation for that group 
    max_frac = df1.groupby('bed_readID')['overlap_frac_anno'].transform(max)
    df2 = df1[df1['overlap_frac_anno'] == max_frac]
    print("df2 shape (max overlap_frac_annotation):", df2.shape)  # (158644, 16)
    print("unique read ids: ", df2['bed_readID'].nunique()) # 158608
    # print(df2[keys])

    # readID is df2['bed_readID'] 
    df2.to_csv(f"{out_prefix}_bed_readID_to_annotation_geneID.tsv", header=True, index=False, sep='\t')

    # Create a DataFrame with rows containing geneID in `targets``
    df_ontarget = df2[df2['anno_geneID'].isin(targets)]
    print("df_ontarget shape:", df_ontarget.shape)  # (74008, 16)
    print("unique read ids: ", df_ontarget['bed_readID'].nunique()) # 73986
    df_ontarget.to_csv(f"{out_prefix}_Ontarget_bed_readID_to_annotation_geneID.tsv", header=True, index=False, sep='\t')

    readID_count = df_ontarget.groupby('anno_geneID')['bed_readID'].nunique().reset_index()
    readID_count.columns = ['anno_geneID', f'{out_prefix}_mapped_readID_count']
    readID_count.to_csv(f"{out_prefix}_readPerOntargetGene.tsv", header=True, index=False, sep='\t')




def parse_args(sys_argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", help="Mode to run the script in", required=True, type=str)
    parser.add_argument("--input", help="timebin, hour range, allinfo", required=True, type=str)
    parser.add_argument("--out_prefix", help="output prefix", required=True, type=str)
    parser.add_argument("--targets", help="targeted gene ids", required=False, type=str)
    args = parser.parse_args(sys_argv)
    return args


def main(sys_argv):
    args = parse_args(sys_argv)

    if args.mode == "mapped_to_annotation_geneID": 
        input = args.input
        targets = pd.read_csv(args.targets, sep="\t", header=None)
        targets = targets[0].tolist()  

        print(f"Processing {input}")
        mapped_to_annotation_geneID(args.input, args.out_prefix, targets)
        print(f"====== Finished processing {input} ======")


    elif args.mode == "allinfo_to_annotation_exonblock":
        input_file = args.input
        targets = pd.read_csv(args.targets, sep="\t", header=None)
        targets = targets[0].tolist()  

        print(f"Processing {input_file}")
        allinfo_exonblock_to_annotation_exonblock(args.input, args.out_prefix, targets)
        print(f"====== Finished processing {input_file} ======")

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





# def uniq_exons_per_gene(input, out_prefix):
#     # input: ${sample}_AllInfo_Ontarget_geneID_readID_exonChain_timebin.txt
#     df = pd.read_table(input, sep='\t', header=None, 
#                        names=["gene_id","read_id", "exon_chain", "time_bin", "hour_range"])
#     # print(df.head())

#     # Split exon chains into individual exons, remove empty strings, and create a set of unique exons for each gene
#     df['exon_blocks'] = df['exon_chain'].apply(lambda x: set(filter(None, x.split(';%;'))))
#     df['exon_block_count'] = df['exon_blocks'].apply(len)
#     # print(df[['read_id', 'exon_blocks', 'exon_block_count']].head())
    
#     # Save intermediate data showing the extracted exon blocks
#     # intermediate_file = f"{out_prefix}_intermediate_exon_blocks.tsv"
#     # df[['read_id', 'gene_id', 'exon_blocks', 'exon_block_count']].to_csv(intermediate_file, sep='\t', index=False)
#     # print(f"Intermediate data saved to {intermediate_file}")
    
#     # Calculate unique exons per gene
#     print(df.head())
#     exon_counts = df.groupby('gene_id')['exon_blocks'].apply(lambda x: len(set.union(*x)) )
#     print(exon_counts.head())
#     exon_counts = exon_counts.reset_index()
#     exon_counts.columns = ['gene_id', 'unique_exon_count']
#     exon_counts['sample'] = f'{out_prefix}_splice_allinfo'
#     print(exon_counts.head())

#     output_file = f"{out_prefix}_AllInfo_unique_exons_per_target_gene.tsv"
#     exon_counts.to_csv(output_file, sep='\t', index=False)
#     print(f"Results saved to {output_file}")
