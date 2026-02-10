import os
import sys
import argparse
from traceback import print_exc


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

base_path="/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/all_plots/supplementary/"

def plot_heatmap(exp1, exp2, name1, name2, title):
    # Load the data from both files
    columns = ['gene_id', 'accept_type', 'accept_count','control_type','control_count','enrich_ratio']

    ds1 = pd.read_csv(f"{base_path}/{exp1}", sep="\t", header=None, names=columns, index_col=0)
    ds2 = pd.read_csv(f"{base_path}/{exp2}", sep="\t", header=None, names=columns, index_col=0)

    # Find overlapping gene IDs
    overlapping_ids = ds1.index.intersection(ds2.index)

    # Filter both dataframes to keep only rows with overlapping IDs
    ds1_filtered = ds1.loc[overlapping_ids].drop(columns=['accept_type', 'control_type'])
    ds2_filtered = ds2.loc[overlapping_ids].drop(columns=['accept_type', 'control_type'])

    # Merge the filtered dataframes based on the overlapping gene IDs
    data = pd.merge(ds1_filtered, ds2_filtered, left_index=True, right_index=True, suffixes=('_ds1', '_ds2'))
    print(data)


    title = title
    xlab = 'Dataset'
    ylab = 'Target genes'

    data = data[['enrich_ratio_ds1', 'enrich_ratio_ds2']]
    # Rename datasets 
    data = data.rename(columns={'enrich_ratio_ds1': name1, 'enrich_ratio_ds2': name2})
    cmap='RdYlBu_r'

    ### clustered heatmap
    clustermap = sns.clustermap(
        data, cmap=cmap, figsize=(8, 15), 
        col_cluster=True, row_cluster=True, dendrogram_ratio=(0.2, 0.02),
        cbar_pos=(-0.05, .5, .03, .4),
        yticklabels=False,
    )
    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=0, fontsize=8)
    clustermap.fig.suptitle(title)
    clustermap.ax_heatmap.set_xlabel(xlab)
    clustermap.ax_heatmap.set_ylabel(ylab)
    clustermap.savefig(f'heatmap_genewise_ratio_{title}_clustermap.png')
    plt.show()


    ### log2 transform 
    import numpy as np
    data_log2 = np.log2(data + 1)
    sorted_genes = data_log2.max(axis=1).sort_values(ascending=False).index
    data_log2 = data_log2.loc[sorted_genes]

    clustermap = sns.clustermap(
        data_log2, cmap=cmap, figsize=(8, 15), 
        col_cluster=True, row_cluster=True, dendrogram_ratio=(0.2, 0.02),
        cbar_pos=(-0.05, .5, .03, .4),
        yticklabels=False,
    )
    # disable color cbar_pos=None #cbar_kws={'orientation': 'vertical', 'shrink': 0.5}

    # Customize the color bar
    cbar = clustermap.ax_heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_title('Log2(ratio)', fontsize=15, pad=10, loc='left')
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(1)

    # Customize the clustermap title and labels
    #plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), fontsize=15)
    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=0, fontsize=8)

    clustermap.fig.suptitle(f'{title} (Log2 Transformed)', y=1.02,x=0.55,fontsize=20)
    clustermap.ax_heatmap.set_xlabel(xlab, fontsize=10)
    clustermap.ax_heatmap.set_ylabel(ylab, fontsize=10)

    clustermap.savefig(f'heatmap_genewise_ratio_{title}_log2_clustermap.png', dpi=300)
    clustermap.savefig(f'heatmap_genewise_ratio_{title}_log2_clustermap.pdf', dpi=300)



def parse_args(sys_argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", help="Mode to run the script in", required=True, type=str)
    parser.add_argument("--ds1", help="ds1", required=True, type=str)
    parser.add_argument("--ds2", help="ds2", required=True, type=str)
    parser.add_argument("--name1", help="name1", required=True, type=str)
    parser.add_argument("--name2", help="name2", required=True, type=str)
    parser.add_argument("--title", help="title", required=True, type=str)
    args = parser.parse_args(sys_argv)
    return args

def main(sys_argv):
    args = parse_args(sys_argv)

    if args.mode == "plot": 
        plot_heatmap(args.ds1, args.ds2, args.name1, args.name2, args.title)


        
if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



###### Scisorseq Output ######

"""
ds_m2_junction=/home/xil4009/store_tilgnerlab/targeted_enrich_experiments_for_thesis/P56M2HPC_transcript_enrich_all_junction_20240125/20240125_1928_MN34448_FAX66344_43d9cd4d_pass.fastq

accept_Gene2Reads=/athena/tilgnerlab/scratch/weh4002/readfish_20240125_exon_enrichment/scisorseqr/Accept_Gene2Reads_mapBestTranscriptWise_freq.txt
control_Gene2Reads=/athena/tilgnerlab/scratch/weh4002/readfish_20240125_exon_enrichment/scisorseqr/Control_Gene2Reads_mapBestTranscriptWise_freq.txt
target_geneID=/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt 
3377

ds_m2_lowExp=/home/xil4009/store_tilgnerlab/targeted_enrich_experiments_for_thesis/P56M2HPC_transcript_enrich_low_expression_junction_20240202/20240202_2142_MN34448_FAY58876_bf73ec39_pass.fastq

accept_Gene2Reads=/home/xil4009/store_tilgnerlab/targeted_enrich_experiments_for_thesis/P56M2HPC_transcript_enrich_low_expression_junction_20240202/analysis_target_gene_coverage/accept/Accept_Gene2Reads_mapBestTranscriptWise_freq.txt
control_Gene2Reads=/home/xil4009/store_tilgnerlab/targeted_enrich_experiments_for_thesis/P56M2HPC_transcript_enrich_low_expression_junction_20240202/analysis_target_gene_coverage/control/Control_Gene2Reads_mapBestTranscriptWise_freq.txt
target_geneID=/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/P56M2HPC_Missed_JunctionProbeGenes.txt
1573

Calculate Enrich Ratio: 
# 按照 gene ID 合并 accept control 两个表, 并填充缺失的 read count
cd analysis_target_gene_coverage

awk -F"\t" -v control_file="$control_Gene2Reads" -v accept_file="$accept_Gene2Reads" '
BEGIN {
    # Read the pass control file and store counts
    while ((getline < control_file) > 0) {
        clt_count[$3] = $2
    }

    # Read the accept file and store counts
    while ((getline < accept_file) > 0) {
        accept_count[$3] = $2
    }
}
{  
    gene_id = $1
    accept_reads = (gene_id in accept_count) ? accept_count[gene_id] : "0"
    clt_reads = (gene_id in clt_count) ? clt_count[gene_id] : "0"
    print gene_id "\t" "accept" "\t" accept_reads "\t" "clt" "\t" clt_reads
}
' $target_geneID > geneID_acceptReads_cltReads


# 计算 ennrich factor 
awk -F"\t" '{if ($3==0 || $5==0) {print $0 "\t" 0} else {print $0 "\t" $3/$5} }' geneID_acceptReads_cltReads > geneID_acceptReads_cltReads_enrichRatio

# 合并不同 dataset 的 enrich ratio, 留overlapping genes
ds_m2_junction=/home/xil4009/store_tilgnerlab/targeted_enrich_experiments_for_thesis/P56M2HPC_transcript_enrich_all_junction_20240125/analysis_target_gene_coverage/geneID_acceptReads_cltReads_enrichRatio
ds_m2_lowExp=/home/xil4009/store_tilgnerlab/targeted_enrich_experiments_for_thesis/P56M2HPC_transcript_enrich_low_expression_junction_20240202/analysis_target_gene_coverage1/geneID_acceptReads_cltReads_enrichRatio
# common elements = lowExp genes
comm -12 lowExp_first_col.txt junction_first_col.txt > overlap_genes_junction_lowExp
join -t $'\t' -1 1 -2 1 $ds_m2_junction $ds_m2_lowExp | awk -F"\t" 'BEGIN { print "GeneID \t P56M2HPC_allJunction_19h \t P56M2HPC_lowExp_41h"} {print $1, $6, $11}' OFS="\t" > P56M2HPC_junction_lowExp_geneID_enrichRatio

# 画图 
cd /home/xil4009/store_tilgnerlab/enrichment_factor_gene_wise
conda activate py310
python
"""

