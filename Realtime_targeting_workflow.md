# Real-time Targeting Workflow


### ONT Library preparation and sequencing 

For all samples, 100–200 fmol cDNA processed by double-TSO removal and long-fragment selection underwent ONT library construction by using the Ligation Sequencing Kit (Catalog No. SQK-LSK114, Oxford Nanopore Technologies, Oxford, UK), according to the manufacturer’s protocol. The ONT libraries were loaded onto either MinION sequencer with MinION flow cells (Catalog No.  FLO-MIN114, R10.4.1, Oxford Nanopore Technologies), or PromethION 2 solo sequencer with PromethION flow cells (Catalog No.  FLO-PRO114M, R10.4.1, Oxford Nanopore Technologies). Libraries were sequenced for up to 72h. MinKNOW 23.07.12 (MinKNOW core 5.7.5, Guppy 7.1.4, Bream 7.7.6) was used to control the runs. All experiments were run on a desktop with Intel Xeon Silver 4214R 2.4GHz central processing units (CPUs) and Nvidia RTX A4000, 16GB, 4DP graphics cards (GPUs). During runs, MinKNOW basecall option was disabled to maximize computing resources available for real-time basecalling. 


### Configuration 

```bash
# To enable readfish basecalling 
sudo chmod 774 /tmp/.guppy/5555 

# To change settings
break_reads_after_seconds = 0.4 (internal config: sequencing_PRO114_DNA_e8_2_400K.toml)
max_chunks = 1 (in TOML)

# Before you start, remeber to 
Reload script!
turn off basecall in MinKNOW 
```

### Target gene sets
Listed under `resources/gene_target_lists`
* SYNnND (3377)
* SYN (698)
* LA-SYNnND VIS2 (1480)
* LA-SYNnND HPC (1518)


The low-abundance (LA) genes were selected based on previous runs. 

PromethION LA-SYNnND-VIS2

```bash
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/ 

# 1535 genes with < 10 reads in SYNnND-VIS2 result 
target_cutoff10=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/**P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716**/rf_enrichment_reference/geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10.txt

# Additional filter based on SYNnND-VIS2 result 
M9junction=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/**P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716**/scisorseqr/accept_bcFilter/mmalign/LRProcessingOutput/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz

# filter out genes of problematicRead  
zcat $M9junction | grep problematicRead | awk '{split($2,a,"@;@"); print a[1]"\n"a[2]}' | sort -u | awk '{split($1,a,"."); print a[1]}' > P56M9VIS_prom_probe_junction_pc_20240716_**problematicRead**_ids

# 54 /home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/P56M9VIS_prom_probe_junction_pc_20240716_**problematicRead**_ids

# remove 154 genes
comm -23 $target_cutoff10 P56M9VIS_prom_probe_junction_pc_20240716_problematicRead_ids > **real_geneID_**P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes 
```

MinION LA-SYNnND-HPC target list 

```bash
# minion SYNnND-HPC results
SYNnND_minion_result=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M2HPC_minion_probe_all_junction_20240125/scisorseqr/accept_bcFilter/mmalign/LRProcessingOutput/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz 

# filter out genes of problematicRead  
zcat $SYNnND_minion_result | grep problematicRead | awk '{split($2,a,"@;@"); print a[1]"\n"a[2]}' | sort -u | awk '{split($1,a,"."); print a[1]}' > P56M2HPC_minion_probe_all_junction_20240125_problematicRead_ids

# 147 P56M2HPC_minion_probe_all_junction_20240125_problematicRead_ids

# remove 147 genes
grep -vFf P56M2HPC_minion_probe_all_junction_20240125_problematicRead_ids $probe > P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly

# final MinION list 1518 genes 
new_probe=/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly
```

All gene sets we tested are listed in **supplementary table 3.**  



### Prepare reference for a gene set

```bash
# grep pc transcript sequences of desired genes 
seqkit grep -nrf $genes $ref > genes.fa

# convert fasta to mmi 
/minimap2 -x map-ont -d genes.mmi genes.fa
```


### Run Readfish

```bash
conda activate readfish_v2023

toml='mouse_genelist_split_region_prometh.toml'
readfish targets --device P2S-01746-A \
              --experiment-name "P56M9VIS_synaptic" \
              --toml $toml \
              --log-file rf_P56M9VIS_exome_prometh_synaptic.log

conda activate readfish_v2023
toml='mouse_junction_probe_splict_region_minion.toml'
readfish targets --device MN34448 \
              --experiment-name "P56M8VIS_lowExp_junction_pc_chunk04_fast" \
              --toml $toml \
              --log-file rf_P56M8VIS_exome_minion.log
```

### Basecalling

```bash
# Basecall on desktop
dorado='/home/tilgnerlab/readuntil/dorado-0.5.0-linux-x64/bin/dorado' # old pc
dorado='/home/tilgnerlab/Downloads/dorado-0.5.3-linux-x64/bin/dorado' # new pc
$dorado basecaller hac pod5 --emit-fastq --min-qscore 9 > basecall_hac.fastq

$dorado basecaller -b 640 hac $pod5_files --min-qscore 9 --emit-fastq > basecall_hac_20240626_1045_P2S-01746-A_PAU50570_d2c25ecb.fastq

# Upload to server for analysis 
scp $file xil4009@aphrodite.med.cornell.edu:${goal}
```

# Data preprocessing

### Split Fastq by condition 

All the raw pod5 files were basecalled post sequencing with high-accuracy mode of MinKNOW 23.07. Each run would generate a table with read-level information. Column “Condition” indicated the region in which the read was sequenced in. Reads belonged to the “control” condition were denoted as control reads. Reads from “analysis” condition were split into “accepted” (decision `stop_receiving`) or “rejected” (decision `unblock`). “Rejected” reads were not carried forward for gene- or isoform-level analysis.

```bash
file="readfish.tsv" # output file generated by readfish 

############### split reads IDs by condition ###############
mkdir decision_readID

cat $file | awk -F"\t" '$10=="control" {print $3}' | sort | uniq > decision_readID/control_decision_readID.txt
cat $file | awk -F"\t" '$10=="analysis" && $9=="stop_receiving" {print $3}' | sort | uniq > decision_readID/accept_decision_readID.txt
cat $file | awk -F"\t" '$10=="analysis" && $7=="1" && $8=="no_map" {print $3}' | sort | uniq > decision_readID/noMap_1st_decision_readID.txt
grep -Fxvf decision_readID/accept_decision_readID.txt decision_readID/noMap_1st_decision_readID.txt > decision_readID/reject_decision_readID.txt

wc -l decision_readID/*_readID.txt

############### split Fastq by condition ###############
fastq="*.fastq.gz"
sbatch scripts/general/fastq_by_condition.sh "$fastq" "accept"
sbatch scripts/general/fastq_by_condition.sh "$fastq" "control"
sbatch scripts/general/fastq_by_condition.sh "$fastq" "reject"
```

### Barcode to celltype mapping 

Cell types for 10x barcodes were derived from matched Illumina short-read data as previously described [39]. In brief, these annotations were transferred to the long‑read data via the barcode-to-cell mappings produced by Scisorseqr (v0.1.9) [55], which outputs a per-read record of barcode, cell type, and splice structure. 

Listed under `resources/barcode2celltype`
* HPC
* VIS1
* VIS2

### Scisorseq

To retain high-confidence spliced molecules, the following filters were applied uniformly to both control and enriched reads. First, barcodes were identified by exact string matching against the 10x whitelist from the corresponding short‑read library. Next, reads were aligned to the GRCm39 reference genome (GCF_000001635.27) using Minimap2 [54] and retained only if mapped. Alignments were compared to GENCODE vM34 annotations and only reads spanning at least one exon-exon junction were kept.

All steps were performed with Scisorseqr under default parameters, except that minTimesIsoObserve in InfoPerLongRead was set to 1 to retain every spliced read in the final output “Allinfo.gz”. 

```bash
############### run Scisorseqr ####################
mkdir scisorseqr && cd $_
script="scripts/scisorseq_scripts/scisor_workflow.sh"

bcFile=/home/xil4009/store_tilgnerlab/barcode_list/P56M9VIS_barcode2celltype.txt
#bcFile=/athena/tilgnerlab/scratch/weh4002/readfish_20240125_exon_enrichment/fastq_byCondition/P56M2HPC_barcode2celltype.txt
#bcFile=/athena/tilgnerlab/scratch/weh4002/readfish_20240211_exome_enrichment/P56M8VIS_barcode2celltype.txt

sbatch "$script" "control" "$bcFile"
sbatch "$script" "accept" "$bcFile"

############### summarize scisorseqr Deconv & LRProcessing summary ###############
conda activate readfish
python scripts/scisorseq_scripts/summary_merge_DeconvBC.py
python scripts/scisorseq_scripts/summary_merge_LRProcessingOutput.py
```

### Scisorseq stats

```bash
############### gene coverage stats ###############
cd ..

probe=Mouse.junction.probe.geneID.txt
scripts/scisorseq_scripts/scisor_stats.sh "control_bcFilter" "$probe"
scripts/scisorseq_scripts/scisor_stats.sh "accept_bcFilter" "$probe"

############### summarize stats across sample ###############
cd analysis_target_gene_coverage

dir="."
output_file="summary_scisorseq_stats.txt"
file_paths=($(find "$dir" -type f -name 'scisorseq_stats_summary_*'))

awk 'BEGIN {OFS="\t"} {print $1 "\t" $2}' "${file_paths[0]}" > "$output_file"
for ((i=1; i<${#file_paths[@]}; i++)); do
  awk 'BEGIN {OFS="\t"} {print $2}' "${file_paths[i]}" | paste -d "\t" "$output_file" - > temp_file && mv temp_file "$output_file"
done
```

### UMI correction of AllInfo
Filtered reads were further classified by their gene assignments: Reads mapping to any target gene were labeled “on-target”; all others were labeled “non-target”.  For gene-level analysis, the filtered dataset was used directly. For exon‑level analyses (Figure 4C–J), unique molecular identifier (UMI)‑based deduplication was performed to mitigate PCR amplification bias. For each gene-cell combination, UMIs with an edit distance less than 4 were considered technical variants of the same molecular tag and were collapsed, retaining only unique molecular events. 

Use `scripts/count-umis-from-all-info.py`

```bash
conda activate isoquant 

exp_list=(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
)

for exp in "${exp_list[@]}"; do
  correction="scripts/count-umis-from-all-info.py"

  input=${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo.gz
  python $correction $input ${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz
  
  input=${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo.gz
  python $correction $input ${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz

done
```
