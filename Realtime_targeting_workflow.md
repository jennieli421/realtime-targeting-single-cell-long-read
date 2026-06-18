# Real-time Targeting Workflow

> This document describes the complete experimental workflow for real-time adaptive sequencing using Readfish, from ONT library preparation to data preprocessing.

## Table of Contents

- [ONT Library Preparation and Sequencing](#ont-library-preparation-and-sequencing)
- [Readfish Configuration](#readfish-configuration)
- [Target Gene Sets](#target-gene-sets)
- [Reference Preparation](#reference-preparation)
- [Run Readfish](#run-readfish)
- [Basecalling](#basecalling)
- [Data Preprocessing](#data-preprocessing)
  - [Split Fastq by Condition](#split-fastq-by-condition)
  - [Barcode to Cell Type Mapping](#barcode-to-cell-type-mapping)
  - [Scisorseqr Pipeline](#scisorseqr-pipeline)
  - [Scisorseqr Statistics](#scisorseqr-statistics)
  - [UMI Correction of AllInfo](#umi-correction-of-allinfo)

---

## ONT Library Preparation and Sequencing

For all samples, 100–200 fmol cDNA processed by double-TSO removal and long-fragment selection underwent ONT library construction using the Ligation Sequencing Kit (Catalog No. SQK-LSK114, Oxford Nanopore Technologies, Oxford, UK), according to the manufacturer's protocol.

The ONT libraries were loaded onto either:
- **MinION sequencer** with MinION flow cells (Catalog No. FLO-MIN114, R10.4.1, Oxford Nanopore Technologies)
- **PromethION 2 solo sequencer** with PromethION flow cells (Catalog No. FLO-PRO114M, R10.4.1, Oxford Nanopore Technologies)

Libraries were sequenced for up to 72 h. **MinKNOW 23.07.12** (MinKNOW core 5.7.5, Guppy 7.1.4, Bream 7.7.6) was used to control the runs.

All experiments were run on a desktop with:
- Intel Xeon Silver 4214R 2.4 GHz CPUs
- Nvidia RTX A4000, 16 GB, 4DP GPUs

During runs, the MinKNOW basecall option was disabled to maximize computing resources available for real-time basecalling.

---

## Readfish Configuration

```bash
# Enable readfish basecalling
sudo chmod 774 /tmp/.guppy/5555

# Change settings in TOML config
break_reads_after_seconds = 0.4   # internal config: sequencing_PRO114_DNA_e8_2_400K.toml
max_chunks = 1                    # in TOML

# Before starting, remember to:
# 1. Reload script
# 2. Turn off basecall in MinKNOW
```

---

## Target Gene Sets

Gene sets are listed under `resources/gene_target_lists`:

| Gene Set | Size | Description |
|----------|------|-------------|
| SYNnND | 3,377 | Synaptic + neurodevelopmental genes |
| SYN | 698 | Synaptic genes only |
| LA-SYNnND VIS2 | 1,480 | Low-abundance subset for VIS-2 |
| LA-SYNnND HPC | 1,518 | Low-abundance subset for HPC |

The low-abundance (LA) genes were selected based on previous runs.

### PromethION LA-SYNnND-VIS2 target list

```bash
cd <BASE_PATH>/Exome_Enrich_additional/check_probelmatic_genes/

# 1,535 genes with < 10 reads in SYNnND-VIS2 result
target_cutoff10=<BASE_PATH>/Exome_Enrich_Major_Datasets/P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716/rf_enrichment_reference/geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10.txt

# Additional filter based on SYNnND-VIS2 result
M9junction=<BASE_PATH>/Exome_Enrich_Major_Datasets/P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716/scisorseqr/accept_bcFilter/mmalign/LRProcessingOutput/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz

# Filter out genes with problematic reads
zcat "$M9junction" | grep problematicRead | \
  awk '{split($2,a,"@;@"); print a[1]"\n"a[2]}' | sort -u | \
  awk '{split($1,a,"."); print a[1]}' > P56M9VIS_prom_probe_junction_pc_20240716_problematicRead_ids

# 54 problematic genes
wc -l P56M9VIS_prom_probe_junction_pc_20240716_problematicRead_ids

# Remove 54 genes from target list
comm -23 "$target_cutoff10" P56M9VIS_prom_probe_junction_pc_20240716_problematicRead_ids > \
  real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes
```

### MinION LA-SYNnND-HPC target list

```bash
# MinION SYNnND-HPC results
SYNnND_minion_result=<BASE_PATH>/Exome_Enrich_Major_Datasets/P56M2HPC_minion_probe_all_junction_20240125/scisorseqr/accept_bcFilter/mmalign/LRProcessingOutput/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz

# Filter out genes with problematic reads
zcat "$SYNnND_minion_result" | grep problematicRead | \
  awk '{split($2,a,"@;@"); print a[1]"\n"a[2]}' | sort -u | \
  awk '{split($1,a,"."); print a[1]}' > P56M2HPC_minion_probe_all_junction_20240125_problematicRead_ids

# 147 problematic genes
wc -l P56M2HPC_minion_probe_all_junction_20240125_problematicRead_ids

# Remove 147 genes from probe list
grep -vFf P56M2HPC_minion_probe_all_junction_20240125_problematicRead_ids "$probe" > \
  P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly

# Final MinION list: 1,518 genes
new_probe=<BASE_PATH>/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly
```

All gene sets tested are listed in **Supplementary Table 3**.

---

## Reference Preparation

```bash
# Extract protein-coding transcript sequences of desired genes
seqkit grep -nrf "$genes" "$ref" > genes.fa

# Convert FASTA to minimap2 index
minimap2 -x map-ont -d genes.mmi genes.fa
```

---

## Run Readfish

```bash
conda activate readfish_v2023

# PromethION example
toml='mouse_genelist_split_region_prometh.toml'
readfish targets --device P2S-01746-A \
  --experiment-name "P56M9VIS_synaptic" \
  --toml "$toml" \
  --log-file rf_P56M9VIS_exome_prometh_synaptic.log

# MinION example
toml='mouse_junction_probe_split_region_minion.toml'
readfish targets --device MN34448 \
  --experiment-name "P56M8VIS_lowExp_junction_pc_chunk04_fast" \
  --toml "$toml" \
  --log-file rf_P56M8VIS_exome_minion.log
```

---

## Basecalling

```bash
# Basecall on desktop
# Old PC
dorado='/home/tilgnerlab/readuntil/dorado-0.5.0-linux-x64/bin/dorado'

# New PC
dorado='/home/tilgnerlab/Downloads/dorado-0.5.3-linux-x64/bin/dorado'

# Run basecalling
"$dorado" basecaller hac pod5 --emit-fastq --min-qscore 9 > basecall_hac.fastq

"$dorado" basecaller -b 640 hac "$pod5_files" --min-qscore 9 --emit-fastq > \
  basecall_hac_20240626_1045_P2S-01746-A_PAU50570_d2c25ecb.fastq

# Upload to server for analysis
scp "$file" xil4009@aphrodite.med.cornell.edu:"${goal}"
```

---

## Data Preprocessing

### Split Fastq by Condition

All raw pod5 files were basecalled post-sequencing with high-accuracy mode of MinKNOW 23.07. Each run generated a table with read-level information. The "Condition" column indicated the region in which the read was sequenced:

- **Control**: Reads from the control condition
- **Accepted** (`stop_receiving`): Reads from the analysis condition that were accepted
- **Rejected** (`unblock`): Reads from the analysis condition that were rejected (not used for gene- or isoform-level analysis)

```bash
file="readfish.tsv"  # Output file generated by readfish

# ---- 1. Split read IDs by condition ----
mkdir decision_readID

cat "$file" | awk -F"\t" '$10=="control" {print $3}' | sort | uniq > \
  decision_readID/control_decision_readID.txt

cat "$file" | awk -F"\t" '$10=="analysis" && $9=="stop_receiving" {print $3}' | sort | uniq > \
  decision_readID/accept_decision_readID.txt

cat "$file" | awk -F"\t" '$10=="analysis" && $7=="1" && $8=="no_map" {print $3}' | sort | uniq > \
  decision_readID/noMap_1st_decision_readID.txt

grep -Fxvf decision_readID/accept_decision_readID.txt \
  decision_readID/noMap_1st_decision_readID.txt > \
  decision_readID/reject_decision_readID.txt

wc -l decision_readID/*_readID.txt

# ---- 2. Split FASTQ by condition ----
fastq="*.fastq.gz"
sbatch scripts/general/fastq_by_condition.sh "$fastq" "accept"
sbatch scripts/general/fastq_by_condition.sh "$fastq" "control"
sbatch scripts/general/fastq_by_condition.sh "$fastq" "reject"
```

### Barcode to Cell Type Mapping

Cell types for 10x barcodes were derived from matched Illumina short-read data as previously described [39]. These annotations were transferred to the long-read data via the barcode-to-cell mappings produced by **Scisorseqr (v0.1.9)** [55], which outputs a per-read record of barcode, cell type, and splice structure.

Barcode-to-cell-type files are listed under `resources/barcode2celltype`:
- HPC
- VIS1
- VIS2

### Scisorseqr Pipeline

To retain high-confidence spliced molecules, the following filters were applied uniformly to both control and enriched reads:

1. **Barcode filtering**: Barcodes identified by exact string matching against the 10x whitelist from the corresponding short-read library
2. **Alignment**: Reads aligned to the GRCm39 reference genome (GCF_000001635.27) using Minimap2 [54]; only mapped reads retained
3. **Splicing filter**: Alignments compared to GENCODE vM34 annotations; only reads spanning at least one exon-exon junction kept

All steps were performed with Scisorseqr under default parameters, except that `minTimesIsoObserve` in `InfoPerLongRead` was set to **1** to retain every spliced read in the final output `AllInfo.gz`.

```bash
# ---- 1. Run Scisorseqr ----
mkdir scisorseqr && cd $_
script="scripts/scisorseq_scripts/scisor_workflow.sh"

bcFile=<BASE_PATH>/barcode_list/P56M9VIS_barcode2celltype.txt
# bcFile=<BASE_PATH>/P56M2HPC_barcode2celltype.txt
# bcFile=<BASE_PATH>/P56M8VIS_barcode2celltype.txt

sbatch "$script" "control" "$bcFile"
sbatch "$script" "accept" "$bcFile"

# ---- 2. Summarize Scisorseqr outputs ----
conda activate readfish
python scripts/scisorseq_scripts/summary_merge_DeconvBC.py
python scripts/scisorseq_scripts/summary_merge_LRProcessingOutput.py
```

### Scisorseqr Statistics

```bash
# ---- 1. Gene coverage stats ----
cd ..

probe=Mouse.junction.probe.geneID.txt
scripts/scisorseq_scripts/scisor_stats.sh "control_bcFilter" "$probe"
scripts/scisorseq_scripts/scisor_stats.sh "accept_bcFilter" "$probe"

# ---- 2. Summarize stats across samples ----
cd analysis_target_gene_coverage

dir="."
output_file="summary_scisorseq_stats.txt"
file_paths=($(find "$dir" -type f -name 'scisorseq_stats_summary_*'))

awk 'BEGIN {OFS="\t"} {print $1 "\t" $2}' "${file_paths[0]}" > "$output_file"
for ((i=1; i<${#file_paths[@]}; i++)); do
  awk 'BEGIN {OFS="\t"} {print $2}' "${file_paths[i]}" | \
    paste -d "\t" "$output_file" - > temp_file && mv temp_file "$output_file"
done
```

### UMI Correction of AllInfo

Filtered reads were further classified by their gene assignments:
- **On-target**: Reads mapping to any target gene
- **Non-target**: All other reads

For gene-level analysis, the filtered dataset was used directly. For exon-level analyses (Figure 4C–J), unique molecular identifier (UMI)-based deduplication was performed to mitigate PCR amplification bias. For each gene-cell combination, UMIs with an edit distance less than 4 were considered technical variants of the same molecular tag and were collapsed, retaining only unique molecular events.

**Script:** `scripts/count-umis-from-all-info.py`

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

correction="scripts/count-umis-from-all-info.py"

for exp in "${exp_list[@]}"; do
  # Control pool
  input="${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo.gz"
  output="${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"
  python "$correction" "$input" "$output"

  # Accept pool
  input="${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo.gz"
  output="${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"
  python "$correction" "$input" "$output"
done
```
