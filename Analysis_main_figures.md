# Analysis for Main Figures

> This document contains the complete analysis pipeline for Figures 2–4 in the manuscript.

## Table of Contents

- [Dependencies](#dependencies)
- [Figure 2A: Evaluation of size-selection protocol](#figure-2a-evaluation-of-size-selection-protocol)
- [Figure 2B–F: Cumulative on-target reads over time](#figure-2bf-cumulative-on-target-reads-over-time)
- [Figure 2G: Overall enrichment ratio heatmap](#figure-2g-overall-enrichment-ratio-heatmap)
- [Figure 3A–B: Enrichment ratio over time](#figure-3ab-enrichment-ratio-over-time)
- [Figure 4A: Gene-level enrichment ratio (stacked bars)](#figure-4a-gene-level-enrichment-ratio-stacked-bars)
- [Figure 4B: Exons per target gene (stacked bars)](#figure-4b-exons-per-target-gene-stacked-bars)
- [Figure 4C: Reads per exon (scatter plot)](#figure-4c-reads-per-exon-scatter-plot)
- [Figure 4D–H: Differential exon inclusion analysis](#figure-4dh-differential-exon-inclusion-analysis)
- [Figure 4I–J: ScisorWiz exon usage plots](#figure-4ij-scisorwiz-exon-usage-plots)

---

## Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| bedtools | v2.26.0 | Exon intersection |
| Python | 3.8.0 | Data processing |
| R | 4.3.1 | Plotting and statistics |
| Conda envs | `isoquant`, `scisorATAC`, `readfish` | Environment management |

**Custom scripts used throughout:**
- `scripts/scisorseq_scripts/exon_coverage.py`
- `scripts/general/read_timestamp.py`
- `scripts/general/R_plot_functions.R`
- `scripts/scisorseq_scripts/scisor_timebin_LRPOutput.sh`
- `scripts/scisorseq_scripts/scisor_timebin_allinfo.sh`
- `scripts/scisorseq_scripts/summary_allReadTypes_ontarget_readsPerTimebin.py`

**Key path variables:**
- `<BASE_PATH>`: Root directory for experimental data
- `<PLOT_ROOT>`: Output directory for figures
- `<ANNO_PATH>`: Directory for annotation files

---

## Figure 2A: Evaluation of size-selection protocol

### Overview
The sequencing data used preserved cDNAs from [our previous study](https://www.nature.com/articles/s41593-024-01616-4), with an additional size-selection procedure to obtain longer fragments. We compared the dataset generated from the current study (sample VIS-1 with targeting SYNnND) to the same sample (P56_M8_VIS) with probe enrichment from the previous study.

Exons covered by each read were identified by overlapping exon chains from Scisorseq `AllInfo.gz` with sorted exon blocks from GENCODE vM34 annotation using `bedtools intersect`. 100k reads from each dataset were randomly sampled, and the exon count per read was determined. The subsampling process was repeated 10 times.

### Step 1: Extract exon coordinates

**Input:**
- `AllInfo.gz` from Scisorseqr (both current and previous study)
- `anno.exon.sorted.bed` (sorted annotation exon blocks)

**Script:** `scripts/scisorseq_scripts/exon_coverage.py`

```bash
anno_exons=<ANNO_PATH>/mouse_vM34_exon_count/anno.exon.sorted.bed

# Convert AllInfo exon coordinates into sorted BED format
# Columns: chr start stop readID_exonblock geneID strand

# Previous study (BICCN)
sample="SC_P56_M8_VIS"
zcat "${biccn_allinfo}" | \
  awk '{ n=split($7,a,";%;"); for(i=2;i<=n;i++){print a[i]"\t"$2"\t"$1} }' | \
  sort -u | \
  awk 'OFS="\t" {split($1,a,"_"); print a[1], a[2], a[3], $3"_"$1, $2, a[4]}' | \
  sort -k1,1 -k2,2n | uniq > ${sample}_allinfo_exon_sorted.bed

# Current study (size-selected)
sample="VIS1_SYNnND_prom_accept_bcFilter"
zcat "${allinfo}" | \
  awk '{ n=split($7,a,";%;"); for(i=2;i<=n;i++){print a[i]"\t"$2"\t"$1} }' | \
  sort -u | \
  awk 'OFS="\t" {split($1,a,"_"); print a[1], a[2], a[3], $3"_"$1, $2, a[4]}' | \
  sort -k1,1 -k2,2n | uniq > ${sample}_allinfo_exon_sorted.bed

# Find shared exons
bedtools intersect -a ${sample}_allinfo_exon_sorted.bed \
                   -b "${anno_exons}" -sorted -wao | \
  sort | uniq > ${sample}_allinfo_exon_intersect_annotation

# Process intersection output
# Output: 1) all exons; 2) exons of on-target genes
targets=resources/gene_target_lists/Mouse.junction.probe.geneID.txt
script=scripts/scisorseq_scripts/exon_coverage.py

python "${script}" \
  --mode "allinfo_to_annotation_exonblock" \
  --input ${sample}_allinfo_exon_intersect_annotation \
  --out_prefix ${sample} \
  --targets "${targets}"
```

### Step 2: Subsample 10 times

```python
import pandas as pd
import numpy as np

# Process each dataset separately
# Dataset 1: Previous study
input_file = "SC_P56_M8_VIS_Ontarget_allinfo_exonblock_to_annotation_exonblock.tsv"
output_prefix = "SC_P56_M8_VIS"

# Dataset 2: Current study (size-selected)
input_file = "VIS1_SYNnND_prom_accept_bcFilter_Ontarget_allinfo_exonblock_to_annotation_exonblock.tsv"
output_prefix = "VIS1_SYNnND_prom_accept_bcFilter"

df = pd.read_csv(input_file, sep="\t")
unique_readIDs = df['readID'].unique()

for i in range(1, 11):
    # Randomly sample 0.1 million unique readIDs
    sampled_readID = np.random.choice(unique_readIDs, size=100000, replace=False)
    filtered_df = df[df['readID'].isin(sampled_readID)]

    # Count exons per read
    exon_per_read = filtered_df.groupby('readID')['allinfo_readID'].nunique().reset_index()
    exon_per_read.columns = ['read_id', 'exon_per_read']
    exon_per_read['sample'] = output_prefix
    exon_per_read.to_csv(
        f"{output_prefix}_exonPerRead_sub100k_rep{i}.tsv",
        header=True, index=False, sep='\t'
    )
```

### Step 3: Plot Figure 2A

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/size_selection_cutoff1/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

base_path <- "<BASE_PATH>/size_selection_compare_M8_cutoff1"

# ---- 1. Data processing ----
all_data <- list()

for (i in 1:10) {
  accept_df <- read.table(
    glue("{base_path}/VIS1_SYNnND_prom_accept_bcFilter_exonPerRead_sub100k_rep{i}.tsv"),
    sep = '\t', header = TRUE
  )
  control_df <- read.table(
    glue("{base_path}/SC_P56_M8_VIS_exonPerRead_sub100k_rep{i}.tsv"),
    sep = '\t', header = TRUE
  )
  label <- paste0("Rep", i)

  accept_df$Type <- "Size-selected"
  control_df$Type <- "Non-size-selected"
  accept_df$Label <- label
  control_df$Label <- label

  all_data[[i]] <- rbind(accept_df, control_df)
}

final_df <- do.call(rbind, all_data)
final_df$Label <- factor(final_df$Label, levels = paste0("Rep", 1:10))
final_df$Type <- factor(final_df$Type, levels = c("Size-selected", "Non-size-selected"))

# Colors
colors <- c("#efba7a", "#f1776b")

# ---- 2. Wilcoxon test ----
p_data <- final_df %>%
  group_by(Label) %>%
  summarise(
    p_val = wilcox.test(exon_per_read ~ Type, paired = TRUE)$p.value,
    max_y = 22,
    .groups = 'drop'
  ) %>%
  mutate(
    sig_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    line_y = max_y + 0.1,
    y_pos = line_y + 0.2,
    x_left = as.numeric(Label) - 0.2,
    x_right = as.numeric(Label) + 0.2
  )

# ---- 3. Plot ----
p <- ggplot(final_df, aes(x = Label, y = exon_per_read, fill = Type, color = Type)) +
  geom_boxplot(
    width = 0.5, position = position_dodge(width = 0.6),
    size = 0.5, alpha = 0.65, outlier.shape = NA
  ) +
  geom_segment(
    data = p_data,
    aes(x = x_left, xend = x_right, y = line_y, yend = line_y),
    inherit.aes = FALSE,
    color = "black", linewidth = 0.5
  ) +
  geom_text(
    data = p_data,
    aes(x = Label, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    size = 2.5, color = "black", vjust = 0
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(
    title = "Exon count per read in cDNA libraries (sample VIS-1)",
    x = NULL,
    y = "Number of exons per read"
  ) +
  ylim(0, 25) +
  theme_jennie(plot_title_size = 9, base_size = 9)

# ---- 4. Save ----
save_pdf(p, glue("{plotDir}/size_selection_exon_per_read"), 7, 2.5)
```

---

## Figure 2B–F: Cumulative on-target reads over time

### Overview
Reads in the enriched and control pools were grouped into time bins based on sequencing timestamps. For each bin, the number of on-target or non-target reads was calculated.

**Formulas:**
- **Overall enrichment ratio** (Figure 2G):  
  `total on-target spliced reads in enriched / total on-target spliced reads in control`
- **Time-bin enrichment ratio** (Figure 3A–B):  
  `on-target spliced reads in time bin t (enriched) / on-target spliced reads in time bin t (control)`
- **Gene-level enrichment ratio** (Figure 4A):  
  `on-target spliced reads assigned to gene A (enriched) / on-target spliced reads assigned to gene A (control)`

### Step 1: Split read IDs by time bins

**Inputs:**
- `fastq.gz` (basecalled reads)
- `readfish.tsv` (readfish output with timestamps)

**Script:** `scripts/general/read_timestamp.py`

```bash
cd <ROOT_DIR>
mkdir timestamp && cd timestamp

conda activate isoquant
script=scripts/general/read_timestamp.py
file=readfish.tsv

# Accept pool
python "${script}" --mode hours_time_bin \
  --read_timestamp accept.readID_timestamp_readLength.tsv \
  --sample accept --hour 1

# Control pool
python "${script}" --mode hours_time_bin \
  --read_timestamp control.readID_timestamp_readLength.tsv \
  --sample control --hour 1

# Reject pool (optional)
python "${script}" --mode read_timestamp \
  --fastq ../fastq_by_condition/reject/reject.fastq.gz \
  --tsv "${file}" --sample reject

python "${script}" --mode hours_time_bin \
  --read_timestamp reject.readID_timestamp_readLength.tsv \
  --sample reject --hour 1
```

### Step 2: Process time-binned AllInfo

**Scripts:**
- `scripts/scisorseq_scripts/scisor_timebin_LRPOutput.sh`
- `scripts/scisorseq_scripts/scisor_timebin_allinfo.sh`
- `scripts/scisorseq_scripts/summary_allReadTypes_ontarget_readsPerTimebin.py`

```bash
x="1"  # Time bin = 1 hour

mkdir allinfo_timebin && cd $_

# Skip header, get sorted file with readID - timebin
cat ../accept_readID_timestamp_readLength_${x}hTimeBin | \
  tail -n +1 | cut -f1,4,5 | sort -u > accept_readid_${x}htimebin_sorted
cat ../control_readID_timestamp_readLength_${x}hTimeBin | \
  tail -n +1 | cut -f1,4,5 | sort -u > control_readid_${x}htimebin_sorted

# Experiment-to-target mappings
declare -A pairs
pairs=(
  ["P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"]="<ANNO_PATH>/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"]="<ANNO_PATH>/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"]="<BASE_PATH>/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes"
  ["P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"]="<BASE_PATH>/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly"
)

conda activate isoquant

for folder in "${!pairs[@]}"; do
    cd "${folder}/timestamp/allinfo_timebin/"
    probe=${pairs[$folder]}
    echo "folder = ${folder}"
    echo "target list = ${probe}"

    # LRP processing
    scripts/scisorseq_scripts/scisor_timebin_LRPOutput.sh "${probe}" "${x}"

    # AllInfo processing
    wc -l *_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt
    scripts/scisorseq_scripts/scisor_timebin_allinfo.sh "${probe}" "${x}"
    wc -l *_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt

    # Summarize
    python scripts/scisorseq_scripts/summary_allReadTypes_ontarget_readsPerTimebin.py

    cd ../../../
done

# Check accept total on-target read count
awk '{sum += $3} END {print sum}'
awk '{sum += $6} END {print sum}'
```

### Step 3: Plot cumulative reads (Figure 2B–F)

**Function:** `ontarget_cumulative_count_accept_control()` (from `scripts/general/R_plot_functions.R`)

```r
conda activate scisorATAC

source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/ontarget_read_count_over_time_allinfo_cutoff1/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

base_path <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/"

exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList",
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1",
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1",
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
)

names <- c(
  "SYNnND VIS-1 (P)",
  "LA-SYNnND VIS-2 (P)",
  "LA-SYNnND HPC (M)",
  "SYNnND VIS-2 (P)",
  "SYNnND HPC (M)",
  "SYN VIS-1 (P)",
  "SYN VIS-2 (P)"
)

plot_list <- list()
for (i in seq_along(exp_list)) {
  p <- ontarget_cumulative_count_accept_control(exp_list[i], names[i])
  plot_list[[i]] <- p
}

p <- wrap_plots(plot_list, ncol = 1)
prefix <- glue("{plotDir}/stacked_bar_pie_chart_all_0326")
save_png(p, prefix, 6.5, 3.5 * length(exp_list))
save_pdf(p, prefix, 6.5, 3.5 * length(exp_list))
```

---

## Figure 2G: Overall enrichment ratio heatmap

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/enrich_ratio_heatmap_allinfo_cutoff1/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

data <- data.frame(
  Gene_set = c("LA-SYNnND", "SYNnND", "SYN"),
  VIS_1 = c(NA, 1.43, 1.82),
  VIS_2 = c(1.39, 1.82, 1.55),
  HPC = c(1.89, 1.57, NA)
)

library(reshape2)
data_long <- melt(data, id.vars = "Gene_set",
                  variable.name = "Sample", value.name = "Value")

data_long$Gene_set <- factor(data_long$Gene_set,
                             levels = c("SYN", "SYNnND", "LA-SYNnND"))

tile_line_width <- 3.5
tile_line_color <- "transparent"

p <- ggplot(data_long, aes(x = Sample, y = Gene_set, fill = Value)) +
  geom_tile(color = tile_line_color, linewidth = tile_line_width) +
  geom_text(
    aes(label = ifelse(is.na(Value), "", sprintf("%.2f", Value))),
    color = "black", size = 5.5, family = "ArialMT"
  ) +
  scale_fill_gradient(
    low = "#e8f7f4",
    high = "#91c8d6",
    na.value = "#F9F9F9",
    name = "Enrichment ratio",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1.2,
      barheight = 10,
      direction = "vertical"
    )
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_jennie(base_size = 14) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    axis.line = element_blank()
  ) +
  coord_fixed()

title <- glue("{plotDir}/Overall_enrich_ratio_3x3")
save_png(p, title, 7, 5)
save_pdf(p, title, 7, 5)
```

---

## Figure 3A–B: Enrichment ratio over time

**Input:** `summary_allReadTypes_ontarget_readsPerTimebin`

```r
conda activate scisorATAC

source("scripts/general/R_plot_functions.R")
base_path <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/"

plotDir <- glue("{<PLOT_ROOT>}/enrich_ratio_over_time_allinfo_cutoff1_1htimebin")

# ---- Panel 1: SYNnND ----
synnnd_exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1",
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"
)

synnnd_names <- c(
  "SYNnND VIS-1 (P)",
  "SYNnND VIS-2 (P)",
  "SYNnND HPC (M)"
)

name2color <- list(
  "SYNnND VIS-1 (P)" = "#1f77b4",  # blue M8
  "SYNnND VIS-2 (P)" = "#ff7f0e",  # orange M9
  "SYNnND HPC (M)" = "#3cb44b"     # green M2
)
name2linetype <- c(
  "SYNnND VIS-1 (P)" = "solid",
  "SYNnND VIS-2 (P)" = "solid",
  "SYNnND HPC (M)" = "solid"
)

df <- format_enrich_ratio_timebin(synnnd_exp_list, synnnd_names)
df <- df %>%
  filter(as.integer(as.character(hour_start)) <= 35) %>%
  mutate(hour_start = as.numeric(as.character(hour_start)))

ymax <- 3
x_breaks <- seq(0, 35, by = 2)

p1 <- ggplot(df, aes(x = hour_start, y = allinfo_ratio,
                     color = name, linetype = name, group = name)) +
  geom_smooth(se = FALSE, size = 0.7) +
  labs(
    title = NULL,
    x = "Sequencing time (hours)",
    y = "Ratio (enriched/control)"
  ) +
  ylim(0, ymax) +
  scale_color_manual(values = name2color) +
  scale_linetype_manual(values = name2linetype) +
  scale_x_continuous(breaks = x_breaks) +
  theme_jennie(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.width = unit(2.5, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(ncol = 2),
    linetype = guide_legend(ncol = 2)
  )

# ---- Panel 2: SYN, LA ----
LA_exp_list <- c(
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1",
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
)

LA_names <- c(
  "SYN VIS-1 (P)",
  "SYN VIS-2 (P)",
  "LA-SYNnND HPC (M)",
  "LA-SYNnND VIS-2 (P)"
)

LA_name2color <- list(
  "SYN VIS-1 (P)" = "#1f77b4",
  "SYN VIS-2 (P)" = "#ff7f0e",
  "LA-SYNnND VIS-2 (P)" = "#ff7f0e",
  "LA-SYNnND HPC (M)" = "#3cb44b"
)
LA_name2linetype <- c(
  "SYN VIS-1 (P)" = "dashed",
  "SYN VIS-2 (P)" = "dashed",
  "LA-SYNnND HPC (M)" = "solid",
  "LA-SYNnND VIS-2 (P)" = "solid"
)

df <- format_enrich_ratio_timebin(LA_exp_list, LA_names)
df <- df %>%
  filter(as.integer(as.character(hour_start)) <= 35) %>%
  mutate(hour_start = as.numeric(as.character(hour_start)))

p2 <- ggplot(df, aes(x = hour_start, y = allinfo_ratio,
                     color = name, linetype = name, group = name)) +
  geom_smooth(se = FALSE, size = 0.7) +
  labs(
    title = NULL,
    x = "Sequencing time (hours)",
    y = "Ratio (enriched/control)"
  ) +
  ylim(0, ymax) +
  scale_color_manual(values = LA_name2color) +
  scale_linetype_manual(values = LA_name2linetype) +
  scale_x_continuous(breaks = x_breaks) +
  theme_jennie(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.width = unit(2.5, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(ncol = 2),
    linetype = guide_legend(ncol = 2)
  )

# ---- Combine panels ----
combined_plot <- (p1 / p2)
title <- ggdraw() + draw_label("Enrich ratio over time", fontface = 'bold', size = 14)
final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.05, 1))

width <- 7
height <- 4
title <- "combined_enrich_ratio_over_time_2panels"
prefix <- glue("{plotDir}/{title}")
save_pdf(final_plot, prefix, width, height * 2)
save_png(final_plot, prefix, width, height * 2)
```

### Summary statistics (mean ± SD)

```r
exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1",
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1",
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1",
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
)

names <- c(
  "SYNnND VIS-1 (P)",
  "SYNnND VIS-2 (P)",
  "SYNnND HPC (M)",
  "SYN VIS-1 (P)",
  "SYN VIS-2 (P)",
  "LA-SYNnND HPC (M)",
  "LA-SYNnND VIS-2 (P)"
)

df <- format_enrich_ratio_timebin(exp_list, names)
df <- df %>%
  filter(as.integer(as.character(hour_start)) <= 35) %>%
  mutate(hour_start = as.numeric(as.character(hour_start)))

# Get mean, SD, max, and min for each experiment
for (exp in unique(df$name)) {
  print(exp)
  tmp <- df[df$name == exp, ] %>% select(allinfo_ratio)
  print(glue("mean ± sd: {round(mean(tmp$allinfo_ratio), 2)} ± {round(sd(tmp$allinfo_ratio), 2)}"))

  tmp_max <- df[df$name == exp, ] %>% select(time_bin, allinfo_ratio) %>% slice_max(allinfo_ratio)
  print(glue("max: allinfo_ratio = {round(tmp_max$allinfo_ratio, 2)}, time_bin = {tmp_max$time_bin}"))

  tmp_min <- df[df$name == exp, ] %>% select(time_bin, allinfo_ratio) %>% slice_min(allinfo_ratio)
  print(glue("min: allinfo_ratio = {round(tmp_min$allinfo_ratio, 2)}, time_bin = {tmp_min$time_bin}"))
}
```

---

## Figure 4A: Gene-level enrichment ratio (stacked bars)

### Step 1: Create gene-to-read count files

**Input:** `AllInfo`

```bash
exp_list=(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"
)

for exp in "${exp_list[@]}"; do
  # Accept pool
  cd "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/analysis_target_gene_coverage/accept_bcFilter/"
  sample="accept_bcFilter"
  cat AllInfo_Ontarget_Gene2reads.txt | cut -f1 | sort | uniq -c | \
    tr -s " " | awk -v sample="$sample" -F" " '{print sample "\t" $1 "\t" $2}' > \
    "${sample}_AllInfo_Ontarget_Gene2reads_freq.txt"

  # Control pool
  cd "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/analysis_target_gene_coverage/control_bcFilter/"
  sample="control_bcFilter"
  cat AllInfo_Ontarget_Gene2reads.txt | cut -f1 | sort | uniq -c | \
    tr -s " " | awk -v sample="$sample" -F" " '{print sample "\t" $1 "\t" $2}' > \
    "${sample}_AllInfo_Ontarget_Gene2reads_freq.txt"
done
```

### Step 2: Calculate enrichment ratios per gene

```bash
cd <BASE_PATH>/Exome_Enrich_Major_Datasets/genewise_enrich_ratio_tables
base_path="<BASE_PATH>/Exome_Enrich_Major_Datasets"

declare -A pairs
pairs=(
  ["P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"]="<ANNO_PATH>/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"]="<ANNO_PATH>/mouse_probe/synaptic"
  ["P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"]="<BASE_PATH>/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes"
  ["P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"]="<BASE_PATH>/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly"
)

for exp in "${!pairs[@]}"; do
  targets=${pairs[$exp]}

  accept_Gene2Reads="${base_path}/${exp}/analysis_target_gene_coverage/accept_bcFilter/accept_bcFilter_AllInfo_Ontarget_Gene2reads_freq.txt"
  control_Gene2Reads="${base_path}/${exp}/analysis_target_gene_coverage/control_bcFilter/control_bcFilter_AllInfo_Ontarget_Gene2reads_freq.txt"

  # Merge by gene ID, fill missing read counts with 0
  awk -F"\t" -v control_file="$control_Gene2Reads" -v accept_file="$accept_Gene2Reads" '
    BEGIN {
      while ((getline < control_file) > 0) { clt_count[$3] = $2 }
      while ((getline < accept_file) > 0) { accept_count[$3] = $2 }
    }
    {
      gene_id = $1
      accept_reads = (gene_id in accept_count) ? accept_count[gene_id] : "0"
      clt_reads = (gene_id in clt_count) ? clt_count[gene_id] : "0"
      print gene_id "\t" "accept" "\t" accept_reads "\t" "clt" "\t" clt_reads
    }
  ' "$targets" | awk -F"\t" '
    { if ($3 == 0 || $5 == 0) { print $0 "\t" 0 } else { print $0 "\t" $3 / $5 } }
  ' > "geneID_acceptReads_cltReads_enrichRatio_${exp}"
done
```

### Step 3: Plot stacked bars

**Function:** `ratio_per_ontarget_gene_bars()` (from `scripts/general/R_plot_functions.R`)

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/exon_coverage_umicorrect")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

base_path <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/"

exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1",
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList",
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList",
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
)

names <- c(
  "SYNnND VIS-1 (P)",
  "SYNnND VIS-2 (P)",
  "SYNnND HPC (M)",
  "LA-SYNnND VIS-2 (P)",
  "LA-SYNnND HPC (M)",
  "SYN VIS-1 (P)",
  "SYN VIS-2 (P)"
)

p <- ratio_per_ontarget_gene_bars(exp_list, names)

title <- "genewise_enrich_ratio_all"
prefix <- glue("{plotDir}/{title}")
save_pdf(p, prefix, 8, 5)
save_png(p, prefix, 8, 5)
```

---

## Figure 4B: Exons per target gene (stacked bars)

For each target gene, we compared the number of unique exons detected in the enriched versus control pool and assigned the gene to one of three categories: (1) more exons in enriched, (2) more in control, or (3) equal numbers.

**Function:** `exon_per_ontarget_gene_bars()` (from `scripts/general/R_plot_functions.R`)

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/exon_coverage_umicorrect")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

base_path <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/"

exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1",
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList",
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList",
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
)

names <- c(
  "SYNnND VIS-1 (P)",
  "SYNnND VIS-2 (P)",
  "SYNnND HPC (M)",
  "LA-SYNnND VIS-2 (P)",
  "LA-SYNnND HPC (M)",
  "SYN VIS-1 (P)",
  "SYN VIS-2 (P)"
)

p <- exon_per_ontarget_gene_bars(exp_list, names)

title <- "exon_per_target_gene_all"
prefix <- glue("{plotDir}/{title}")
save_png(p, prefix, 8, 5)
save_pdf(p, prefix, 8, 5)
```

---

## Figure 4C: Reads per exon (scatter plot)

### Overview
Exonic coordinates for spliced, on-target reads were extracted from the `AllInfo` file and intersected with GENCODE vM34 exon annotations using `bedtools intersect` (v2.26.0), yielding per-exon read counts.

**Input:** UMI-corrected AllInfo (`AllInfo_IncompleteReads.filtered.corrected.gz`)

**Script:** `scripts/scisorseq_scripts/exon_coverage.py`

### Step 1: bedtools intersect

```bash
cd <BASE_PATH>/Exome_Enrich_Major_Datasets

exp_list=(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"
)

anno_exons=<ANNO_PATH>/mouse_vM34_exon_count/anno.exon.sorted.bed

for folder in "${exp_list[@]}"; do
  mkdir -p "${folder}/exon_coverage/"
  cd "${folder}/exon_coverage/"

  samples=("accept_bcFilter" "control_bcFilter")
  for sample in "${samples[@]}"; do
    allinfo="../scisorseqr/${sample}/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"

    # Convert AllInfo exon coords to sorted BED
    zcat "${allinfo}" | \
      awk '{ n=split($7,a,";%;"); for(i=2;i<=n;i++){print a[i]"\t"$2"\t"$1} }' | \
      sort -u | \
      awk 'OFS="\t" {split($1,a,"_"); print a[1], a[2], a[3], $3"_"$1, $2, a[4]}' | \
      sort -k1,1 -k2,2n | uniq > ${sample}_allinfo_exon_sorted.bed

    bedtools intersect -a ${sample}_allinfo_exon_sorted.bed \
                       -b "${anno_exons}" -sorted -wao | \
      sort | uniq > ${sample}_allinfo_exon_intersect_annotation
  done

  cd ../../
done
```

### Step 2: Process intersections

```bash
declare -A pairs
pairs=(
  ["P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"]="<ANNO_PATH>/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"]="<ANNO_PATH>/mouse_probe/synaptic"
  ["P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"]="<BASE_PATH>/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes"
  ["P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"]="<ANNO_PATH>/Mouse.junction.probe.geneID.txt"
  ["P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"]="<BASE_PATH>/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly"
)

for folder in "${!pairs[@]}"; do
  cd "${folder}/exon_coverage/"

  targets=${pairs[$folder]}
  echo "folder = ${folder}"
  echo "target list = ${targets}"

  script=scripts/scisorseq_scripts/exon_coverage.py
  samples=("accept_bcFilter" "control_bcFilter")

  for sample in "${samples[@]}"; do
    input=${sample}_allinfo_exon_intersect_annotation
    python "${script}" \
      --mode "allinfo_to_annotation_exonblock" \
      --input "${input}" \
      --out_prefix "${sample}" \
      --targets "${targets}"
  done

  cd ../../
done
```

### Step 3: Scatter plot

**Inputs:**
- `control_bcFilter_exonPerOntargetGene.tsv`
- `accept_bcFilter_exonPerOntargetGene.tsv`
- `control_bcFilter_readPerOntargetExon.tsv`
- `accept_bcFilter_readPerOntargetExon.tsv`

**Function:** `exon_coverage_scatter_plot()` (from `scripts/general/R_plot_functions.R`)

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/exon_coverage_umicorrect")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

base_path <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/"

# Main figure
exp_list <- c("P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1")
names <- c("SYNnND VIS-1")

plot_list <- list()
for (i in seq_along(exp_list)) {
  p <- exon_coverage_scatter_plot(exp_list[i], names[i])
  plot_list[[i]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 1)
prefix <- glue("{plotDir}/exon_coverage_main")
save_png(combined_plot, prefix, 11, 5 * length(exp_list))
save_pdf(combined_plot, prefix, 11, 5 * length(exp_list))
```

---

## Figure 4D–H: Differential exon inclusion analysis

### Overview
Differential exon inclusion between neuronal (ExciteL23, ExciteL4, ExciteL5, ExciteL6, ExciteNP, InhibNeuron) and non-neuronal (Astrocytes, Microglia, OPCs, DivOPCs, COPs, MFOLs, MOLs) cell types was assessed using ScISOr-ATAC (`casesVcontrols` function).

A 2 x 2 contingency table was required to pass a chi-squared criterion (expected count >= 5 in majority of cells) before Fisher's exact test. P values were corrected using the Benjamini-Yekutieli method.

We compared enriched and control pools at three levels: all exons, exons passing chi-squared criterion, and significant exons (FDR < 0.05). A one-sided paired t-test was applied to paired observations from three datasets (SYNnND VIS-1, VIS-2, LA-SYNnND VIS-2) at each level.

### Step 1: Run scisorATAC

**Parameter:** `min_reads` = 5 (minimum reads for sum of 2 AllInfos per exon; default = 10)

```bash
exp="P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
sample="P56M8VIS_SYNnND_umicorrect"

cd "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}"

conda activate scisorATAC
mkdir scisorATAC_umicorrect && cd $_

# Define inputs
accept_allinfo="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"
control_allinfo="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"

# ---- Separate by cell type ----
mkdir accept_neuron_VS_nonNeuron_input && cd $_

# Neurons
zcat "${accept_allinfo}" | \
  awk '/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP|InhibNeuron/ \
  {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | \
  gzip -c > ${sample}_accept_AllInfo_neuron.gz

# Non-neurons
zcat "${accept_allinfo}" | \
  awk '/Astrocytes|Microglia|OPCs|DivOPCs|COPs|MFOLs|MOLs/ \
  {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | \
  gzip -c > ${sample}_accept_AllInfo_nonNeuron.gz

cd ..

mkdir control_neuron_VS_nonNeuron_input && cd $_

zcat "${control_allinfo}" | \
  awk '/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP|InhibNeuron/ \
  {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | \
  gzip -c > ${sample}_control_AllInfo_neuron.gz

zcat "${control_allinfo}" | \
  awk '/Astrocytes|Microglia|OPCs|DivOPCs|COPs|MFOLs|MOLs/ \
  {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | \
  gzip -c > ${sample}_control_AllInfo_nonNeuron.gz

cd ..

# ---- Run scisorATAC ----
PathToChromosomeFile="<REF_PATH>/mouse.all.chr.1-20.XY"
PathToAnnotation="<ANNO_PATH>/gencode.vM34.annotation.gtf.gz"
PathToCellTypeFile="<BASE_PATH>/SC_P56_M8_VIS/sc_mouse_celltype_file.txt"

NumThreads=12
ci_low=0.05
ci_high=0.95
MinNumReads=5
OL_fraction=0.8

# Accept: neuron vs non-neuron
out_path="accept_neuron_VS_nonNeuron"
mkdir "${out_path}" && cd $_

caseListPath="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/caseListPath_accept_neuron"
controlListPath="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/controlListPath_accept_nonNeuron"

touch "${caseListPath}" "${controlListPath}"
echo "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/${sample}_accept_AllInfo_neuron.gz" > "${caseListPath}"
echo "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/${sample}_accept_AllInfo_nonNeuron.gz" > "${controlListPath}"

echo "/tmp/" > tmpdirs
time bash <SCRIPTS_PATH>/v1.1a_exonInclusion_CTspecific_case_control.sh \
  "${caseListPath}" "${controlListPath}" tmpdirs \
  "${PathToChromosomeFile}" "${NumThreads}" "${PathToAnnotation}" \
  <SCRIPTS_PATH>/other-scripts/ 0.05 0.95 "${MinNumReads}" zcat 0.8 "${PathToCellTypeFile}" &> report

cd ..

# Control: neuron vs non-neuron
out_path="control_neuron_VS_nonNeuron"
mkdir "${out_path}" && cd $_

caseListPath="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/caseListPath_control_neuron"
controlListPath="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/controlListPath_control_nonNeuron"

touch "${caseListPath}" "${controlListPath}"
echo "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/${sample}_control_AllInfo_neuron.gz" > "${caseListPath}"
echo "<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/${sample}_control_AllInfo_nonNeuron.gz" > "${controlListPath}"

echo "/tmp/" > tmpdirs
time bash <SCRIPTS_PATH>/v1.1a_exonInclusion_CTspecific_case_control.sh \
  "${caseListPath}" "${controlListPath}" tmpdirs \
  "${PathToChromosomeFile}" "${NumThreads}" "${PathToAnnotation}" \
  <SCRIPTS_PATH>/other-scripts/ 0.05 0.95 "${MinNumReads}" zcat 0.8 "${PathToCellTypeFile}" &> report
```

### Step 2: Count exons at three levels

```bash
probe=Mouse.junction.probe.geneID.txt

# Level 1: All exons passing min_reads cutoff
# Columns: exon, inc-case, exc-case, inc-control, exc-control, sum-case, sum-control, dPSI
zcat cases_vs_controls.counts.tab.gz | \
  awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$2+$3"\t"$4+$5}' | \
  awk -F"\t" '$6 > 0 && $7 > 0 {print $0"\t"($2/$6)-($4/$7)}' | \
  sort -u > cases_vs_controls.counts.dPSI

# Add gene ID and filter for target genes
cat cases_vs_controls.counts.dPSI | \
  awk -F"\t" 'OFS="\t" {split($1,a,"_"); split(a[4],geneid,"."); print geneid[1],$1,$2,$3,$4,$5,$6,$7,$8}' | \
  sort | join -t $'\t' -1 1 -2 1 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9 "${probe}" - | \
  tr ' ' '\t' > cases_vs_controls.counts.dPSI.byjunctionGeneID

# Level 2: Exons tested by chi-squared
# Columns: exon, inc-case, exc-case, inc-control, exc-control, pval, FDR, LOR, sum-case, sum-control, dPSI
zcat cases_vs_controls.counts.pVal.FDR.LOR.tab.gz | \
  awk -F"\t" '{print $0"\t"$2+$3"\t"$4+$5}' | \
  awk -F"\t" '$9 > 0 && $10 > 0 {print $0"\t"($2/$9)-($4/$10)}' | \
  sort -u > cases_vs_controls.counts.pVal.FDR.LOR.dPSI

cat cases_vs_controls.counts.pVal.FDR.LOR.dPSI | \
  awk 'OFS="\t" {split($1,a,"_"); split(a[4],geneid,"."); print geneid[1],$0}' | \
  sort | join -t $'\t' -1 1 -2 1 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12 "${probe}" - | \
  tr ' ' '\t' > cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID

# Level 3: Significant exons (FDR < 0.05)
# FDR is column 8
cat cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID | \
  awk '{if ($8 <= 0.05) {print $0}}' > cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR

# Summary counts
cd ..
wc -l */cases_vs_controls.counts.dPSI.byjunctionGeneID                    # Total exons
wc -l */cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID       # Exons passing chi2
wc -l */cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR # Significant exons
```

### Step 3: Get statistics (Figure 4E–G)

```r
scisorATAC_stats <- function(enrich_file, control_file, title) {
  enrich_data <- read.table(enrich_file, header = FALSE)
  control_data <- read.table(control_file, header = FALSE)

  enrich_ids <- enrich_data[, 2]
  control_ids <- control_data[, 2]

  enrich_uniq_ids <- setdiff(enrich_ids, control_ids)
  control_uniq_ids <- setdiff(control_ids, enrich_ids)
  overlap_ids <- intersect(enrich_ids, control_ids)

  message(title)
  message("Enrich Control Overlap Enrich_unique Control_unique")
  message(paste(length(enrich_ids), length(control_ids),
                length(overlap_ids), length(enrich_uniq_ids),
                length(control_uniq_ids), sep = " "))
}

# Example usage for one experiment
exp <- "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"

# Total exons tested
enrich <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.dPSI.byjunctionGeneID")
control <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.dPSI.byjunctionGeneID")
scisorATAC_stats(enrich, control, "Num exon tested in scisorATAC:")

# Exons passing chi-sq criterion
enrich <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID")
control <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID")
scisorATAC_stats(enrich, control, "Num exon passed chisq criterion (tested):")

# Significant exons
enrich <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
control <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
scisorATAC_stats(enrich, control, "Num exon tested and are sig:")
```

### Step 4: Boxplot + paired t-test (Figure 4E–G)

```r
library(ggplot2)
library(ggsignif)
library(reshape2)
library(patchwork)
library(glue)

source("scripts/general/R_plot_functions.R")

# Define data
data_list <- list(
  total = list(
    df = data.frame(
      dataset = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"),
      Control = c(4878, 3294, 1576),
      Enriched = c(5096, 3936, 1977)
    ),
    header = "Total exon count"
  ),
  chi2_testable = list(
    df = data.frame(
      dataset = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"),
      Control = c(731, 380, 64),
      Enriched = c(961, 653, 132)
    ),
    header = "Chi2-testable exon count"
  ),
  sigFDR = list(
    df = data.frame(
      dataset = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"),
      Control = c(216, 154, 19),
      Enriched = c(307, 266, 46)
    ),
    header = "Number of sigFDR exons\n(neuron vs non-neuron)"
  )
)

plot_single <- function(df, header) {
  test_result <- t.test(df$Control, df$Enriched, paired = TRUE, alternative = "less")
  p_value <- test_result$p.value
  cat("Paired t-test p-value:", p_value, "\n")

  df$dataset <- factor(df$dataset, levels = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"))
  df_melted <- reshape2::melt(df, id.vars = c("dataset"))

  ggplot(df_melted, aes(x = variable, y = value, fill = variable)) +
    geom_point(aes(color = dataset), size = 2) +
    geom_boxplot(
      alpha = 0.3, outlier.shape = NA, width = 0.5, linewidth = 0.3,
      position = position_dodge(width = 0.3)
    ) +
    geom_signif(
      comparisons = list(c("Control", "Enriched")),
      annotations = glue("p = {round(p_value, 3)}"),
      y_position = max(df_melted$value) * 1.1,
      textsize = 3.5,
      tip_length = 0.03,
      size = 0.4,
      vjust = -0.5
    ) +
    labs(x = NULL, y = header, fill = NULL, color = "Dataset") +
    scale_fill_manual(values = c("Control" = "skyblue", "Enriched" = "#f2a93b")) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    coord_cartesian(ylim = c(0, max(df_melted$value) * 1.2)) +
    theme_jennie(plot_title_size = 10, base_size = 10) +
    guides(fill = "none")
}

# Combine plots
plot_total <- plot_single(data_list$total$df, data_list$total$header)
plot_chi2 <- plot_single(data_list$chi2_testable$df, data_list$chi2_testable$header)
plot_sigFDR <- plot_single(data_list$sigFDR$df, data_list$sigFDR$header)

combined_plot <- plot_total + plot_chi2 + plot_sigFDR

prefix <- glue("{plotDir}/paired_t_test_combined")
save_png(combined_plot, prefix, 14, 4)
save_pdf(combined_plot, prefix, 14, 4)
```

### Step 5: Bar plot of significant exons (Figure 4D)

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/scisorATAC_sig_exons_umicorrect/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

df_sigFDR <- data.frame(
  dataset = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"),
  Enriched = c(307, 266, 46),
  Control = c(216, 154, 19)
)
df_sigFDR$dataset <- factor(df_sigFDR$dataset,
                            levels = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"))

df_melted <- reshape2::melt(df_sigFDR, id.vars = c("dataset"))

colors <- c("Enriched" = "#f2a93b", "Control" = "skyblue")
bar_width <- 0.6

p <- ggplot(df_melted, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = bar_width),
    width = bar_width
  ) +
  labs(
    title = "Number of sigFDR exons\n(neuron vs non-neuron)",
    x = NULL, y = "Exon number", fill = "Pool"
  ) +
  geom_text(
    aes(label = sprintf("%.0f", value)),
    position = position_dodge(width = bar_width),
    vjust = -0.5, size = 4, family = "ArialMT"
  ) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, max(df_melted$value) * 1.1)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_jennie(plot_title_size = 12, base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

prefix <- glue("{plotDir}/neuron_vs_nonNeuron_sigFDR_exons")
save_png(p, prefix, 6, 5)
save_pdf(p, prefix, 6, 5)
```

### Step 6: Venn diagram (Figure 4H)

```r
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

source("scripts/general/R_plot_functions.R")
plotDir <- glue("{<PLOT_ROOT>}/scisorATAC_sig_exons_umicorrect/")

# Define experiment
exp <- "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
prefix <- glue("{plotDir}/venn_SYNnND_VIS2_umi.png")
title <- "sigFDR exons in enriched & control\n(SYNnND VIS-2)"

# Read significant exons
accept <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
control <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")

accept_ids <- as.vector(read.table(accept, header = FALSE)[, 2])
control_ids <- as.vector(read.table(control, header = FALSE)[, 2])

exon_list <- list(
  "Enriched" = accept_ids,
  "Control" = control_ids
)

myCol <- c("#F2A93B", "#87CEEB")

venn.diagram(
  x = exon_list,
  filename = prefix,
  output = TRUE,
  imagetype = "png",
  height = 1000,
  width = 1000,
  resolution = 300,
  lty = 'blank',
  fill = myCol,
  cex = 0.8,
  fontfamily = "Arial",
  fontface = "plain",
  cat.cex = 0.8,
  cat.fontfamily = "Arial",
  cat.fontface = "plain",
  cat.dist = c(0.07, 0.08),
  cat.pos = c(-50, 50),
  margin = 0.2,
  main = title,
  main.cex = 0.8,
  main.fontface = "plain",
  main.fontfamily = "Arial"
)
```

---

## Figure 4I–J: ScisorWiz exon usage plots

### Step 1: Rank unique significant exons

For all accepted unique exons, rank by dPSI and pick top genes for plotting.

```r
exp <- "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"

accept <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
control <- glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")

accept_data <- read.table(accept, header = FALSE)
control_data <- read.table(control, header = FALSE)

accept_ids <- accept_data[, 2]
control_ids <- control_data[, 2]

# Unique to enriched pool
accept_uniq_ids <- setdiff(accept_ids, control_ids)

# Filter and rank by absolute dPSI (column 12)
accept_uniq_ids_data <- accept_data[accept_data[, 2] %in% accept_uniq_ids, ]
accept_uniq_ids_data <- accept_uniq_ids_data[order(-abs(accept_uniq_ids_data[, 12])), ]

write.table(
  accept_uniq_ids_data,
  glue("<BASE_PATH>/Exome_Enrich_Major_Datasets/{exp}/scisorWiz_umicorrect/ranked_unique_accept_exons_{exp}.txt"),
  sep = "\t", row.names = FALSE, col.names = TRUE
)
```

### Step 2: Prepare 11-column AllInfo files

```bash
conda activate scisorATAC

# Example: M9 Prom SYNnND
exp="P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
sample="P56M9VIS_SYNnND_umicorrect"

accept_allinfo="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"

# Rename cell types
zcat "${accept_allinfo}" | awk '{
  gsub(/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP/, "ExcitatoryNeurons");
  gsub(/InhibNeuron/, "InhibitoryNeurons");
  gsub(/OPCs|DivOPCs/, "OPCs");
  gsub(/COPs|MFOLs|MOLs/, "Oligodendrocytes");
  print
}' | gzip -c > ${sample}_accept_AllInfo_rename.gz

# Convert to 11 columns
zcat ${sample}_accept_AllInfo_rename.gz | \
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | \
  gzip -c > ${sample}_accept_AllInfo_rename_11cols.gz

# Repeat for control
control_allinfo="<BASE_PATH>/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz"

zcat "${control_allinfo}" | awk '{
  gsub(/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP/, "ExcitatoryNeurons");
  gsub(/InhibNeuron/, "InhibitoryNeurons");
  gsub(/OPCs|DivOPCs/, "OPCs");
  gsub(/COPs|MFOLs|MOLs/, "Oligodendrocytes");
  print
}' | gzip -c > ${sample}_control_AllInfo_rename.gz

zcat ${sample}_control_AllInfo_rename.gz | \
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | \
  gzip -c > ${sample}_control_AllInfo_rename_11cols.gz
```

### Step 3: Extract reads covering target exons

```bash
# Example: M9 LA-SYNnND, gene Mtfr1
cd "<BASE_PATH>/Exome_Enrich_Major_Datasets/P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList/scisorWiz_umicorrect"

# Accept
zcat P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols.gz | \
  awk -v gene=ENSMUSG00000027601.14 -v exon=chr3_19262655_19262691 '
    BEGIN { split(exon, givenExon, "_") }
    {
      if ($2 != gene) next;
      n = split($9, a, ";%;");
      split(a[2], b, "_");
      if (b[2] > givenExon[3]) next;
      split(a[n], c, "_");
      if (c[3] < givenExon[2]) next;
      print;
    }
  ' | gzip -c > P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols_Mtfr1.gz

# Control (same command with control file)
```

### Step 4: Plot with ScisorWiz

```r
conda activate R_new

library(ScisorWiz)
library(glue)

# Define experiment
scisorwiz <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1/scisorWiz_umicorrect"
gene.query.list <- c("Atp8a1", "Cstf2")

cTypeFile <- "<BASE_PATH>/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212/scisorWiz/ct2color.txt"
annotation.file <- "<ANNO_PATH>/gencode.vM34.annotation.gtf.gz"

for (i in seq_along(gene.query.list)) {
  gene.query <- gene.query.list[i]

  # Enriched pool
  allInfoFile <- glue("{scisorwiz}/P56M9VIS_SYNnND_umicorrect_accept_AllInfo_rename_11cols_{gene.query}.gz")
  outpath <- glue("{scisorwiz}/accept/")
  output.dir1 <- paste0(outpath, gene.query)

  ScisorWiz_AllInfo(
    gencodeAnno = annotation.file,
    AllInfoInput = allInfoFile,
    cellTypeFile = cTypeFile,
    gene = gene.query,
    cluster = 1,
    ci = 0.05,
    mismatchCutoff = 0.05,
    outputDir = output.dir1,
    interactive = FALSE
  )

  # Control pool
  allInfoFile <- glue("{scisorwiz}/P56M9VIS_SYNnND_umicorrect_control_AllInfo_rename_11cols_{gene.query}.gz")
  outpath <- glue("{scisorwiz}/control/")
  output.dir1 <- paste0(outpath, gene.query)

  ScisorWiz_AllInfo(
    gencodeAnno = annotation.file,
    AllInfoInput = allInfoFile,
    cellTypeFile = cTypeFile,
    gene = gene.query,
    cluster = 1,
    ci = 0.05,
    mismatchCutoff = 0.05,
    outputDir = output.dir1,
    interactive = FALSE
  )
}
```
