# Analysis for Supplementary Figures & Tables

> This document contains the analysis pipeline for supplementary figures (S1–S18) and tables (S2–S5).

## Table of Contents

- [Figure S1: Transcript length distributions after cDNA size-selection](#figure-s1-transcript-length-distributions-after-cdna-size-selection)
- [Figure S2: Read length of size-selected vs non-size-selected libraries](#figure-s2-read-length-of-size-selected-vs-non-size-selected-libraries)
- [Figure S3: Read count summary across all sequencing runs](#figure-s3-read-count-summary-across-all-sequencing-runs)
- [Figure S5 & Table S2: Read length distributions between targeted and control pools](#figure-s5--table-s2-read-length-distributions-between-targeted-and-control-pools)
- [Figure S7: Enrichment ratios across experimental conditions](#figure-s7-enrichment-ratios-across-experimental-conditions)
- [Figure S8: Gene-level expression concordance across replicates](#figure-s8-gene-level-expression-concordance-across-replicates)
- [Figure S9: Concordance of gene-level enrichment ratios across replicates](#figure-s9-concordance-of-gene-level-enrichment-ratios-across-replicates)
- [Figure S10: Correlation between enrichment ratio and molecular characteristics](#figure-s10-correlation-between-enrichment-ratio-and-molecular-characteristics)
- [Figure S13: Scatter plot of non-target exons](#figure-s13-scatter-plot-of-non-target-exons)
- [Figure S14: Distribution of PSI and dPSI values in enriched pool](#figure-s14-distribution-of-psi-and-dpsi-values-in-enriched-pool)
- [Figure S15: Distribution of absolute dPSI for significant exons](#figure-s15-distribution-of-absolute-dpsi-for-significant-exons)
- [Figure S16: dPSI correlation between biological replicates](#figure-s16-dpsi-correlation-between-biological-replicates)
- [Figure S18 & Table S5: Cell barcode and UMI statistics](#figure-s18--table-s5-cell-barcode-and-umi-statistics)
- [Table S4: Gene sets and categories](#table-s4-gene-sets-and-categories)

---
## Figure S1: Transcript length distributions after cDNA size-selection

### Step 1: Extract read lengths

```bash
conda activate readfish

# M8 VIS control (size-selected)
cd <BASE_PATH>/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/analysis_read_length

zcat <BASE_PATH>/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo.gz | cut -f1 > SYNnND_VIS1_cutoff1_AllInfo_readIDs

readid=<BASE_PATH>/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/analysis_read_length/SYNnND_VIS1_cutoff1_AllInfo_readIDs
fastq=<BASE_PATH>/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212/fastq_by_condition/control/P56M1VIS_SYNnND_control.fastq.gz

seqkit grep -f "$readid" "$fastq" | awk 'FNR % 4 == 1 {id=$1} {getline; print id, length($0)}' > readID_2_length

# BICCN (non-size-selected)
cd <BASE_PATH>/SC_P56_M8_VIS/analysis_read_length

readid=<BASE_PATH>/SC_P56_M8_VIS/analysis_read_length/SC_P56_M8_VIS_cutoff1_AllInfo_readIDs
fastq=<BASE_PATH>/SC_P56_M8_VIS/fastq_by_condition/pass_control/P56_SC_M8_VIS_ONT.fastq.gz

seqkit grep -f "$readid" "$fastq" | awk 'FNR % 4 == 1 {id=$1} {getline; print id, length($0)}' > readID_2_length
```

### Step 2: Density plot

```r
source("scripts/general/R_plot_functions.R")

plotDir <- glue("{<PLOT_ROOT>}/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

anno <- "<ANNO_PATH>/gencode.vM34.annotation.gtf.gz"

# Parse GTF for protein-coding transcript lengths
gtf_data <- read.table(anno, sep = "\t", comment.char = "#", header = FALSE,
                       col.names = c("seqname", "source", "feature", "start", "end",
                                     "score", "strand", "frame", "attribute"))

df <- gtf_data %>%
  filter(feature == "exon") %>%
  mutate(
    transcript_id = str_match(attribute, "transcript_id ([^;]+);")[, 2],
    transcript_type = str_match(attribute, "transcript_type ([^;]+);")[, 2]
  ) %>%
  filter(transcript_type == "protein_coding")

pc_spliced_lengths <- df %>%
  group_by(transcript_id) %>%
  summarise(length = sum(end - start + 1), .groups = "drop")

# Load sample data
name <- "SYNnND VIS-1 (P)"
file_path <- "<BASE_PATH>/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/analysis_read_length/readID_2_length_dedup"
data_SYNnND_VIS1 <- read.csv(file_path, sep = " ", col.names = c("read_id", "seq_len")) %>%
  mutate(name = name) %>%
  select(seq_len, name) %>%
  mutate(seq_len = as.numeric(seq_len))

name <- "BICCN M8"
file_path <- "<BASE_PATH>/SC_P56_M8_VIS/analysis_read_length/readID_2_length_dedup"
data_BICCN_M8 <- read.csv(file_path, sep = " ", col.names = c("read_id", "seq_len")) %>%
  mutate(name = name) %>%
  select(seq_len, name) %>%
  mutate(seq_len = as.numeric(seq_len))

# Combine
df1 <- data.frame(length = pc_spliced_lengths$length, group = "Protein-coding transcripts")
df2 <- data.frame(length = data_SYNnND_VIS1$seq_len, group = "Size-selected VIS-1")
df3 <- data.frame(length = data_BICCN_M8$seq_len, group = "Non-size-selected VIS-1")

df <- rbind(df1, df2, df3)
df$group <- factor(df$group, levels = c("Protein-coding transcripts", "Non-size-selected VIS-1", "Size-selected VIS-1"))

# Density plot
p <- ggplot(df, aes(x = length, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(x = "Transcript length (bp)", y = "Density",
       title = "Density plot of transcript lengths",
       color = "Group", fill = "Group") +
  theme_jennie() +
  xlim(0, 6000)

save_png(p, glue("{plotDir}/transcript_length_density"), 10, 6)
save_pdf(p, glue("{plotDir}/transcript_length_density"), 10, 6)

# Summary statistics
for (name in c("Protein-coding transcripts", "Non-size-selected VIS-1", "Size-selected VIS-1")) {
  print(name)
  tmp <- df[df$group == name, ]
  quantiles <- quantile(tmp$length, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  print(glue("Median: {quantiles['50%']}"))
  print(glue("IQR: {quantiles['25%']} - {quantiles['75%']}"))
  print(glue("Middle 95% range: {quantiles['2.5%']} - {quantiles['97.5%']}"))
  print(glue("Mean: {mean(tmp$length)}"))
}
```

**Output summary:**

| Group | Median | IQR | Mean |
|-------|--------|-----|------|
| Protein-coding transcripts | 1842 | 834–3312 | 2444 |
| Non-size-selected VIS-1 | 913 | 745–1158 | 987 |
| Size-selected VIS-1 | 1817 | 1572–2071 | 1841 |

---

## Figure S2: Read length of size-selected vs non-size-selected libraries

```r
# Non-size-selected (BICCN)
sc_p56_m8_len <- read.csv(
  glue("<BASE_PATH>/SC_P56_M8_VIS-analysis_read_length/control_readLength.txt"),
  sep = "\t", header = FALSE
)

df_sc_p56_m8 <- data.frame(
  Sample = "SC_P56_M8_VIS",
  Type = "Non-size-selected",
  read_length = sc_p56_m8_len$V2
)

df_sc_p56_m8_sub <- df_sc_p56_m8 %>%
  group_by(Sample, Type) %>%
  slice_sample(n = 1000000, replace = FALSE) %>%
  ungroup()

# Size-selected
exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212",
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625",
  "P56M2HPC_minion_probe_all_junction_20240125",
  "P56M2HPC_minion_low_expression_junction_20240202"
)

name_list <- c(
  "SYNnND VIS-1 (P)",
  "SYN VIS-1 (P)",
  "SYNnND VIS-2 (P)",
  "LA-SYNnND VIS-2 (P)",
  "SYN VIS-2 (P)",
  "SYNnND HPC (M)",
  "LA-SYNnND HPC (M)"
)

root <- "<BASE_PATH>/Exome_Enrich_Major/"

df_all <- purrr::map2_dfr(exp_list, name_list, function(exp, sample_name) {
  raw_accept_len <- read.csv(glue("{root}/{exp}/analysis_read_length/accept_readLength.txt"), sep = "\t", header = FALSE)
  raw_control_len <- read.csv(glue("{root}/{exp}/analysis_read_length/control_readLength.txt"), sep = "\t", header = FALSE)

  bind_rows(
    data.frame(Sample = sample_name, Type = "Accepted (Raw)", read_length = raw_accept_len$V2),
    data.frame(Sample = sample_name, Type = "Control (Raw)", read_length = raw_control_len$V2)
  )
})

df_all$Sample <- factor(df_all$Sample, levels = name_list)

# Plot
color_mapping <- c(
  "Non-size-selected" = "#CE93D8",
  "Accepted (Raw)" = "#FFA8B8",
  "Control (Raw)" = "#80CBC4"
)

combined_df <- rbind(
  df_all[, c("Sample", "Type", "read_length")],
  df_sc_p56_m8[, c("Sample", "Type", "read_length")]
)

p_combined <- ggplot(combined_df, aes(x = Sample, y = read_length, fill = Type)) +
  geom_violin(position = position_dodge(width = 0.75), alpha = 0.5, width = 0.65, trim = TRUE, color = NA) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.15, outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = color_mapping) +
  labs(x = "", y = "Read length (bp)") +
  coord_cartesian(ylim = c(0, 5000)) +
  theme_jennie(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )

save_pdf(p_combined, "R1Q4_sequencing_length_combined", 12, 7)
```

---

## Figure S3: Read count summary across all sequencing runs

### Step 1: Extract read counts

```bash
cd <BASE_PATH>/P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList

# Use fastq_by_condition.sh to split by condition (comment out read length parts)
scripts/general/fastq_by_condition.sh "$fastq" "accept"
scripts/general/fastq_by_condition.sh "$fastq" "control"
```

### Step 2: Bar plot

```r
plotDir <- glue("{<PLOT_ROOT>}/supplementary/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

source("scripts/general/R_plot_functions.R")

format_read_length_files <- function(exp_list, names) {
  data_list <- list()

  for (i in seq_along(exp_list)) {
    exp <- exp_list[i]
    name <- names[i]

    accept_file_path <- file.path(base_path, exp, "analysis_read_length", "accept_readLength.txt")
    control_file_path <- file.path(base_path, exp, "analysis_read_length", "control_readLength.txt")

    if (file.exists(accept_file_path)) {
      temp <- read.table(accept_file_path, sep = "\t", header = FALSE)
      accept_read_count <- nrow(temp)

      temp <- read.table(control_file_path, sep = "\t", header = FALSE)
      control_read_count <- nrow(temp)

      data <- data.frame(
        name = name,
        accept_read_count = accept_read_count,
        accept_read_count_m = accept_read_count / 1000000,
        control_read_count = control_read_count,
        control_read_count_m = control_read_count / 1000000
      )
      data_list[[exp]] <- data
    } else {
      warning(paste("File not found:", accept_file_path))
    }
  }
  combined_data <- bind_rows(data_list)
  return(combined_data)
}

base_path <- "<BASE_PATH>/Exome_Enrich_Major_Datasets/"

exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212",
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625",
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList",
  "P56M2HPC_minion_probe_all_junction_20240125",
  "P56M2HPC_minion_low_expression_junction_20240202"
)

names <- c(
  "SYNnND VIS-1 (P)",
  "SYN VIS-1 (P)",
  "SYNnND VIS-2 (P)",
  "SYN VIS-2 (P)",
  "LA-SYNnND VIS-2 (P)",
  "SYNnND HPC (M)",
  "LA-SYNnND HPC (M)"
)

combined_data <- format_read_length_files(exp_list, names)
print(combined_data)

long_data <- pivot_longer(combined_data,
                          cols = c(accept_read_count_m, control_read_count_m),
                          names_to = "Pool",
                          values_to = "read_count") %>%
  mutate(pool = str_replace(Pool, "_read_count_m", ""))

colors <- c("accept_read_count_m" = "#f2a93b", "control_read_count_m" = "skyblue")
bar_width <- 0.8

p <- ggplot(long_data, aes(x = name, y = read_count, fill = Pool)) +
  geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
  labs(title = "Total number of enriched and control reads", x = NULL, y = "Read count (million)") +
  geom_text(aes(label = sprintf("%.1f", read_count)),
            position = position_dodge(width = bar_width), vjust = -0.5, size = 4, family = "ArialMT") +
  coord_cartesian(ylim = c(0, max(long_data$read_count) * 1.1)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_jennie(plot_title_size = 14, base_size = 14) +
  scale_fill_manual(values = colors,
                    labels = c("accept_read_count_m" = "Enriched", "control_read_count_m" = "Control")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

title <- "barplot_read_count_sevenRuns"
prefix <- glue("{plotDir}/{title}")
save_png(p, prefix, 9, 6)
save_pdf(p, prefix, 9, 6)
```

---

## Figure S5 & Table S2: Read length distributions between targeted and control pools

```r
experiments <- list(
  list(name = "SYNnND VIS-1 (P)", exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"),
  list(name = "SYNnND VIS-2 (P)", exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1")
)

root <- "<BASE_PATH>/Exome_Enrich_Major/"

process_experiment <- function(exp_name, exp_id) {
  raw_accept_len <- read.csv(glue("{root}/{exp_id}/timestamp/accept.readID_timestamp_readLength.tsv"), sep = "\t", header = TRUE)
  raw_control_len <- read.csv(glue("{root}/{exp_id}/timestamp/control.readID_timestamp_readLength.tsv"), sep = "\t", header = TRUE)

  filtered_accept_len <- read.csv(glue("{root}/{exp_id}/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_timebin"), sep = "\t",
    col.names = c("hour", "time_bin", "read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons"))

  filtered_control_len <- read.csv(glue("{root}/{exp_id}/timestamp/allinfo_timebin/control_bcFilter_AllInfo_timebin"), sep = "\t",
    col.names = c("hour", "time_bin", "read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons"))

  filtered_accept_len <- filtered_accept_len %>%
    select(read_id, hour, time_bin, gene_id, cellType, barcode, UMI) %>%
    left_join(raw_accept_len, by = "read_id")

  filtered_control_len <- filtered_control_len %>%
    select(read_id, hour, time_bin, gene_id, cellType, barcode, UMI) %>%
    left_join(raw_control_len, by = "read_id")

  df_plot <- bind_rows(
    data.frame(experiment = exp_name, type = "Accepted (Raw)", read_length = raw_accept_len$seq_len),
    data.frame(experiment = exp_name, type = "Accepted (Filtered)", read_length = filtered_accept_len$seq_len),
    data.frame(experiment = exp_name, type = "Control (Raw)", read_length = raw_control_len$seq_len),
    data.frame(experiment = exp_name, type = "Control (Filtered)", read_length = filtered_control_len$seq_len)
  )

  return(df_plot)
}

# Process all data
all_data <- list()
for (i in seq_along(experiments)) {
  exp_info <- experiments[[i]]
  df <- process_experiment(exp_info$name, exp_info$exp)
  all_data[[i]] <- df
}

combined_df <- bind_rows(all_data)
write.csv(combined_df, "R2Q3_part1.csv")

# Or load pre-computed
combined_df <- read.csv("R2Q3_part1_VIS1VIS2_accept_control_length.csv")

combined_df$type <- factor(combined_df$type,
  levels = c("Accepted (Raw)", "Control (Raw)", "Accepted (Filtered)", "Control (Filtered)"))

color_mapping <- c(
  "Accepted (Raw)" = "chocolate1",
  "Accepted (Filtered)" = "darkgoldenrod1",
  "Control (Raw)" = "#2ECC71",
  "Control (Filtered)" = "#3498DB"
)

p <- ggplot(combined_df, aes(x = type, y = read_length, fill = type)) +
  geom_violin(alpha = 0.5, width = 0.7, trim = TRUE, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Read length distribution across experiments", x = "", y = "Read length (bp)", fill = "Read Type") +
  coord_cartesian(ylim = c(0, 5000)) +
  theme_jennie() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~experiment, ncol = 2)

save_png(p, "R2Q3_length_distributions_combined", 10, 6, 300)
save_pdf(p, "R2Q3_length_distributions_combined", 10, 6)

# Export statistics table
tmp <- combined_df %>%
  group_by(experiment, type) %>%
  summarise(
    n_reads = n(),
    median_len = median(read_length, na.rm = TRUE),
    iqr_len = IQR(read_length, na.rm = TRUE),
    mean_len = mean(read_length, na.rm = TRUE),
    sd_len = sd(read_length, na.rm = TRUE),
    .groups = 'drop'
  )

write_csv(tmp, "R2Q3_part1_summary_VIS1VIS2.csv")
```

---

## Figure S7: Enrichment ratios across experimental conditions

```r
df <- tribble(
  ~Sample, ~GeneSet, ~Platform, ~EnrichRatio,
  "HPC",  "SYNnND",     "MinION",       1.57,
  "HPC",  "LA-SYNnND",  "MinION",       1.89,
  "VIS1", "SYNnND",     "PromethION",   1.43,
  "VIS2", "SYNnND",     "PromethION",   1.82,
  "VIS2", "LA-SYNnND",  "PromethION",   1.39,
  "VIS1", "SYN",        "PromethION",   1.82,
  "VIS2", "SYN",        "PromethION",   1.55
) %>% mutate(Label = paste(Sample, GeneSet, Platform, sep = "_"))

synnnd_df <- tribble(
  ~Sample, ~GeneSet, ~Platform, ~EnrichRatio,
  "HPC",  "SYNnND",     "MinION",       1.57,
  "VIS1", "SYNnND",     "PromethION",   1.43,
  "VIS2", "SYNnND",     "PromethION",   1.82,
) %>% mutate(Label = paste(Sample, GeneSet, Platform, sep = "_"))

vis2_df <- tribble(
  ~Sample, ~GeneSet, ~Platform, ~EnrichRatio,
  "VIS2", "SYNnND",     "PromethION",   1.82,
  "VIS2", "LA-SYNnND",  "PromethION",   1.39,
  "VIS2", "SYN",        "PromethION",   1.55
) %>% mutate(Label = paste(Sample, GeneSet, Platform, sep = "_"))

all_df <- bind_rows(
  df %>% mutate(Group = "All runs"),
  synnnd_df %>% mutate(Group = "SYNnND"),
  vis2_df %>% mutate(Group = "VIS2")
)

summary_df <- all_df %>%
  group_by(Group) %>%
  summarise(
    mean_enrich = mean(EnrichRatio),
    sd_enrich = sd(EnrichRatio),
    se_enrich = sd_enrich / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(xlab = paste0(Group, "\n", "n = ", n))

bar_fill <- "steelblue"
line_color <- "darkred"

combined_plot <- ggplot(summary_df, aes(x = xlab, y = mean_enrich)) +
  geom_col(width = 0.5, fill = bar_fill, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_enrich - sd_enrich, ymax = mean_enrich + sd_enrich),
                width = 0.25, size = 0.3, color = line_color) +
  geom_point(data = all_df %>% left_join(summary_df %>% select(Group, xlab), by = "Group"),
             aes(x = xlab, y = EnrichRatio),
             position = position_jitter(width = 0.12), size = 0.6, color = "black", alpha = 0.6,
             inherit.aes = FALSE) +
  labs(x = "", y = "Enrichment Ratio") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_jennie(base_size = 12) +
  theme(
    axis.line.x = element_line(size = 0.2),
    axis.line.y = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2),
    axis.title.y = element_text(size = 11)
  )

save_png(combined_plot, "R1Q4_enrich_ratio_error_bar", 4, 4)
save_pdf(combined_plot, "R1Q4_enrich_ratio_error_bar", 4, 4)
```

---

## Figure S8: Gene-level expression concordance across replicates

```r
base_path <- "genewise_enrich_ratio_tables"

# Read and normalize TPM data
read_count_TPM_df <- function(exp, sample_name) {
  df <- read.table(glue("{base_path}/geneID_acceptReads_cltReads_enrichRatio_{exp}"), sep = "\t", header = FALSE)
  colnames(df) <- c("gene_id", "accept_label", "accept_reads", "control_label", "control_reads", "enrichRatio")

  total_accept_reads <- sum(df$accept_reads, na.rm = TRUE)

  df2 <- df %>%
    transmute(gene_id, TPM = accept_reads / total_accept_reads * 1e6)

  colnames(df2)[2] <- sample_name
  df2
}

df_list <- purrr::map2(exp_list, names, read_count_TPM_df)
names(df_list) <- names

# Prepare comparison datasets
synnnd_df <- df_list[c("SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")
syn_df <- df_list[c("SYN VIS-1 (P)", "SYN VIS-2 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")
vis1_df <- df_list[c("SYN VIS-1 (P)", "SYNnND VIS-1 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")
vis2_df <- df_list[c("SYN VIS-2 (P)", "SYNnND VIS-2 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")

# Correlation plot function
create_corr_plot <- function(data, x_col, y_col, title, x_lab, y_lab) {
  ggplot(data, aes(x = log2(!!sym(x_col) + 1), y = log2(!!sym(y_col) + 1))) +
    geom_point(size = 0.8, color = "cyan3", alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 0.5, color = "gray50") +
    ggpubr::stat_cor(method = "pearson", aes(label = paste(..r.label..)),
                     label.x.npc = 0.01, label.y.npc = 0.9, color = "black", size = 3.5) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_jennie(plot_title_size = 11, base_size = 10)
}

p1 <- create_corr_plot(synnnd_df, "SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)", "Between samples (gene set = SYNnND)", "VIS-1", "VIS-2")
p2 <- create_corr_plot(syn_df, "SYN VIS-1 (P)", "SYN VIS-2 (P)", "Between samples (gene set = SYN)", "VIS-1", "VIS-2")
p3 <- create_corr_plot(vis1_df, "SYN VIS-1 (P)", "SYNnND VIS-1 (P)", "Between gene sets (sample = VIS-1)", "SYN", "SYNnND")
p4 <- create_corr_plot(vis2_df, "SYN VIS-2 (P)", "SYNnND VIS-2 (P)", "Between gene sets (sample = VIS-2)", "SYN", "SYNnND")

final_fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = margin(5, 5, 5, 5, "mm"),
        plot.tag = element_text(size = 10, face = "bold"))

plotDir <- "supplementary_figures_tables"
save_png(final_fig, glue("{plotDir}/R1Q2_replicates_read_count_TPM"), 6.5, 7)
save_pdf(final_fig, glue("{plotDir}/R1Q2_replicates_read_count_TPM"), 6.5, 7)
```

---

## Figure S9: Concordance of gene-level enrichment ratios across replicates

```r
read_enrich_df <- function(exp, sample_name) {
  df <- read.table(glue("{base_path}/geneID_acceptReads_cltReads_enrichRatio_{exp}"), sep = "\t", header = FALSE)
  colnames(df) <- c("gene_id", "accept_label", "accept_reads", "control_label", "control_reads", "enrichRatio")
  df2 <- df[, c("gene_id", "enrichRatio")]
  colnames(df2) <- c("gene_id", sample_name)
  df2
}

df_list <- purrr::map2(exp_list, names, read_enrich_df)
names(df_list) <- names

synnnd_df <- df_list[c("SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")
syn_df <- df_list[c("SYN VIS-1 (P)", "SYN VIS-2 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")
vis1_df <- df_list[c("SYN VIS-1 (P)", "SYNnND VIS-1 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")
vis2_df <- df_list[c("SYN VIS-2 (P)", "SYNnND VIS-2 (P)")] %>% purrr::reduce(inner_join, by = "gene_id")

# Quadrant statistics
q1_stats <- function(df, x, y, cutoff = 1) {
  n_total <- nrow(df)
  n_q1 <- sum(df[[x]] > cutoff & df[[y]] > cutoff, na.rm = TRUE)
  tibble::tibble(n_q1 = n_q1, pct_q1 = n_q1 / n_total,
                 label = sprintf("Q1: n = %d (%.1f%%)", n_q1, pct_q1 * 100))
}

q4_stats <- function(df, x, y, cutoff = 1) {
  n_total <- nrow(df)
  n_q4 <- sum(df[[x]] < cutoff & df[[y]] < cutoff, na.rm = TRUE)
  tibble::tibble(n_q4 = n_q4, pct_q4 = n_q4 / n_total,
                 label = sprintf("Q4: n = %d (%.1f%%)", n_q4, pct_q4 * 100))
}

create_plot <- function(data, x_col, y_col, title, axis_lim, label_x, label_y, text_y_offset) {
  q1_df <- q1_stats(data, x = x_col, y = y_col, cutoff = 1)
  q4_df <- q4_stats(data, x = x_col, y = y_col, cutoff = 1)

  cat(paste0(title, " sum pct: ", round(q1_df$pct_q1 + q4_df$pct_q4, 2), "\n"))

  ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(size = 0.8, color = "cyan3", alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 0.5, color = "cyan4") +
    ggpubr::stat_cor(method = "pearson",
                     aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                     label.x = label_x, label.y = label_y, color = "black", size = 3) +
    geom_text(data = q1_df, aes(x = label_x, y = label_y - text_y_offset[1], label = label),
              inherit.aes = FALSE, size = 3, color = "black", hjust = 0, vjust = 1) +
    geom_text(data = q4_df, aes(x = label_x, y = label_y - text_y_offset[2], label = label),
              inherit.aes = FALSE, size = 3, color = "black", hjust = 0, vjust = 1) +
    labs(title = title) +
    theme_jennie(plot_title_size = 11, base_size = 10) +
    xlim(0, axis_lim) + ylim(0, axis_lim) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7)
}

p1 <- create_plot(synnnd_df, "SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)", "Between samples (gene set = SYNnND)", axis_lim = 15, label_x = 6, label_y = 14, text_y_offset = c(1, 2.5))
p2 <- create_plot(syn_df, "SYN VIS-1 (P)", "SYN VIS-2 (P)", "Between samples (gene set = SYN)", 10, 4, 9.5, c(0.8, 1.6))
p3 <- create_plot(vis1_df, "SYN VIS-1 (P)", "SYNnND VIS-1 (P)", "Between gene sets (sample = VIS-1)", 9, 4, 8.5, c(0.8, 1.6))
p4 <- create_plot(vis2_df, "SYN VIS-2 (P)", "SYNnND VIS-2 (P)", "Between gene sets (sample = VIS-2)", 9, 4, 8.5, c(0.8, 1.6))

final_fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = margin(9, 9, 9, 9, "pt"),
        plot.tag = element_text(size = 10, face = "bold"))

save_png(final_fig, "R1Q2_replicates", 6.5, 7)
save_pdf(final_fig, "R1Q2_replicates", 6.5, 7)
```

---

## Figure S10: Correlation between enrichment ratio and molecular characteristics

```r
plotDir <- "supplementary_figures_tables"

# Annotation CDS length
gtf <- import("gencode.vM34.annotation.gtf")
cds_df <- gtf[gtf$type == "CDS"] %>%
  as.data.frame() %>%
  mutate(gene_id = sub("\\..*$", "", gene_id),
         transcript_id = sub("\\..*$", "", transcript_id))

tx_cds_lengths <- cds_df %>%
  group_by(gene_id, transcript_id) %>%
  summarise(cds_length = sum(width), .groups = "drop")

gene_mean_cds <- tx_cds_lengths %>%
  group_by(gene_id) %>%
  summarise(mean_cds_length = mean(cds_length), .groups = "drop")

# Set experiment
name <- "SYNnND VIS-1 (P)"
exp <- "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"

# Enrichment ratios
merged_df <- read.table(glue('genewise_enrich_ratio_tables/geneID_acceptReads_cltReads_enrichRatio_{exp}'), sep = '\t', header = FALSE)
colnames(merged_df) <- c("gene_id", "accept_label", "accept_reads", "control_label", "control_reads", "enrichRatio")
merged_df$gene_id <- sub("\\..*$", "", merged_df$gene_id)

# Average read length of gene targets
readlen <- read.csv(glue("{root}/{exp}/timestamp/control.readID_timestamp_readLength.tsv"), sep = "\t", header = TRUE)

# UMI-corrected gene expression
allinfo <- glue("{root}/{exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz")
allinfo <- read.csv(gzfile(allinfo), header = FALSE, sep = "\t")
colnames(allinfo) <- c("read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons")

# Scatter plot function
create_scatter_plot <- function(data, x_var, x_label, x_lim = NULL) {
  p <- ggplot(data, aes(x = !!sym(x_var), y = enrichRatio)) +
    geom_point(size = 0.8, color = "cyan3", alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 0.5, color = "cyan4") +
    stat_cor(method = "pearson",
             aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             label.x.npc = 0.5, label.y.npc = 0.9, color = "black", size = 3) +
    labs(title = name, x = x_label, y = "Enrich ratio", color = "Sample") +
    theme_jennie(plot_title_size = 9, base_size = 10) +
    theme(legend.position = "none")

  if (!is.null(x_lim)) {
    p <- p + xlim(x_lim[1], x_lim[2])
  }
  return(p)
}

# A) Average read length
read_len_df <- allinfo %>%
  dplyr::select(read_id, gene_id) %>%
  left_join(readlen, by = "read_id") %>%
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
  group_by(gene_id) %>%
  summarise(avg_read_length = mean(seq_len, na.rm = TRUE), .groups = "drop") %>%
  left_join(merged_df, by = "gene_id")

p1 <- create_scatter_plot(read_len_df, "avg_read_length", "Mean read length (bp)", x_lim = c(0, 4000))

# B) Annotation CDS length
annot_df <- merged_df %>%
  left_join(gene_mean_cds, by = "gene_id") %>%
  filter(!is.na(mean_cds_length))

p2 <- create_scatter_plot(annot_df, "mean_cds_length", "Mean transcript CDS length (bp)")

# C) UMI-corrected gene expression
gene_umi_counts <- allinfo %>%
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
  group_by(gene_id) %>%
  summarise(umi_count = n_distinct(UMI), .groups = "drop")

umi_df <- merged_df %>%
  left_join(gene_umi_counts, by = "gene_id")

p3 <- create_scatter_plot(umi_df, "umi_count", "UMI count")

# Combine
combined_plot <- (p1 / p2 / p3) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = margin(5, 5, 5, 5, "pt"),
        plot.tag = element_text(size = 10, face = "bold"))

save_png(combined_plot, glue("{plotDir}/R2Q3_enrichRatio_{name}"), 3.5, 8)
save_pdf(combined_plot, glue("{plotDir}/R2Q3_enrichRatio_{name}"), 3.5, 8)
```

---

## Figure S13: Scatter plot of non-target exons

```r
exon_coverage_offtarget_scatter_plot <- function(exp, name, subtitle) {
  parts <- strsplit(exp, "_")[[1]]
  device <- parts[2]
  sample <- parts[1]
  target <- paste(parts[3], parts[4], parts[5], sep = "_")
  title <- glue("{name}")

  # Read per off-target exon
  control_df <- read.table(glue('{base_path}/{exp}/exon_coverage2/control_bcFilter_readPerOfftargetExon.tsv'), sep = '\t', header = TRUE)
  accept_df <- read.table(glue('{base_path}/{exp}/exon_coverage2/accept_bcFilter_readPerOfftargetExon.tsv'), sep = '\t', header = TRUE)

  merged_df <- merge(control_df, accept_df, by = 'anno_exonblock', all.x = TRUE, all.y = TRUE) %>%
    replace_na(list(control_bcFilter_unique_readID_count = 0, accept_bcFilter_unique_readID_count = 0))

  larger_in_accept <- sum(merged_df$accept_bcFilter_unique_readID_count > merged_df$control_bcFilter_unique_readID_count)
  larger_in_control <- sum(merged_df$accept_bcFilter_unique_readID_count < merged_df$control_bcFilter_unique_readID_count)
  equal <- sum(merged_df$accept_bcFilter_unique_readID_count == merged_df$control_bcFilter_unique_readID_count)
  total <- nrow(merged_df)

  cor_test <- cor.test(merged_df$control_bcFilter_unique_readID_count, merged_df$accept_bcFilter_unique_readID_count)
  cor_coeff <- cor_test$estimate
  p_value <- cor_test$p.value

  # Dot frequency for coloring
  frequency_df <- merged_df %>%
    group_by(control_bcFilter_unique_readID_count, accept_bcFilter_unique_readID_count) %>%
    summarise(freq = n(), .groups = 'drop')

  merged_df <- merged_df %>%
    left_join(frequency_df, by = c("control_bcFilter_unique_readID_count", "accept_bcFilter_unique_readID_count"))

  # Bin labels and colors
  if (device == "minion") {
    bin_breaks <- c(-Inf, 10, 20, 50, Inf)
    bin_labels <- c("<=10", "10-20", "20-50", ">50")
  } else {
    bin_breaks <- c(-Inf, 5, 10, 50, Inf)
    bin_labels <- c("<=5", "5-10", "10-50", ">50")
  }
  bin_colors <- c("#66FFFF", "#FFCCFF", "#FF66FF", "#993FFF")
  names(bin_colors) <- bin_labels

  merged_df <- merged_df %>%
    mutate(freq_bin = cut(freq, breaks = bin_breaks, labels = bin_labels))

  max_limit <- max(merged_df$control_bcFilter_unique_readID_count, merged_df$accept_bcFilter_unique_readID_count)

  # Full scatter
  p1 <- ggplot(merged_df, aes(x = control_bcFilter_unique_readID_count,
                              y = accept_bcFilter_unique_readID_count, color = freq_bin)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
    labs(title = title, x = 'Reads per non-target exon in control',
         y = 'Reads per non-target exon in enriched', color = 'Freq') +
    xlim(0, max_limit) + ylim(0, max_limit) +
    scale_color_manual(values = bin_colors) +
    theme_jennie(plot_title_size = 10, base_size = 9) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none")

  # Zoomed scatter
  if (device == "minion") {
    max_limit_control <- calculate_max_exclude_99outliers(merged_df$control_bcFilter_unique_readID_count)
    max_limit_accept <- calculate_max_exclude_99outliers(merged_df$accept_bcFilter_unique_readID_count)
    max_limit <- max(max_limit_control, max_limit_accept)
  } else {
    max_limit_control <- calculate_max_excluding_outliers(merged_df$control_bcFilter_unique_readID_count)
    max_limit_accept <- calculate_max_excluding_outliers(merged_df$accept_bcFilter_unique_readID_count)
    max_limit <- max(max_limit_control, max_limit_accept)
  }

  p2 <- ggplot(merged_df, aes(x = control_bcFilter_unique_readID_count,
                              y = accept_bcFilter_unique_readID_count, color = freq_bin)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +
    labs(title = glue("enriched > control: {larger_in_accept} ({sprintf('%.1f', larger_in_accept/total * 100)}%)\nenriched = control: {equal} ({sprintf('%.1f', equal/total * 100)}%)\nenriched < control: {larger_in_control} ({sprintf('%.1f', larger_in_control/total * 100)}%)"),
         x = 'Reads per non-target exon in control',
         y = 'Reads per non-target exon in enriched', color = 'Freq') +
    xlim(0, max_limit) + ylim(0, max_limit) +
    scale_color_manual(values = bin_colors) +
    theme_jennie(plot_title_size = 9, base_size = 9) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    coord_fixed(ratio = 1) +
    theme(plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray40"))

  combined_plot <- (p1 + p2)
  return(combined_plot)
}

# Plot
base_path <- "<BASE_PATH>/Exome_Enrich_Major"

exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
)
names <- c("SYNnND VIS-1", "SYNnND VIS-2")

plot_list <- list()
for (i in seq_along(exp_list)) {
  p <- exon_coverage_offtarget_scatter_plot(exp_list[i], names[i])
  plot_list[[i]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(5, 5, 5, 5, "pt"),
        plot.tag = element_text(size = 10, face = "bold"))

prefix <- glue("{plotDir}/R3Q1_exon_coverage_offtarget")
save_pdf(combined_plot, prefix, 6, 3.5 * length(exp_list))
save_png(combined_plot, prefix, 6, 3.5 * length(exp_list))
```

---

## Figure S14: Distribution of PSI and dPSI values in enriched pool

### PSI distribution

```r
all_df_gene <- all_df %>%
  mutate(gene = str_extract(exon, "ENSMUSG\\d+"))

# Enriched pool PSI
df_enrich <- all_df_gene %>%
  filter(Pool == "Enriched") %>%
  mutate(case_PSI = Case_include / (Case_include + Case_exclude),
         control_PSI = Control_include / (Control_include + Control_exclude))

# Randomly select one exon per (Sample, type, gene)
set.seed(123)
df_sampled <- df_enrich %>%
  group_by(Sample, type, gene) %>%
  sample_n(size = 1) %>%
  ungroup()

# Pivot long
df_long <- df_sampled %>%
  select(Sample, type, gene, case_PSI, control_PSI) %>%
  pivot_longer(cols = c(case_PSI, control_PSI), names_to = "group", values_to = "PSI") %>%
  mutate(group = ifelse(group == "case_PSI", "Case", "Control"))

# PSI quantile stats
psi_quantiles <- df_long %>%
  group_by(Sample, type, group) %>%
  summarise(
    q25 = quantile(PSI, 0.25, na.rm = TRUE),
    q50 = median(PSI, na.rm = TRUE),
    q75 = quantile(PSI, 0.75, na.rm = TRUE),
    Count_total = n(),
    Count_0.25_0.75 = sum(PSI >= 0.25 & PSI <= 0.75, na.rm = TRUE),
    Percent_0.25_0.75 = Count_0.25_0.75 / Count_total * 100,
    .groups = "drop"
  ) %>%
  mutate(group = recode(group, "Case" = "Neuron", "Control" = "Glia")) %>%
  mutate(across(c(q25, q50, q75), ~ round(.x, 2))) %>%
  mutate(Percent_0.25_0.75 = round(Percent_0.25_0.75, 1))

print(psi_quantiles)
write_csv(psi_quantiles, "Q4-3_PSI_quantiles.csv")

# PSI plot
p1 <- ggplot(df_long, aes(x = PSI, fill = group)) +
  geom_density(alpha = 0.5, linewidth = 0.2) +
  facet_grid(Sample ~ type,
             labeller = labeller(type = c("all" = "Total exons",
                                          "chi2-testable" = "Chi2-testable exons",
                                          "sigFDR" = "Significant exons"))) +
  scale_fill_manual(values = c("Case" = "#80CBC4", "Control" = "#FFC2D6"),
                    labels = c("Case" = "Neuron", "Control" = "Glia")) +
  scale_y_continuous(breaks = function(lim) seq(floor(lim[1]), ceiling(lim[2]), by = 0.5)) +
  labs(x = "PSI", y = "Density", fill = "") +
  theme_jennie(plot_title_size = 14, axis_text_size = 10) +
  theme(legend.position = "bottom", legend.justification = "left", strip.text = element_text(hjust = 0.5))
```

### dPSI distribution

```r
set.seed(123)
df_sampled <- df_enrich %>%
  group_by(Sample, Pool, type, gene) %>%
  sample_n(size = 1) %>%
  ungroup()

p2 <- ggplot(df_sampled, aes(x = dPSI, fill = type)) +
  geom_density(alpha = 0.5, linewidth = 0.2) +
  facet_grid(Sample ~ Pool) +
  scale_fill_manual(values = c("all" = "#3498db", "chi2-testable" = "salmon", "sigFDR" = "#FFCC80"),
                    labels = c("all" = "Total exons", "chi2-testable" = "Chi2-testable exons",
                               "sigFDR" = "Significant exons (FDR<0.05)")) +
  labs(x = "dPSI", y = "Density", fill = "") +
  theme_jennie(plot_title_size = 14, axis_text_size = 10) +
  theme(legend.position = "bottom", legend.justification = "right", strip.text = element_text(hjust = 0.5))

# Combine
combined_plot <- (p1 + p2) +
  plot_layout(widths = c(3.5, 1)) +
  plot_annotation(title = "Distribution in enriched pool (random one exon per gene)",
                  theme = theme(plot.title = element_text(hjust = 0.5)), tag_levels = "A") &
  theme(plot.margin = margin(5, 5, 5, 5, "pt"),
        plot.tag = element_text(size = 14, face = "bold"))

save_png(combined_plot, glue("{plotDir}/Q4-3_detecting_partial_included_exons"), 9, 5, 300)
save_pdf(combined_plot, glue("{plotDir}/Q4-3_detecting_partial_included_exons"), 9, 5)
```

---

## Figure S15: Distribution of absolute dPSI for significant exons

```r
# Group 1: newly-discovered exons only sig after targeting
# Group 2: exons already significant before targeting

prepare_group_density_data <- function(sample_name, df_enrich, seed = 123) {
  # Chi2-testable data for this sample
  sample_chi2 <- all_df_gene %>%
    filter(type == "chi2-testable") %>%
    filter(Sample == sample_name)

  # Pivot wide
  sample_wide <- sample_chi2 %>%
    dplyr::select(exon, Pool, FDR) %>%
    pivot_wider(names_from = Pool, values_from = FDR, names_prefix = "FDR_") %>%
    filter(!is.na(FDR_Control) & !is.na(FDR_Enriched))

  # Define groups
  group1_exons <- sample_wide %>%
    filter(FDR_Control >= 0.05, FDR_Enriched < 0.05) %>%
    pull(exon)

  group2_exons <- sample_wide %>%
    filter(FDR_Control < 0.05) %>%
    pull(exon)

  # Randomly select one exon per gene
  set.seed(123)
  df_sampled <- df_enrich %>%
    filter(type == "chi2-testable", Sample == sample_name) %>%
    group_by(Sample, Pool, type, gene) %>%
    sample_n(size = 1) %>%
    ungroup()

  # Label groups
  group1_df <- df_sampled %>% filter(exon %in% group1_exons) %>% mutate(group = "Significant after targeting")
  group2_df <- df_sampled %>% filter(exon %in% group2_exons) %>% mutate(group = "Significant before targeting")

  combined_df <- bind_rows(group1_df, group2_df)
  return(combined_df)
}

# Process each sample
combined_vis1 <- prepare_group_density_data("VIS-1", df_enrich)
combined_vis2 <- prepare_group_density_data("VIS-2", df_enrich)

all_combined <- bind_rows(combined_vis1, combined_vis2)

# Box plot
all_combined$group <- factor(all_combined$group,
  levels = c("Significant before targeting", "Significant after targeting"))

p3 <- ggplot(all_combined, aes(x = group, y = abs(dPSI), fill = group)) +
  geom_boxplot(alpha = 0.6, linewidth = 0.2, outlier.size = 0.8) +
  facet_wrap(Sample ~ Pool) +
  scale_fill_manual(values = c("Significant before targeting" = "skyblue",
                               "Significant after targeting" = "#f2a93b")) +
  labs(x = NULL, y = "Absolute dPSI", fill = "") +
  theme_jennie(plot_title_size = 14, axis_text_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(hjust = 0.5),
        axis.text.x = element_blank(), plot.margin = margin(9, 9, 9, 9, "pt")) +
  geom_signif(comparisons = list(c("Significant before targeting", "Significant after targeting")),
              test = "wilcox.test", map_signif_level = FALSE, step_increase = 0.1,
              textsize = 3, vjust = -0.2)

save_png(p3, glue("{plotDir}/Q4-3_dPSI_before_after_targeting-box"), 5, 5, 300)
save_pdf(p3, glue("{plotDir}/Q4-3_dPSI_before_after_targeting-box"), 5, 5)
```

---

## Figure S16: dPSI correlation between biological replicates

### Scatter plots

```r
# 1. sigFDR exons
df_wide <- df_sig %>%
  select(exon, Sample, Pool, dPSI) %>%
  pivot_wider(names_from = Sample, values_from = dPSI) %>%
  filter(!is.na(`VIS-1`) & !is.na(`VIS-2`))

cor_val <- cor(df_wide$`VIS-1`, df_wide$`VIS-2`, method = "pearson")
cor_p <- cor.test(df_wide$`VIS-1`, df_wide$`VIS-2`, method = "pearson")$p.value

p1 <- ggplot(df_wide, aes(x = `VIS-1`, y = `VIS-2`)) +
  geom_point(alpha = 0.6, size = 2, color = "#2c7bb6") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 0.5, color = "cyan4") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 2,
           label = paste0(" r = ", round(cor_val, 3), "\nP = ", format(cor_p, digits = 3, scientific = TRUE))) +
  labs(title = "dPSI correlation (Enriched pool, FDR<0.05)") +
  theme_jennie(plot_title_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_fixed()

save_png(p1, glue("{plotDir}/Q4-4_sigFDR_correlation"), 6, 6, 300)

# 2. chi2-testable exons
df_wide <- df_chi %>%
  select(exon, Sample, Pool, dPSI) %>%
  pivot_wider(names_from = Sample, values_from = dPSI) %>%
  filter(!is.na(`VIS-1`) & !is.na(`VIS-2`))

cor_val <- cor(df_wide$`VIS-1`, df_wide$`VIS-2`, method = "pearson")
cor_p <- cor.test(df_wide$`VIS-1`, df_wide$`VIS-2`, method = "pearson")$p.value

p2 <- ggplot(df_wide, aes(x = `VIS-1`, y = `VIS-2`)) +
  geom_point(alpha = 0.6, size = 2, color = "#2c7bb6") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 0.5, color = "cyan4") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 2,
           label = paste0(" r = ", round(cor_val, 3), "\nP = ", format(cor_p, digits = 3, scientific = TRUE))) +
  labs(title = "dPSI correlation (Enriched pool, chi2-testable)") +
  theme_jennie(plot_title_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_fixed()

save_png(p2, glue("{plotDir}/Q4-4_chi2_correlation"), 6, 6, 300)

combined_plot <- (p2 + p1) +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5)), tag_levels = "A") &
  theme(plot.margin = margin(5, 5, 5, 5, "pt"),
        plot.tag = element_text(size = 14, face = "bold"))

save_png(combined_plot, glue("{plotDir}/Q4-4_dPSI_correlation"), 10, 5, 300)
save_pdf(combined_plot, glue("{plotDir}/Q4-4_dPSI_correlation"), 10, 5)
```

### Venn diagram

```r
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
myCol <- c("#F2A93B", "#87CEEB")

# 1. sigFDR exons
sig_vis1 <- df_sig %>% filter(Sample == "VIS-1") %>% pull(exon) %>% unique()
sig_vis2 <- df_sig %>% filter(Sample == "VIS-2") %>% pull(exon) %>% unique()

venn.diagram(
  x = list(VIS1 = sig_vis1, VIS2 = sig_vis2),
  category.names = c("VIS-1", "VIS-2"),
  filename = "Q4-4_venn_sigFDR.png",
  output = TRUE, imagetype = "png", height = 1000, width = 1000, resolution = 300,
  lty = 'blank', fill = myCol,
  cex = 0.8, fontfamily = "Arial", fontface = "plain",
  cat.cex = 0.8, cat.fontfamily = "Arial", cat.fontface = "plain",
  cat.dist = c(0.07, 0.08), cat.pos = c(-50, 50),
  margin = 0.2,
  main = "Overlap of significant exons (Enriched pool)",
  main.cex = 0.8, main.fontface = "plain", main.fontfamily = "Arial"
)

# 2. chi2-testable exons
sig_vis1 <- df_chi %>% filter(Sample == "VIS-1") %>% pull(exon) %>% unique()
sig_vis2 <- df_chi %>% filter(Sample == "VIS-2") %>% pull(exon) %>% unique()

venn.diagram(
  x = list(VIS1 = sig_vis1, VIS2 = sig_vis2),
  category.names = c("VIS-1", "VIS-2"),
  filename = "Q4-4_venn_chi2-testable.png",
  output = TRUE, imagetype = "png", height = 1000, width = 1000, resolution = 300,
  lty = 'blank', fill = myCol,
  cex = 0.8, fontfamily = "Arial", fontface = "plain",
  cat.cex = 0.8, cat.fontfamily = "Arial", cat.fontface = "plain",
  cat.dist = c(0.07, 0.08), cat.pos = c(-50, 50),
  margin = 0.2,
  main = "Overlap of chi2-testable exons (Enriched pool)",
  main.cex = 0.8, main.fontface = "plain", main.fontfamily = "Arial"
)
```

---

## Figure S18 & Table S5: Cell barcode and UMI statistics

```r
# Use allinfo.gz as input (not UMI-corrected)
root <- "<BASE_PATH>/Exome_Enrich_Major"
exp <- "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"

allinfo_VIS1 <- glue("{root}/{exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo.gz")
allinfo_VIS1 <- read.csv(gzfile(allinfo_VIS1), sep = "\t",
  col.names = c("read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons"))

exp <- "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
allinfo_VIS2 <- glue("{root}/{exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo.gz")
allinfo_VIS2 <- read.csv(gzfile(allinfo_VIS2), sep = "\t",
  col.names = c("read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons"))

# Stats per barcode
stats_VIS1 <- allinfo_VIS1 %>%
  group_by(barcode) %>%
  summarise(genes_per_cell = n_distinct(gene_id), umi_per_cell = n_distinct(UMI), .groups = "drop")

stats_VIS2 <- allinfo_VIS2 %>%
  group_by(barcode) %>%
  summarise(genes_per_cell = n_distinct(gene_id), umi_per_cell = n_distinct(UMI), .groups = "drop")

combined_stats <- bind_rows(
  stats_VIS1 %>% mutate(group = "SYNnND VIS-1"),
  stats_VIS2 %>% mutate(group = "SYNnND VIS-2")
)

# Summary table
summary_table <- data.frame(
  Metric = c("Total Reads", "Unique Cell Types", "Unique Barcodes", "Unique Barcode+UMI",
             "Genes per Barcode (mean ± sd)", "UMIs per Barcode (mean ± sd)"),
  `SYNnND VIS-1` = c(
    nrow(allinfo_VIS1),
    dplyr::n_distinct(allinfo_VIS1$cellType),
    dplyr::n_distinct(allinfo_VIS1$barcode),
    dplyr::n_distinct(paste(allinfo_VIS1$barcode, allinfo_VIS1$UMI)),
    sprintf("%.2f ± %.2f", mean(stats_VIS1$genes_per_cell), sd(stats_VIS1$genes_per_cell)),
    sprintf("%.2f ± %.2f", mean(stats_VIS1$umi_per_cell), sd(stats_VIS1$umi_per_cell))
  ),
  `SYNnND VIS-2` = c(
    nrow(allinfo_VIS2),
    dplyr::n_distinct(allinfo_VIS2$cellType),
    dplyr::n_distinct(allinfo_VIS2$barcode),
    dplyr::n_distinct(paste(allinfo_VIS2$barcode, allinfo_VIS2$UMI)),
    sprintf("%.2f ± %.2f", mean(stats_VIS2$genes_per_cell), sd(stats_VIS2$genes_per_cell)),
    sprintf("%.2f ± %.2f", mean(stats_VIS2$umi_per_cell), sd(stats_VIS2$umi_per_cell))
  ),
  stringsAsFactors = FALSE
)

summary_table
write_csv(summary_table, glue("{plotDir}/R2Q4_barcode_umi_stats.csv"))

# Plot
stats_long <- combined_stats %>%
  pivot_longer(cols = c(genes_per_cell, umi_per_cell), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("genes_per_cell", "umi_per_cell"),
                         labels = c("Genes per Cell", "UMI per Cell")))

colors <- c("SYNnND VIS-1" = "steelblue", "SYNnND VIS-2" = "salmon")

p <- ggplot(stats_long, aes(x = metric, y = value, fill = group)) +
  geom_violin(alpha = 0.5, width = 0.7, trim = TRUE, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black") +
  facet_wrap(~ group, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Count", fill = "Group") +
  theme_jennie(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

save_png(p, glue("{plotDir}/R2Q4_barcode_umi_stats"), 6.5, 4)
save_pdf(p, glue("{plotDir}/R2Q4_barcode_umi_stats"), 6.5, 4)
```

---

## Table S4: Gene sets and categories

```bash
# Create labeled gene set files
awk '{print $0 "\t" "SYNnND"}' "<ANNO_PATH>/Mouse.junction.probe.geneID.txt" > GENCODE_GRCm39_vM34_mouse_junction_probe_geneID.txt

awk '{print $0 "\t" "syngo"}' "<ANNO_PATH>/syngo_mouse_geneID" > GENCODE_GRCm39_vM34_syngo_mouse_geneID.txt

# LA-SYNnND VIS2 (1480 genes)
awk '{print $0 "\t" "LA_SYNnND_VIS2"}' "<BASE_PATH>/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes" > LA_SYNnND_VIS2_real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes.txt

# LA-SYNnND HPC (1518 genes)
awk '{print $0 "\t" "LA_SYNnND_HPC"}' "<BASE_PATH>/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly" > LA_SYNnND_HPC_P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly.txt
```

### Combine as Excel with multiple sheets

```r
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

dir_path <- "<BASE_PATH>/gene_sets_category"
txt_files <- list.files(path = dir_path, pattern = "\\.txt$", full.names = TRUE)

wb <- createWorkbook()

for (file in txt_files) {
  data <- read.table(file, header = TRUE, sep = "\t")
  sheet_name <- tools::file_path_sans_ext(basename(file))
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}

saveWorkbook(wb, "<BASE_PATH>/combined_gene_sets.xlsx", overwrite = TRUE)
```
