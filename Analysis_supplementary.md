# Analysis for supplementary figures & tables

## Figure S1 Transcript length distributions after cDNA size-selection 

```bash
conda activate readfish 

##### M8 VIS control allinfo reads 
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/analysis_read_length

zcat /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/scisorseqr/**control_bcFilter**/mmalign/LongReadInfo/AllInfo.gz | cut -f1 > SYNnND_VIS1_cutoff1_AllInfo_readIDs

readid=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/analysis_read_length/SYNnND_VIS1_cutoff1_AllInfo_readIDs 
fastq=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212/fastq_by_condition/control/P56M1VIS_SYNnND_control.fastq.gz

seqkit grep -f $readid $fastq | awk 'FNR % 4 == 1 {id=$1} {getline; print id, length($0)}' > readID_2_length

##### BICCN 
cd /home/xil4009/scratch_tilgnerlab/SC_P56_M8_VIS/analysis_read_length

readid=/home/xil4009/scratch_tilgnerlab/SC_P56_M8_VIS/analysis_read_length/SC_P56_M8_VIS_cutoff1_AllInfo_readIDs
fastq=/home/xil4009/scratch_tilgnerlab/SC_P56_M8_VIS/fastq_by_condition/pass_control/P56_SC_M8_VIS_ONT.fastq.gz

seqkit grep -f $readid $fastq | awk 'FNR % 4 == 1 {id=$1} {getline; print id, length($0)}' > readID_2_length

```

### Density plot

```r

source("scripts/general/R_plot_functions.R")
plotDir = glue("{plotRoot}/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
###########
  
anno = "/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/gencode.vM34.annotation.gtf.gz" 
plotDir = "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/all_plots/"

# Get gtf
gtf_data <- read.table(anno, sep = "\t", comment.char = "#", header = FALSE, 
                       col.names = c("seqname", "source", "feature", "start", "end", 
                                     "score", "strand", "frame", "attribute"))
# filter for exons of protein-coding transcript, extract trascript ID
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

exon_data <- gtf_data %>%
  filter(feature == "exon") %>%
  mutate(transcript_id = str_match(attribute, "transcript_id ([^;]+);")[,2])

head(exon_data)

# sum of lengths of exons for each tx 
spliced_lengths <- exon_data %>%
  group_by(transcript_id) %>%
  summarise(length = sum(end - start + 1), .groups = "drop") 

head(spliced_lengths)
mean(spliced_lengths$length)
median(spliced_lengths$length)

mean(lengths)
median(lengths)

mean(pc_spliced_lengths$length)
median(pc_spliced_lengths$length)

### SYNnND VIS1
name = "SYNnND VIS-1 (P)"
file_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/analysis_read_length/readID_2_length_dedup"
data_SYNnND_VIS1 <- read.csv(file_path, sep=" ", col.names=c("read_id", "seq_len")) %>% 
    mutate(name = name) %>% 
    select(seq_len, name) %>%
    mutate(seq_len = as.numeric(seq_len)) # %>% slice_sample(n = sample_size) 

# 1813264
lengths_SYNnND_VIS1 = data_SYNnND_VIS1$seq_len

###
name = "BICCN M8"
file_path <- "/home/xil4009/scratch_tilgnerlab/SC_P56_M8_VIS/analysis_read_length/readID_2_length_dedup"
data_BICCN_M8 <- read.csv(file_path, sep=" ", col.names=c("read_id", "seq_len")) %>%  # , nrows=10
    mutate(name = name) %>% 
    select(seq_len, name) %>%
    mutate(seq_len = as.numeric(seq_len)) # %>% slice_sample(n = sample_size) 

# 11967088
lengths_BICCN_M8 = data_BICCN_M8$seq_len

### 
df1 <- data.frame(length = pc_spliced_lengths$length, group = "Protein-coding transcripts")
df2 <- data.frame(length = lengths_SYNnND_VIS1, group = "Size-selected VIS-1")
df3 <- data.frame(length = lengths_BICCN_M8, group = "Non-size-selected VIS-1")

df <- rbind(df1, df2, df3)
head(df)

df$group <- factor(df$group, levels = c("Protein-coding transcripts", "Non-size-selected VIS-1", "Size-selected VIS-1"))

# density plot
p <- ggplot(df, aes(x = length, color = group, fill = group)) +
  geom_density(alpha = 0.5) +  
  labs(x = "Transcript length (bp)", y = "Density",
       title = "Density plot of transcript lengths",
       color = "Group", fill = "Group") +
  theme_jennie() + xlim(0,6000)
save_png(p, glue("{plotDir}/transcript_length_密度图"), 10, 6)
save_pdf(p, glue("{plotDir}/transcript_length_密度图"), 10, 6)

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

[1] "Protein-coding transcripts"
Median: 1842
IQR: 834 - 3312
Mean: 2443.96936906087

[1] "Non-size-selected VIS-1"
Median: 913
IQR: 745 - 1158
Mean: 987.090388655954

[1] "Size-selected VIS-1"
Median: 1817
IQR: 1572 - 2071
Mean: 1841.45561319256



## Figure S2 Read length of runs after long-molecule size-selection compared to non-size selected libraries

```r
### naive reads --------
sc_p56_m8_len <- read.csv(glue("/Users/jennie2000/Desktop/analysis_realtime_targeting/data_tmp/SC_P56_M8_VIS-analysis_read_length/control_readLength.txt"), sep="\t", header=FALSE)

df_sc_p56_m8 <- bind_rows(
	data.frame(
		Sample = "SC_P56_M8_VIS",
		Type = "Non-size-selected",
		read_length = sc_p56_m8_len$V2
	)
)

df_sc_p56_m8_sub <- df_sc_p56_m8 %>%
	group_by(Sample, Type) %>%
	slice_sample(
		n = 1000000,
		replace = FALSE
	) %>%
	ungroup()

### size-selected reads ------
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

root = "Exome_Enrich_Major/"

#  read and combine all exps
df_all <- purrr::map2_dfr(exp_list, name_list, function(exp, sample_name) {
	print(exp)
	raw_accept_len <- read.csv(glue("{root}/{exp}/analysis_read_length/accept_readLength.txt"), sep="\t", header=FALSE)
	raw_control_len <- read.csv(glue("{root}/{exp}/analysis_read_length/control_readLength.txt"), sep="\t", header=FALSE)

	df <- bind_rows(
		data.frame(
			Sample = sample_name,
			Type = "Accepted (Raw)",
			read_length = raw_accept_len$V2
		),
		data.frame(
			Sample = sample_name,
			Type = "Control (Raw)",
			read_length = raw_control_len$V2
		)
	)
	return(df)
})

df_all$Sample <- factor(df_all$Sample, levels = name_list)

# subsampling for testing visualization 
# df_sub <- df_all %>%
# 	group_by(Sample, Type) %>%
# 	slice_sample(
# 		n = 1000,
# 		replace = FALSE
# 	) %>%
# 	ungroup()
# 
# df_sub$Sample <- factor(df_sub$Sample, levels = name_list)

### plot 8 run distributions-------
combined_df <- rbind(
	df_all[, c("Sample", "Type", "read_length")],
	df_sc_p56_m8[, c("Sample", "Type", "read_length")]
)

color_mapping <- c(
	"Non-size-selected" = "#CE93D8",     # 玫粉色
	"Accepted (Raw)" = "#FFA8B8",     # 玫粉色
	"Control (Raw)" = "#80CBC4"      # 海绿色
)

p_combined <- ggplot(combined_df, aes(x = Sample, y = read_length, fill = Type)) +
	geom_violin(
		position = position_dodge(width = 0.75),
		alpha = 0.5,
		width = 0.65,
		trim = TRUE,
		color = NA
	) +
	geom_boxplot(
		position = position_dodge(width = 0.75),
		width = 0.15,
		outlier.shape = NA,
		alpha = 0.8
	) +
	scale_fill_manual(values = color_mapping) +
	labs(
		x = "",
		y = "Read length (bp)"
	) +
	coord_cartesian(ylim = c(0, 5000)) +
	theme_jennie(base_size = 14) +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position = "top",
		legend.title = element_blank(), 
		plot.margin = margin(2, 2, 2, 2, "pt")
	)  
p_combined
# save_png(p_combined, "R1Q4_sequencing_length_combined", 12, 7)
```


## Figure S3 Read count summary across all sequencing runs

Use read length output

```bash
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList

# comment out some part of this script, only use the part for read length 
scripts/general/fastq_by_condition.sh "$fastq" "accept"
scripts/general/fastq_by_condition.sh "$fastq" "control"
```

### bar plot

```r
plotDir = glue("{plotRoot}/supplementary/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

source("scripts/general/R_plot_functions.R")

# Combine accept_readLength.txt for each run
format_read_length_files <- function(exp_list, names) {
  data_list <- list()

  for (i in seq_along(exp_list)) {
    exp <- exp_list[i]
    name <- names[i]

    accept_file_path <- file.path(base_path, exp, "**analysis_read_length**", "accept_readLength.txt")
    control_file_path <- file.path(base_path, exp, "**analysis_read_length**", "control_readLength.txt")

****    if (file.exists(accept_file_path)) {
      temp <- read.table(accept_file_path, sep = "\t", header = FALSE)
            accept_read_count <- nrow(temp)

            temp <- read.table(control_file_path, sep = "\t", header = FALSE)
            control_read_count <- nrow(temp)

      data <- data.frame( name = name,
                          accept_read_count = accept_read_count,
                          accept_read_count_m = accept_read_count/1000000,
                          control_read_count = control_read_count,
                          control_read_count_m = control_read_count/1000000
                            )
      data_list[[exp]] <- data
    } else { warning(paste("File not found:", file_path)) }
  }
  combined_data <- bind_rows(data_list)
  return(combined_data)
}

####
# directories <- list.dirs(path = "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets", recursive = FALSE)

base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

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
  scale_y_continuous(expand = c(0, 0)) + # Set y-axis to start from 0
  theme_jennie(plot_title_size = 14, base_size = 14) +
  scale_fill_manual(values = colors,
                    labels = c("accept_read_count_m" = "Enriched", "control_read_count_m" = "Control") ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

title = "barplot_read_count_sevenRuns"
prefix = glue("{plotDir}/{title}")
save_png(p, prefix, 9, 6)
save_pdf(p, prefix, 9, 6)
```

## Figure S5, Table S2
Read length distributions are comparable between targeted and control pools

```r
experiments <- list(
	list(
		name = "SYNnND VIS-1 (P)",
		exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
	),
	list(
		name = "SYNnND VIS-2 (P)",
		exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716"
	)
)

root <- "Exome_Enrich_Major/"

process_experiment <- function(exp_name, exp_id) {
	raw_accept_len <- read.csv(
		glue("{root}/{exp_id}/timestamp/accept.readID_timestamp_readLength.tsv"),
		sep = "\t",
		header = TRUE
	)
	
	raw_control_len <- read.csv(
		glue("{root}/{exp_id}/timestamp/control.readID_timestamp_readLength.tsv"),
		sep = "\t",
		header = TRUE
	)
	
	filtered_accept_len <- read.csv(
		glue("{root}/{exp_id}/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_timebin"),
		sep = "\t",
		col.names = c("hour", "time_bin", "read_id", "gene_id", "cellType", 
					  "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons")
	)
	
	filtered_control_len <- read.csv(
		glue("{root}/{exp_id}/timestamp/allinfo_timebin/control_bcFilter_AllInfo_timebin"),
		sep = "\t",
		col.names = c("hour", "time_bin", "read_id", "gene_id", "cellType", 
					  "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons")
	)

	filtered_accept_len <- filtered_accept_len %>%
		select(read_id, hour, time_bin, gene_id, cellType, barcode, UMI) %>%
		left_join(raw_accept_len, by = "read_id")
	
	filtered_control_len <- filtered_control_len %>%
		select(read_id, hour, time_bin, gene_id, cellType, barcode, UMI) %>%
		left_join(raw_control_len, by = "read_id")
	
	# df for plot 
	df_plot <- bind_rows(
		data.frame(
			experiment = exp_name,
			type = "Accepted (Raw)",
			read_length = raw_accept_len$seq_len
		),
		data.frame(
			experiment = exp_name,
			type = "Accepted (Filtered)",
			read_length = filtered_accept_len$seq_len
		),
		data.frame(
			experiment = exp_name,
			type = "Control (Raw)",
			read_length = raw_control_len$seq_len
		),
		data.frame(
			experiment = exp_name,
			type = "Control (Filtered)",
			read_length = filtered_control_len$seq_len
		)
	)
	
	return(df_plot)
}

# Loop to process all data 
all_data <- list()

for (i in seq_along(experiments)) {
	exp_info <- experiments[[i]]
	df <- process_experiment(exp_info$name, exp_info$exp)
	all_data[[i]] <- df
}

combined_df <- bind_rows(all_data)
write.csv(combined_df, "R2Q3_part1_.csv")

combined_df = read.csv("R2Q3_part1_VIS1&VIS2_accept_control_length.csv")

combined_df$type <- factor(combined_df$type, 
						   levels = c("Accepted (Raw)", "Control (Raw)", "Accepted (Filtered)", "Control (Filtered)"))

# color mapping 
color_mapping <- c(
	"Accepted (Raw)" = "chocolate1",
	"Accepted (Filtered)" = "darkgoldenrod1",
	"Control (Raw)" = "#2ECC71",
	"Control (Filtered)" = "#3498DB"
)

# box + violin plot 
p <- ggplot(combined_df, aes(x = type, y = read_length, fill = type)) +
	geom_violin(alpha = 0.5, width = 0.7, trim = TRUE, color = NA) +
	geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
	scale_fill_manual(values = color_mapping) +
	labs(
		title = "Read length distribution across experiments",
		x = "",
		y = "Read length (bp)",
		fill = "Read Type"
	) +
	coord_cartesian(ylim = c(0, 5000)) +
	theme_jennie() +
	theme(
		legend.position = "none",
		strip.text.x = element_text(hjust = 0.5),
		axis.text.x = element_text(angle = 45, hjust = 1),
	) +
	facet_wrap(~experiment, ncol = 2)
p

save_png(p, "R2Q3_length_distributions_combined", 10, 6, 300)
save_pdf(p, "R2Q3_length_distributions_combined", 10, 6) 

# export statistics table 
tmp = combined_df %>%
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


## Figure S7 Gene-level expression concordance across replicates 
```r
base_path <- "genewise_enrich_ratio_tables"

# Read and normalize TPM data for each experiment
read_count_TPM_df <- function(exp, sample_name) {
	
	df <- read.table(
		glue("{base_path}/geneID_acceptReads_cltReads_enrichRatio_{exp}"),
		sep = "\t", header = FALSE
	)
	
	colnames(df) <- c(
		"gene_id",
		"accept_label", "accept_reads",
		"control_label", "control_reads",
		"enrichRatio"
	)
	print(head(df))
	
	# TPM normalization
	total_accept_reads <- sum(df$accept_reads, na.rm = TRUE)
	print(total_accept_reads)
	
	df2 <- df %>%
		transmute(
			gene_id,
			TPM = accept_reads / total_accept_reads * 1e6
		)
	print(head(df2))
	
	colnames(df2)[2] <- sample_name
	df2
}

df_list <- purrr::map2(exp_list, names, read_count_TPM_df)
names(df_list) <- names

# Prepare data for each comparison
synnnd_df <- df_list[c("SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")  # Merge two SYNnND samples

syn_df <- df_list[c("SYN VIS-1 (P)", "SYN VIS-2 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")  # Merge two SYN samples

vis1_df <- df_list[c("SYN VIS-1 (P)", "SYNnND VIS-1 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")  # Compare gene sets in VIS-1

vis2_df <- df_list[c("SYN VIS-2 (P)", "SYNnND VIS-2 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")  # Compare gene sets in VIS-2

# Function to create correlation plot with consistent formatting
create_corr_plot <- function(data, x_col, y_col, title, x_lab, y_lab) {
	ggplot(data, aes(x = log2(!!sym(x_col) + 1), y = log2(!!sym(y_col) + 1))) +
		geom_point(size = 0.8, color = "cyan3", alpha = 0.3) +
		geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 0.5, color = "gray50") +
		ggpubr::stat_cor(  # Add correlation coefficient and p-value
			method = "pearson",
			aes(label = paste(..r.label..)),
			label.x.npc = 0.01,  # Position at left
			label.y.npc = 0.9,   # Position near top
			color = "black",
			size = 3.5
		) + 
		labs(title = title, 
			 x = x_lab,
			 y = y_lab) +
		theme_jennie(plot_title_size=11, base_size=10)
}

# Create four correlation plots
p1 <- create_corr_plot(  # SYNnND between samples
	synnnd_df, "SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)",
	"Between samples (gene set = SYNnND)",
	"VIS-1", "VIS-2"
)

p2 <- create_corr_plot(  # SYN between samples
	syn_df, "SYN VIS-1 (P)", "SYN VIS-2 (P)",
	"Between samples (gene set = SYN)",
	"VIS-1", "VIS-2"
)

p3 <- create_corr_plot(  # Gene set comparison in VIS-1
	vis1_df, "SYN VIS-1 (P)", "SYNnND VIS-1 (P)",
	"Between gene sets (sample = VIS-1)",
	"SYN", "SYNnND"
)

p4 <- create_corr_plot(  # Gene set comparison in VIS-2
	vis2_df, "SYN VIS-2 (P)", "SYNnND VIS-2 (P)",
	"Between gene sets (sample = VIS-2)",
	"SYN", "SYNnND"
)

# Arrange plots in 2×2 grid with panel labels
final_fig <- (p1 | p2) /
	(p3 | p4) + 
	plot_annotation( tag_levels = 'A') &     
	theme(
		plot.margin = margin(5, 5, 5, 5, "mm"),
		plot.tag = element_text(size = 10, face = "bold") 
	)

plotDir = "supplementary_figures_tables"
save_png(final_fig, glue("{plotDir}/R1Q2_replicates_read_count_TPM"), 6.5, 7)
save_pdf(final_fig, glue("{plotDir}/R1Q2_replicates_read_count_TPM"), 6.5, 7)
```



## Figure S8 Concordance of gene‑level enrichment ratios across replicates

```r
read_enrich_df <- function(exp, sample_name) {
	df <- read.table(
		glue("{base_path}/geneID_acceptReads_cltReads_enrichRatio_{exp}"),
		sep = "\t", header = FALSE
	)
	
	colnames(df) <- c(
		"gene_id",
		"accept_label", "accept_reads",
		"control_label", "control_reads",
		"enrichRatio"
	)
	
	df2 <- df[, c("gene_id", "enrichRatio")]
	colnames(df2) <- c("gene_id", sample_name)
	
	df2
}


df_list <- purrr::map2(exp_list, names, read_enrich_df)
names(df_list) <- names

# Prepare datasets
synnnd_df <- df_list[c("SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")

syn_df <- df_list[c("SYN VIS-1 (P)", "SYN VIS-2 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")

vis1_df <- df_list[c("SYN VIS-1 (P)", "SYNnND VIS-1 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")

vis2_df <- df_list[c("SYN VIS-2 (P)", "SYNnND VIS-2 (P)")] %>%
	purrr::reduce(inner_join, by = "gene_id")


# function: calculate quantile 1 
q1_stats <- function(df, x, y, cutoff = 1) {
	n_total <- nrow(df)
	n_q1 <- sum(df[[x]] > cutoff & df[[y]] > cutoff, na.rm = TRUE)
	
	tibble::tibble(
		n_q1  = n_q1,
		pct_q1 = n_q1 / n_total,
		label = sprintf("Q1: n = %d (%.1f%%)", n_q1, pct_q1 * 100)
	)
}

# function: calculate quantile 4
q4_stats <- function(df, x, y, cutoff = 1) {
	n_total <- nrow(df)
	n_q4 <- sum(df[[x]] < cutoff & df[[y]] < cutoff, na.rm = TRUE)
	
	tibble::tibble(
		n_q4  = n_q4,
		pct_q4 = n_q4 / n_total,
		label = sprintf("Q4: n = %d (%.1f%%)", n_q4, pct_q4 * 100)
	)
}

# Function to create plot with consistent structure
create_plot <- function(data, x_col, y_col, title, axis_lim, 
						label_x, label_y, text_y_offset) {
	
	# Calculate quadrant statistics
	q1_df <- q1_stats(data, x = x_col, y = y_col, cutoff = 1)
	q4_df <- q4_stats(data, x = x_col, y = y_col, cutoff = 1)
	
	# Print sum percentage
	cat(paste0(title, " sum pct: ", round(q1_df$pct_q1 + q4_df$pct_q4, 2), "\n"))
	
	# Create plot
	ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
		geom_point(size = 0.8, color = "cyan3", alpha = 0.5) +
		geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 0.5, color = "cyan4") +
		ggpubr::stat_cor(
			method = "pearson",
			aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
			label.x = label_x,
			label.y = label_y,
			color = "black",
			size = 3
		) + 
		geom_text(
			data = q1_df,
			aes(x = label_x, y = label_y - text_y_offset[1], label = label),
			inherit.aes = FALSE, size = 3, color = "black", hjust = 0, vjust = 1
		) + 
		geom_text(
			data = q4_df,
			aes(x = label_x, y = label_y - text_y_offset[2], label = label),
			inherit.aes = FALSE, size = 3, color = "black", hjust = 0, vjust = 1
		) + 
		labs(title = title) +
		theme_jennie(plot_title_size = 11, base_size = 10) + 
		xlim(0, axis_lim) + ylim(0, axis_lim) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
		geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7)
}


# Create plots with customized parameters
p1 <- create_plot(
	synnnd_df, "SYNnND VIS-1 (P)", "SYNnND VIS-2 (P)",
	"Between samples (gene set = SYNnND)",
	axis_lim=15, label_x=6, label_y=14, text_y_offset=c(1, 2.5)
)

p2 <- create_plot(
	syn_df, "SYN VIS-1 (P)", "SYN VIS-2 (P)",
	"Between samples (gene set = SYN)",
	10, 4, 9.5, c(0.8, 1.6)
)

p3 <- create_plot(
	vis1_df, "SYN VIS-1 (P)", "SYNnND VIS-1 (P)",
	"Between gene sets (sample = VIS-1)",
	9, 4, 8.5, c(0.8, 1.6)
)

p4 <- create_plot(
	vis2_df, "SYN VIS-2 (P)", "SYNnND VIS-2 (P)",
	"Between gene sets (sample = VIS-2)",
	9, 4, 8.5, c(0.8, 1.6)
)

# Arrange plots in 2×2 grid
final_fig <- (p1 | p2) /
	(p3 | p4) + 
	plot_annotation(tag_levels = 'A') &  
	theme(plot.margin = margin(9, 9, 9, 9, "pt"),
		  plot.tag = element_text(size = 10, face = "bold") 
	)  

save_png(final_fig, "R1Q2_replicates", 6.5, 7)
save_pdf(final_fig, "R1Q2_replicates", 6.5, 7)
```


## Figure S9 The correlation between enrichment ratio and several molecular characteristics 

```r
plotDir = "supplementary_figures_tables"

# anno cds length 
gtf <- import("gencode.vM34.annotation.gtf")
cds_df <- gtf[gtf$type == "CDS"] %>%
	as.data.frame() %>%
	mutate(
		gene_id = sub("\\..*$", "", gene_id),
		transcript_id = sub("\\..*$", "", transcript_id)
	)

tx_cds_lengths <- cds_df %>%
	group_by(gene_id, transcript_id) %>%
	summarise(cds_length = sum(width), .groups = "drop")

gene_mean_cds <- tx_cds_lengths %>%
	group_by(gene_id) %>%
	summarise(mean_cds_length = mean(cds_length), .groups = "drop")


# set experiment 
name = "SYNnND VIS-1 (P)"
exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"

name = "SYNnND VIS-2 (P)"
exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"

# enrich ratios
merged_df <- read.table(
	glue('genewise_enrich_ratio_tables/geneID_acceptReads_cltReads_enrichRatio_{exp}'),
	sep = '\t', header = FALSE
)
colnames(merged_df) <- c("gene_id", "accept_label", "accept_reads",
						 "control_label", "control_reads", "enrichRatio")
merged_df$gene_id <- sub("\\..*$", "", merged_df$gene_id)

# average length of gene targets 
root = "/Volumes/Untitled/实验备份/Exome_Enrich_Major"
readlen <- read.csv(glue("{root}/{exp}/timestamp/control.readID_timestamp_readLength.tsv"),
					sep="\t", header=TRUE)

# umi corrected gene expression
allinfo = glue("{root}/{exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz")
allinfo <- read.csv(gzfile(allinfo), header = FALSE, sep = "\t")
colnames(allinfo) = c("read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons") 

# Function to create scatter plot with correlation
create_scatter_plot <- function(data, x_var, x_label, x_lim = NULL) {
	
	p <-ggplot(data, aes(x = !!sym(x_var), y = enrichRatio)) +
		geom_point(size = 0.8, color = "cyan3", alpha = 0.5) +
		geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 0.5, alpha = 0.7) +
		geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 0.5, color = "cyan4") +
		stat_cor(
			method = "pearson",
			aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
			label.x.npc = 0.5,
			label.y.npc = 0.9,
			color = "black",
			size = 3
		) +
		labs(
			title = name,
			x = x_label,
			y = "Enrich ratio",
			color = "Sample"
		) +
		theme_jennie(plot_title_size = 9, base_size = 10) +
		theme( legend.position = "none" )
	
	# Apply x-axis limits if provided
	if (!is.null(x_lim)) {
		p <- p + xlim(x_lim[1], x_lim[2])
	}
	
	return(p)
}

# A) Average read length of gene targets

# read_len_df <- allinfo %>%
# 	dplyr::select(read_id, gene_id) %>%
# 	left_join(readlen, by = "read_id") %>%
# 	mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
# 	group_by(gene_id) %>%
# 	summarise(
# 		avg_read_length = mean(seq_len, na.rm = TRUE),
# 		.groups = "drop"
# 	) %>%
# 	left_join(merged_df, by = "gene_id")
# write.csv(read_len_df, glue("{plotDir}/R2Q3_read_len_df_{name}.csv"))

read_len_df = read.csv(glue("{plotDir}/R2Q3_read_len_df_{name}.csv"))
p1 <- create_scatter_plot(
	read_len_df, "avg_read_length", 
	"Mean read length (bp)",
	x_lim = c(0, 4000)
)

# B) Annotation CDS length

# annot_df <- merged_df %>%
# 	left_join(gene_mean_cds, by = "gene_id") %>%
# 	filter(!is.na(mean_cds_length))
# write.csv(annot_df, glue("{plotDir}/R2Q3_annot_df_{name}.csv"))

annot_df = read.csv(glue("{plotDir}/R2Q3_annot_df_{name}.csv"))
p2 <- create_scatter_plot(
	annot_df, "mean_cds_length",
	"Mean transcript CDS length (bp)"
)

# C) UMI corrected gene expression

# gene_umi_counts <- allinfo %>%
# 	mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
# 	group_by(gene_id) %>%
# 	summarise(
# 		umi_count = n_distinct(UMI),
# 		.groups = "drop"
# 	)
# umi_df <- merged_df %>%
# 	left_join(gene_umi_counts, by = "gene_id")
# write.csv(umi_df, glue("{plotDir}/R2Q3_umi_df_{name}.csv"))

umi_df = read.csv(glue("{plotDir}/R2Q3_umi_df_{name}.csv"))
p3 <- create_scatter_plot(
	umi_df, "umi_count",
	"UMI count"
)

## Plot A+B+C --------
if (name == "SYNnND VIS-1 (P)") {
	combined_plot <- (p1 / p2 / p3) + 
		plot_annotation( tag_levels = 'A') &  
		theme( plot.margin = margin(5, 5, 5, 5, "pt"),
			   plot.tag = element_text(size = 10, face = "bold") 
		)  
	save_png(combined_plot, glue("{plotDir}/R2Q3_enrichRatio_{name}"), 3.5, 8)
	save_pdf(combined_plot, glue("{plotDir}/R2Q3_enrichRatio_{name}"), 3.5, 8)
} else {
	combined_plot_2 <- (p1 / p2 / p3) + 
		plot_annotation( tag_levels = 'A') &  
		theme( plot.margin = margin(5, 5, 5, 5, "pt"),
			   plot.tag = element_text(size = 10, face = "bold", color = "white") 
		)
	save_png(combined_plot_2, glue("{plotDir}/R2Q3_enrichRatio_{name}"), 3.5, 8)
	save_pdf(combined_plot_2, glue("{plotDir}/R2Q3_enrichRatio_{name}"), 3.5, 8)
}
```

## Table S4, Figure S13 Cell barcode and UMI 

```r
# use allinfo.gz as input（not UMI-corrected）

root = "/Volumes/Untitled/实验备份/Exome_Enrich_Major"
exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"

allinfo_VIS1 = glue("{root}/{exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo.gz")

allinfo_VIS1 <- read.csv(gzfile(allinfo_VIS1), sep = "\t",
						 col.names = c("read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons")) 


exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"

allinfo_VIS2 = glue("{root}/{exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo.gz")
allinfo_VIS2 <- read.csv(gzfile(allinfo_VIS2), sep = "\t", 
						 col.names = c("read_id", "gene_id", "cellType", "barcode", "UMI", "intronChain", "exonChain", "status", "numIntrons")) 

# stats
stats_VIS1 <- allinfo_VIS1 %>%
	group_by(barcode) %>%
	summarise(
		genes_per_cell = n_distinct(gene_id),
		umi_per_cell = n_distinct(UMI),
		.groups = "drop"
	)

# stats  
stats_VIS2 <- allinfo_VIS2 %>%
	group_by(barcode) %>%
	summarise(
		genes_per_cell = n_distinct(gene_id),
		umi_per_cell = n_distinct(UMI),
		.groups = "drop"
	)

combined_stats <- bind_rows(
	stats_VIS1 %>% mutate(group = "SYNnND VIS-1"),
	stats_VIS2 %>% mutate(group = "SYNnND VIS-2")
)

summary_table = data.frame(
	Metric = c("Total Reads",
			   "Unique Cell Types",
			   "Unique Barcodes",
			   "Unique Barcode+UMI",
			   "Genes per Barcode (mean ± sd)",
			   "UMIs per Barcode (mean ± sd)"),
	
	`SYNnND VIS-1` = c(
		nrow(allinfo_VIS1),
		dplyr::n_distinct(allinfo_VIS1$cellType),
		dplyr::n_distinct(allinfo_VIS1$barcode),
		dplyr::n_distinct(paste(allinfo_VIS1$barcode, allinfo_VIS1$UMI)),
		sprintf("%.2f ± %.2f",
				mean(stats_VIS1$genes_per_cell),
				sd(stats_VIS1$genes_per_cell)),
		sprintf("%.2f ± %.2f",
				mean(stats_VIS1$umi_per_cell),
				sd(stats_VIS1$umi_per_cell))
	),
	
	`SYNnND VIS-2` = c(
		nrow(allinfo_VIS2),
		dplyr::n_distinct(allinfo_VIS2$cellType),
		dplyr::n_distinct(allinfo_VIS2$barcode),
		dplyr::n_distinct(paste(allinfo_VIS2$barcode, allinfo_VIS2$UMI)),
		sprintf("%.2f ± %.2f",
				mean(stats_VIS2$genes_per_cell),
				sd(stats_VIS2$genes_per_cell)),
		sprintf("%.2f ± %.2f",
				mean(stats_VIS2$umi_per_cell),
				sd(stats_VIS2$umi_per_cell))
	),
	
	stringsAsFactors = FALSE
)
summary_table

write_csv(summary_table, glue("{plotDir}/R2Q4_barcode_umi_stats.csv"))


# df for plotting
stats_long <- combined_stats %>%
	pivot_longer(
		cols = c(genes_per_cell, umi_per_cell),
		names_to = "metric",
		values_to = "value"
	) %>%
	mutate(
		metric = factor(metric, 
						levels = c("genes_per_cell", "umi_per_cell"),
						labels = c("Genes per Cell", "UMI per Cell"))
	)


colors <- c("SYNnND VIS-1" = "steelblue", "SYNnND VIS-2" = "salmon")

# plot boxplot + violin
p = ggplot(stats_long, aes(x = metric, y = value, fill = group)) +
	geom_violin(alpha = 0.5, width = 0.7, trim = TRUE, color = NA) +
	geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black") +
	facet_wrap(~ group, scales = "free_y", nrow = 1) +
	scale_fill_manual(values = colors) +
	labs(
		x = "",
		y = "Count",
		fill = "Group"
	) +
	theme_jennie(base_size = 10) +
	theme(
		legend.position = "none",
		strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
		axis.text.x = element_text(angle = 0, hjust = 0.5)
	)
save_png(p, glue("{plotDir}/R2Q4_barcode_umi_stats"), 6.5, 4)
save_pdf(p, glue("{plotDir}/R2Q4_barcode_umi_stats"), 6.5, 4)
```



--------

## Table S2: Gene sets and category

```bash
# probe=/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/syngo_mouse_geneID
# probe="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
# probe="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
# probe="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/als"
# probe="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/ad"
# probe="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/autism"

awk '{print $0 "\t" "SYNnND"}' "/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt" > GENCODE_GRCm39_vM34_mouse_junction_probe_geneID.txt

awk '{print $0 "\t" "syngo"}' "/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/syngo_mouse_geneID" > GENCODE_GRCm39_vM34_syngo_mouse_geneID.txt

LA-SYNnND VIS2 (1480): 
/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes

awk '{print $0 "\t" "LA_SYNnND_VIS2"}' "/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes" > LA_SYNnND_VIS2_real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes.txt

LA-SYNnND HPC (1518): 
/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly

awk '{print $0 "\t" "LA_SYNnND_HPC"}' "/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly" > LA_SYNnND_HPC_P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly.txt

```

Combine as one Excel multiple sheets

```r
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# Gene set txt files 
dir_path <- "/Users/jennie2000/Downloads/gene_sets_category"
txt_files <- list.files(path = dir_path, pattern = "\\.txt$", full.names = TRUE)

# Init new workbook 
wb <- createWorkbook()

for (file in txt_files) {
  data <- read.table(file, header = TRUE, sep = "\t")
  
  # Create sheet name
  sheet_name <- tools::file_path_sans_ext(basename(file))
  addWorksheet(wb, sheet_name)
  
  writeData(wb, sheet_name, data)
}

saveWorkbook(wb, "/Users/jennie2000/Downloads/combined.xlsx", overwrite = TRUE)
```
