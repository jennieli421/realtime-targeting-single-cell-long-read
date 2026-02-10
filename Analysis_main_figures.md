# Analysis for main figures

## Evaluation of size-selection protocol 
The sequencing data in this study used the preserved cDNAs from our previous study [39] but the libraries were subject to an additional size-selection procedure to obtain longer fragments. To evaluate the effect of this protocol, we compared the dataset generated from current and from previous study regarding number of exons per read. Specifically, we compared data of sample VIS-1 with targeting SYNnND to data of the same sample (P56_M8_VIS) with probe enrichment from previous study [39]. 
The exons covered by each read were identified by overlapping exon chains from Scisorseq “AllInfo.gz” with the sorted exon blocks from GENCODE vM34 annotation, using “bedtools intersect” [56]. 
100k reads from each dataset were randomly sampled and the exon count per read was determined. The subsampling process was repeated 10 times and the distributions of the two datasets for each iteration were presented in the box plot shown in Figure 2A. 

Use `scripts/scisorseq_scripts/exon_coverage.py`

```bash
anno_exons=resources/annotations/mouse_vM34_exon_count/anno.exon.sorted.bed # sorted annotation exon blocks
biccn_allinfo=""
allinfo=""

# convert allinfo exon coords into sorted bed format
# column: chr start stop readID_exonblock geneID strand

sample="SC_P56_M8_VIS"
zcat "$biccn_allinfo" | awk '{ n=split($7,a,";%;"); for(i=2;i<=n;i++){print a[i]"\t"$2"\t"$1} }' | sort -u | awk 'OFS="\t" {split($1,a,"_"); print a[1], a[2], a[3], $3"_"$1, $2, a[4]}' | sort -k1,1 -k2,2n | uniq > ${sample}_allinfo_exon_sorted.bed # forgot to remove version

sample="VIS1_SYNnND_prom_accept_bcFilter"
zcat "$allinfo" | awk '{ n=split($7,a,";%;"); for(i=2;i<=n;i++){print a[i]"\t"$2"\t"$1} }' | sort -u | awk 'OFS="\t" {split($1,a,"_"); print a[1], a[2], a[3], $3"_"$1, $2, a[4]}' | sort -k1,1 -k2,2n | uniq > ${sample}_allinfo_exon_sorted.bed # forgot to remove version

###### shared exons
bedtools intersect -a ${sample}_allinfo_exon_sorted.bed -b "$anno_exons" -sorted -wao | sort | uniq > ${sample}_allinfo_exon_intersect_annotation

###### after getting the allinfo_intersect_annotation output
###### output: 1. all exons; 2. exons of ontarget genes
targets=resources/gene_target_lists/Mouse.junction.probe.geneID.txt
script=scripts/scisorseq_scripts/exon_coverage.py
python $script --mode "allinfo_to_annotation_exonblock" --input ${sample}_allinfo_exon_intersect_annotation --out_prefix ${sample} --targets $targets
```

**subsample 10 times**

```python
cd /home/xil4009/scratch_tilgnerlab/size_selection_compare_M8_cutoff1
python

import pandas as pd
import numpy as np

input_file = "SC_P56_M8_VIS_Ontarget_allinfo_exonblock_to_annotation_exonblock.tsv"
output_prefix="SC_P56_M8_VIS"

input_file = "VIS1_SYNnND_prom_accept_bcFilter_Ontarget_allinfo_exonblock_to_annotation_exonblock.tsv"
output_prefix="VIS1_SYNnND_prom_accept_bcFilter"

df = pd.read_csv(input_file, sep="\t")
unique_readIDs = df['readID'].unique()

for i in range(1, 11):
        # Randomly sample 0.1 million unique readIDs
        sampled_readID = np.random.choice(unique_readIDs, size=100000, replace=False)

        # Filter the DataFrame to include only rows with the sampled readIDs
        filtered_df = df[df['readID'].isin(sampled_readID)]

        # Part 1: subsample read, count exons per gene
        # Group by 'anno_geneID' and count unique 'anno_exonblock'
        #anno_exonblock_count = filtered_df.groupby('anno_geneID')['anno_exonblock'].nunique().reset_index()
        #anno_exonblock_count.columns = ['anno_geneID', 'unique_anno_exonblock_count']
        #output_filename = f"{output_prefix}_exonPerOntargetGene_sub100k_rep{i}.tsv"
        #anno_exonblock_count.to_csv(output_filename, header=True, index=False, sep='\t')

        # Part 2: subsample read, count exons per read
        exon_per_read=filtered_df.groupby('readID')['allinfo_readID'].nunique().reset_index()
        exon_per_read.columns = ['read_id', 'exon_per_read']
        exon_per_read['sample'] = output_prefix
        exon_per_read.to_csv(f"{output_prefix}_exonPerRead_sub100k_rep{i}.tsv", header=True, index=False, sep='\t')
```

**Figure 2A** 

```r
source("scripts/general/R_plot_functions.R")

plotDir = glue("{plotRoot}/size_selection_cutoff1/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
base_path <- "/home/xil4009/scratch_tilgnerlab/size_selection_compare_M8_cutoff1"

#### Box plot #######
all_data <- list()

# Loop through the file pairs and read the data
for (i in 1:10) {
  accept_df <- read.table(glue("{base_path}/VIS1_SYNnND_prom_accept_bcFilter_exonPerRead_sub100k_rep{i}.tsv"), sep = '\t', header = TRUE)
  control_df <- read.table(glue("{base_path}/SC_P56_M8_VIS_exonPerRead_sub100k_rep{i}.tsv"), sep = '\t', header = TRUE)
  label <- paste("Rep", i, sep = "")

  # Add label to data frames and combine
  accept_df$Type <- "Size-selected"
  control_df$Type <- "Control"
  accept_df$Label <- label
  control_df$Label <- label

  all_data[[i]] <- rbind(accept_df, control_df)
}

# Combine all data into a single data frame
final_df <- do.call(rbind, all_data)

# Convert 'Label' to factor with levels in the desired order
final_df$Label <- factor(final_df$Label, levels = c(paste0("Rep", 1:10)))

title="VIS1_SYNnND_prom_accept_sub100k"

colors = c("#efba7a", "#f1776b") # red orange

# vertical
p <- ggplot(final_df, aes(x = exon_per_read, y = Label, fill = Type, color=Type)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.6), size = 0.6, alpha = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +  # Matching outline colors
  labs(title = glue("Exon count per on-target read\n(subsampled 100k)"), x = "Number of exons per read", y = NULL) +
  xlim(0, 20) +
  theme_jennie(plot_title_size = 9, base_size = 9) +
  theme(legend.title = element_blank(), legend.position = "top")

save_png(p, glue("{plotDir}/size_selection_exon_per_read_{title}"), 3, 7)
save_pdf(p, glue("{plotDir}/size_selection_exon_per_read_{title}"), 3, 7)

# horizontal
p <- ggplot(final_df, aes(x = Label, y = exon_per_read, fill=Type)) +
  geom_boxplot(width=0.5, position = position_dodge(width = 0.6), size = 0.6, outlier.shape = NA) + 
  # scale_fill_manual(values=colors) +
  labs(title = glue("Exon number per on-target read {title} (biccn vs junction_pc)"), x = NULL, y = "Exon number per read") +
  ylim(0,20) +
  theme_jennie(plot_title_size=14, base_size=12)
```
    

## Analysis of gene-level enrichment ratio

To track enrichment dynamics, reads in the enriched and control pools were grouped into time bins based on their sequencing timestamps. For each bin, the number of on-target or non-target reads was calculated. The cumulative read counts over time are shown in Figure 2B–F. 

The overall enrichment ratio for each sequencing run is shown in Figure 2G, and was calculated as:
Overall enrichment ratio =  (total number of on-target spliced reads in enriched pool)/(total number of on-target spliced reads in control pool)

The time-bin enrichment ratios are plotted in Figure 3A–B, and were calculated as:
Enrichment ratio of time bin t = (number of on-target spliced reads in time bin t in enriched pool)/(number of on-target spliced reads in time bin t in control  pool)

The enrichment ratio for each gene was computed.  Genes were then classified into three categories: (1) higher enrichment ratio in the enriched pool, (2) higher enrichment ratio in the control pool, or (3) equal in both pools. The distribution of genes across these categories is summarized in Figure 4A. 
Enrichment ratio per gene = (number of on-target spliced reads assigned to gene A in enriched pool)/(number of on-target spliced reads assigned to gene A in control pool)


### Timestamp reads

Split read IDs by time bins

Use: `fastq.gz` `readfish.tsv` `scripts/general/read_timestamp.py`

```bash
cd rootDir

mkdir timestamp
cd timestamp

conda activate isoquant
script=scripts/general/read_timestamp.py
file=readfish.tsv

############## accept

# time stamp each read id
#fastq=../fastq_by_condition/accept/accept.fastq.gz
#python $script --mode read_timestamp --fastq $fastq --tsv $file --sample accept

# assign read ids by time bin
python $script --mode hours_time_bin \
               --read_timestamp accept.readID_timestamp_readLength.tsv --sample accept --hour 1

############## control
#fastq=../fastq_by_condition/control/control.fastq.gz
#python $script --mode read_timestamp --fastq $fastq --tsv $file --sample control

# assign read ids by time bin
python $script --mode hours_time_bin \
               --read_timestamp control.readID_timestamp_readLength.tsv --sample control --hour 1

############## reject
fastq=../fastq_by_condition/reject/reject.fastq.gz
python $script --mode read_timestamp --fastq $fastq --tsv $file --sample reject

# assign read ids by time bin
python $script --mode hours_time_bin \
               --read_timestamp reject.readID_timestamp_readLength.tsv --sample reject --hour 1
```

Use: 
`scripts/scisorseq_scripts/scisor_timebin_LRPOutput.sh` 
`scripts/scisorseq_scripts/scisor_timebin_allinfo.sh`
`scripts/scisorseq_scripts/summary_allReadTypes_ontarget_readsPerTimebin.py`

```bash
x="1" # time bin is 1

# readid to timebin
mkdir allinfo_timebin && cd $_

# skip header, get sorted file with all readID - timebin
cat ../accept_readID_timestamp_readLength_${x}hTimeBin | tail -n +1 | cut -f1,4,5 | sort -u > **accept_readid_${x}htimebin_sorted**
cat ../control_readID_timestamp_readLength_${x}hTimeBin | tail -n +1 | cut -f1,4,5 | sort -u > **control_readid_${x}htimebin_sorted**

######################### Use Scisorseq Allinfo File
declare -A pairs
pairs=(
  ["P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"]="/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/**real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes**"
  ["P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"]="/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/**P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly**"
)

########################
conda activate isoquant 

for folder in "${!pairs[@]}"; do
    cd "$folder/timestamp/allinfo_timebin/"   # Change directory
    probe=${pairs[$folder]}
    echo "folder = $folder"
    echo "target list = $probe"

    # LRP
    scripts/scisorseq_scripts/scisor_timebin_LRPOutput.sh $probe $x

    # ALLINFO
    wc -l *_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt
    scripts/scisorseq_scripts/scisor_timebin_allinfo.sh $probe $x
    wc -l *_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt

    # summarize
    python scripts/scisorseq_scripts/summary_allReadTypes_ontarget_readsPerTimebin.py

    cd ../../../  # Move back to the parent directory
done

# check accept total on-target read count 
awk '{sum += $3} END {print sum}'
awk '{sum += $6} END {print sum}' 
```


### Figure2 B–F: cumulative on-target reads

Use: function `ontarget_cumulative_count_accept_control`

```r
conda activate scisorATAC 

source("scripts/general/R_plot_functions.R")

plotDir = glue("{plotRoot}/ontarget_read_count_over_time_allinfo_cutoff1/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

exp_list = c(
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

source("scripts/general/R_plot_functions.R")
plot_list <- list()
for (i in 1:length(exp_list)) {
        p <- ontarget_cumulative_count_accept_control(exp_list[i], names[i])
        plot_list[[i]] <- p
}

p <- wrap_plots(plot_list, ncol = 1)
prefix = glue("{plotDir}/stacked_bar_pie_chart_all_0326")
save_png(p, prefix, 6.5, 3.5 * length(exp_list))
save_pdf(p, prefix, 6.5, 3.5 * length(exp_list))
```


### Figure 2G: Overall enrich ratio of runs 

3x3 heatmap:

```r
source("scripts/general/R_plot_functions.R")

plotDir = glue("{plotRoot}/enrich_ratio_heatmap_allinfo_cutoff1/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

data <- data.frame(
  Gene_set = c("LA-SYNnND", "SYNnND", "SYN"),
  VIS_1 = c(NA, 1.43, 1.82),
  VIS_2 = c(1.39, 1.82, 1.55),
  HPC = c(1.89, 1.57, NA)
)

library(reshape2)
data_long <- melt(data, id.vars = "Gene_set", variable.name = "Sample", value.name = "Value")

# Ensure the desired row order is preserved
data_long$Gene_set <- factor(data_long$Gene_set, levels = c("SYN", "SYNnND", "LA-SYNnND"))

tile_line_width = 3.5
tile_line_color = "transparent" 

p <- ggplot(data_long, aes(x = Sample, y = Gene_set, fill = Value)) +
	geom_tile(color = tile_line_color, linewidth = tile_line_width) + 
	geom_text(aes(label = ifelse(is.na(Value), "", sprintf("%.2f", Value))), color = "black", size = 8, family = "ArialMT") +
	scale_fill_gradient(
		low = "#e8f7f4",  # Lightest color 
		high = "#91c8d6", # Maximum color
		na.value = "#F9F9F9", 
		name = "Ratio",
		guide = guide_colorbar(  # control colorbar position
			title.position = "top",
			title.hjust = 0.5,
			barwidth = 15,    
			barheight = 0.8,      
			direction = "horizontal"
		)
	) + 
	labs(title = NULL, x = NULL, y = NULL) +
	theme_jennie(base_size=18) +
	scale_x_discrete(expand = c(0, 0)) +  scale_y_discrete(expand = c(0, 0)) +   # Remove gap on y-axis
	theme(
		legend.position = "top",  # 将图例放在顶部
		legend.justification = "center",  # 居中对齐
		axis.line = element_blank(),  # Remove axis lines
		axis.text.x = element_text(face = "bold"),  # Make x-axis labels bold
		axis.text.y = element_text(face = "bold")   # Make y-axis labels bold
	) +
	coord_fixed()

title = "overall_enrich_ratio_3x3"
prefix = glue("{plotDir}/{title}")
save_png(p, prefix, 8, 6)
save_pdf(p, prefix, 8, 6)
```

Error bars:
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
	vis2_df   %>% mutate(Group = "VIS2")
)

summary_df <- all_df %>%
	group_by(Group) %>%
	summarise(
		mean_enrich = mean(EnrichRatio),
		sd_enrich   = sd(EnrichRatio),
		se_enrich   = sd_enrich / sqrt(n()),
		n = n(),
		.groups = "drop"
	) %>%
	mutate( xlab = paste0(Group, "\n", "n = ", n) )

bar_fill <- "steelblue"   
line_color <- "darkred" 

combined_plot <- ggplot(summary_df, aes(x = xlab, y = mean_enrich)) +
	geom_col(width = 0.5, fill = bar_fill, alpha = 0.7) +
	geom_errorbar(
		aes(
			ymin = mean_enrich - sd_enrich,
			ymax = mean_enrich + sd_enrich
		),
		width = 0.25, size = 0.3, color = line_color
	) +
	geom_point(
		data = all_df %>%
			left_join(summary_df %>% select(Group, xlab), by = "Group"),
		aes(x = xlab, y = EnrichRatio),
		position = position_jitter(width = 0.12),
		size = 0.6,
		color = "black",
		alpha = 0.6,
		inherit.aes = FALSE
	) +
	labs( x = "", y = "Enrichment Ratio"
	) +
	scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
	theme_jennie(base_size=12) +
	theme(
		axis.line.x = element_line(size = 0.2),
		axis.line.y = element_line(size = 0.2),
		axis.ticks  = element_line(size = 0.2),
		axis.title.y = element_text(size = 11)
	)

save_png(combined_plot, "R1Q4_enrich_ratio_error_bar", 4, 4)
save_pdf(combined_plot, "R1Q4_enrich_ratio_error_bar", 4, 4)
```



### Figure 3A–B: Enrich ratio over time 

use input `summary_allReadTypes_ontarget_readsPerTimebin`

```r
conda activate scisorATAC 

source("scripts/general/R_plot_functions.R")
base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

#plotDir = glue("{plotRoot}/enrich_ratio_over_time/")
plotDir = glue("{plotRoot}/enrich_ratio_over_time_allinfo_cutoff1_1htimebin")

##### Panel 1: SYNnND

synnnd_exp_list <- c(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1",
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1",
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"
)

names <- c(
    "SYNnND VIS-1 (P)",
    "SYNnND VIS-2 (P)",
    "SYNnND HPC (M)"
)

name2color <- list(
  "SYNnND VIS-1 (P)" = "#1f77b4", # blue M8
  "SYNnND VIS-2 (P)" = "#ff7f0e", # orange M9
  "SYNnND HPC (M)" = "#3cb44b" # green M2
)
name2linetype <- c(
  "SYNnND VIS-1 (P)" = "solid",
  "SYNnND VIS-2 (P)" = "solid",
  "SYNnND HPC (M)" = "solid"
)

df <- format_enrich_ratio_timebin(synnnd_exp_list, names)

# Filter the dataframe to include only the first 50 hours
df <- df %>%
  filter(as.integer(as.character(hour_start)) <= 35) %>%
  mutate(hour_start = as.numeric(as.character(hour_start)))

ymax = 3
x_breaks <- seq(0, 35, by = 2) # Define the breaks you want to show on the x-axis

p1 <- ggplot(df, aes(x = hour_start, y = allinfo_ratio, color = name, linetype = name, group = name)) +
          geom_smooth(se = FALSE, size=0.7) +
          labs(title = NULL, x = "Sequencing time (hours)", y = "Ratio (enriched/control)") +  ylim(0,ymax) +
  scale_color_manual(values = name2color) +
  scale_linetype_manual(values = name2linetype) +
  **scale_x_continuous**(breaks = x_breaks) +  # Specify x-axis breaks
  theme_jennie(base_size=12) +
  theme(legend.title = element_blank(), legend.position = "top",
        legend.key.width = unit(2.5, "cm"), legend.key.height = unit(0.6, "cm"),
            legend.text = element_text(size = 10)  ) +
  guides( color = guide_legend(ncol = 2),  # Arrange the legend in two columns
          linetype = guide_legend(ncol = 2) ) # Arrange the legend in two columns

##### Panel 2: SYN, LA
LA_exp_list <- c(
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1",
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1",
    "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList",
    "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
)

names <- c(
    "SYN VIS-1 (P)",
    "SYN VIS-2 (P)",
    "LA-SYNnND HPC (M)",
    "LA-SYNnND VIS-2 (P)"
)

name2color <- list(
  "SYN VIS-1 (P)" = "#1f77b4", # blue M8
  "SYN VIS-2 (P)" = "#ff7f0e", # orange M9
  "LA-SYNnND VIS-2 (P)" = "#ff7f0e", # orange M9
  "LA-SYNnND HPC (M)" = "#3cb44b" # green M2
)
name2linetype <- c(
  "SYN VIS-1 (P)" = "dashed",
  "SYN VIS-2 (P)" = "dashed",
  "LA-SYNnND HPC (M)" = "solid",
  "LA-SYNnND VIS-2 (P)" = "solid"
)

df <- format_enrich_ratio_timebin(LA_exp_list, names)

# Filter the dataframe to include only the first 35 hours
df <- df %>%
  filter(as.integer(as.character(hour_start)) <= 35) %>%
  mutate(hour_start = as.numeric(as.character(hour_start)))

p2 <- ggplot(df, aes(x = hour_start, y = allinfo_ratio, color = name, linetype = name, group = name)) +
          geom_smooth(se = FALSE, size=0.7) +
          labs(title = NULL, x = "Sequencing time (hours)", y = "Ratio (enriched/control)") +  ylim(0,ymax) +
  scale_color_manual(values = name2color) +
  scale_linetype_manual(values = name2linetype) +
  **scale_x_continuous**(breaks = x_breaks) +  # Specify x-axis breaks
  theme_jennie(base_size=12) +
  theme(legend.title = element_blank(), legend.position = "top",
        legend.key.width = unit(2.5, "cm"), legend.key.height = unit(0.6, "cm"),
            legend.text = element_text(size = 10)  ) +
  guides( color = guide_legend(ncol = 2),  # Arrange the legend in two columns
          linetype = guide_legend(ncol = 2) ) # Arrange the legend in two columns

##### A single multi-panel figure
combined_plot <- (p1 / p2)
title <- ggdraw() + draw_label("Enrich ratio over time", fontface = 'bold', size = 14)
final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.05, 1))

width = 7
height = 4
title = "combined_enrich_ratio_over_time_2panels"
prefix = glue("{plotDir}/{title}")
save_pdf(final_plot, prefix, width, height*2)
save_png(final_plot, prefix, width, height*2)
```

Get ratio mean ± deviation

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
# Filter the dataframe to include only the first x hours
df <- df %>%
  filter(as.integer(as.character(hour_start)) <= 35) %>%
  mutate(hour_start = as.numeric(as.character(hour_start)))
  
  
# get mean and sd in the first 35h 
# get max and min ratio 
for (exp in unique(df$name)) {
print(exp)
tmp = df[df$name == exp, ] %>% select(allinfo_ratio)
print(glue({"mean ± sd: { round(mean(tmp$allinfo_ratio),2) } ± {round(sd(tmp$allinfo_ratio),2)} \n"}))
tmp = df[df$name == exp, ] %>% select(time_bin, allinfo_ratio) %>% slice_max(allinfo_ratio)
print(glue("max: allinfo_ratio: { round(tmp$allinfo_ratio,2) }, time_bin = {tmp$time_bin}"))
tmp = df[df$name == exp, ] %>% select(time_bin, allinfo_ratio) %>% slice_min(allinfo_ratio)
print(glue("min: allinfo_ratio: {round(tmp$allinfo_ratio,2)}, time_bin = {tmp$time_bin}"))
}

[1] "SYNnND VIS-1 (P)"
mean ± sd: 1.37 ± 0.05 
max: allinfo_ratio: 1.49, time_bin = 18
min: allinfo_ratio: 1.26, time_bin = 0
[1] "SYNnND VIS-2 (P)"
mean ± sd: 1.8 ± 0.16 
max: allinfo_ratio: 2.04, time_bin = 6
min: allinfo_ratio: 1.52, time_bin = 26
[1] "SYNnND HPC (M)"
mean ± sd: 1.57 ± 0.1 
max: allinfo_ratio: 1.73, time_bin = 0
min: allinfo_ratio: 1.35, time_bin = 17
[1] "SYN VIS-1 (P)"
mean ± sd: 1.82 ± 0.16 
max: allinfo_ratio: 2.11, time_bin = 32
min: allinfo_ratio: 1.53, time_bin = 1
[1] "SYN VIS-2 (P)"
mean ± sd: 1.58 ± 0.25 
max: allinfo_ratio: 2.29, time_bin = 9
min: allinfo_ratio: 1.22, time_bin = 35
[1] "LA-SYNnND HPC (M)"
mean ± sd: 1.94 ± 0.32 
max: allinfo_ratio: 2.74, time_bin = 6
min: allinfo_ratio: 1.46, time_bin = 25
[1] "LA-SYNnND VIS-2 (P)"
mean ± sd: 1.41 ± 0.08 
max: allinfo_ratio: 1.61, time_bin = 2
min: allinfo_ratio: 1.19, time_bin = 35
```


### Fig 4A: enrich ratio per gene (stacked bars)

Use input `AllInfo`

create a file mapping AllInfo ontarget gene to read number

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
  cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/analysis_target_gene_coverage/accept_bcFilter/
  sample=accept_bcFilter
  cat AllInfo_Ontarget_Gene2reads.txt | cut -f1 | sort | uniq -c | tr -s " " | awk -v sample="$sample" -F" " '{print sample "\t" $1 "\t" $2}' > "${sample}_AllInfo_Ontarget_Gene2reads_freq.txt"

  cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/analysis_target_gene_coverage/control_bcFilter/
  sample=control_bcFilter
  cat  AllInfo_Ontarget_Gene2reads.txt | cut -f1 | sort | uniq -c | tr -s " " | awk -v sample="$sample" -F" " '{print sample "\t" $1 "\t" $2}' > "${sample}_AllInfo_Ontarget_Gene2reads_freq.txt"
done 
```

Frequency statistics

```bash
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/genewise_enrich_ratio_tables
base_path="/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets"

declare -A pairs
pairs=(
  ["P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
  
  ["P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
  ["P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"]="/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes"
    
  ["P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"]="/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly"
)

for exp in "${!pairs[@]}"; do
targets=${pairs[$exp]}

# use allinfo，calculate Gene2Reads
accept_Gene2Reads="/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/analysis_target_gene_coverage/accept_bcFilter/accept_bcFilter_AllInfo_Ontarget_Gene2reads_freq.txt"
control_Gene2Reads="/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/analysis_target_gene_coverage/control_bcFilter/control_bcFilter_AllInfo_Ontarget_Gene2reads_freq.txt"

# merge by gene ID, fill in missing read count
awk -F"\t" -v control_file="$control_Gene2Reads" -v accept_file="$accept_Gene2Reads" '
BEGIN {
    # Read control file and store counts
    while ((getline < control_file) > 0) {
        clt_count[$3] = $2
    }
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
' $targets | awk -F"\t" '{if ($3==0 || $5==0) {print $0 "\t" 0} else {print $0 "\t" $3/$5} }' > geneID_acceptReads_cltReads_enrichRatio_${exp}

done
```

Plot: how many genes have ratio > 1, = 1, < 1
Use function `ratio_per_ontarget_gene_bars`

```r
source("scripts/general/R_plot_functions.R")

plotDir = glue("{plotRoot}/exon_coverage_umicorrect")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

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

p <- ratio_per_ontarget_gene_bars(exp_list,names) # bar plot

title = glue("genewise_enrich_ratio_all")
prefix = glue("{plotDir}/{title}")
save_pdf(p, prefix, 8, 5)
save_png(p, prefix, 8, 5)
```

## Exon-level analyses 

### Figure 4B: Exons per target gene (stacked bars)
For each target gene, we compared the number of unique exons detected in the enriched versus control pool and assigned the gene to one of three categories: (1) more exons in the enriched pool, (2) more in the control pool, or (3) equal numbers. The distribution of genes across these categories is presented in Figure 4B. 

Use function `scripts/general/R_plot_functions.R`

```r
source("scripts/general/R_plot_functions.R")

########### ALLINFO CUTOFF1
plotDir = glue("{plotRoot}/exon_coverage_umicorrect")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

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

source("scripts/general/R_plot_functions.R")
p <- exon_per_ontarget_gene_bars(exp_list,names) 

title = glue("exon_per_target_gene_all")
prefix = glue("{plotDir}/{title}")
save_png(p, prefix, 8, 5)
save_pdf(p, prefix, 8, 5)
```

### Figure 4C, Figure S10-12  Reads per exon 
Exonic coordinates for spliced, on‑target reads were extracted from the “AllInfo” file. These coordinates were intersected with GENCODE vM34 exon annotations using bedtools intersect (v2.26.0), yielding per‑exon read counts. The relationship between exon‑level counts in the enriched and control pools is shown in the scatter plot of Figure 4C.

#### Get exon coverage 

Use: 
Input UMI-corrected AllInfo `AllInfo_IncompleteReads.filtered.corrected.gz`
Script `scripts/scisorseq_scripts/exon_coverage.py`

```bash
cd /athena/tilgnerlab/store/xil4009/Exome_Enrich_Major_Datasets
# Attention: use allinfo umi corrected

exp_list=(
  "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
  "P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"
  "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
  "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
  "P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"
  "P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"
  "P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"
)

##### bedtools intersect allinfo & annotation
for folder in "${exp_list[@]}"; do
    mkdir -p "$folder/exon_coverage/"
    cd "$folder/exon_coverage/"   # Change directory

    anno_exons=/home/xil4009/scratch_tilgnerlab/annotation_exon_count/mouse_vM34_exon_count/anno.exon.sorted.bed # sorted annotation exon blocks

    samples=("accept_bcFilter" "control_bcFilter") 
    for sample in "${samples[@]}"; do
        allinfo=../scisorseqr/${sample}/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz

        # make allinfo exon coords into sorted bed format
        # column: chr start stop readID_exonblock geneID strand
        zcat "$allinfo" | awk '{ n=split($7,a,";%;"); for(i=2;i<=n;i++){print a[i]"\t"$2"\t"$1} }' | sort -u | awk 'OFS="\t" {split($1,a,"_"); print a[1], a[2], a[3], $3"_"$1, $2, a[4]}' | sort -k1,1 -k2,2n | uniq > ${sample}_allinfo_exon_sorted.bed

        echo "perform anno.exon intersect using ${sample}_allinfo_exon_sorted.bed"
        bedtools intersect -a ${sample}_allinfo_exon_sorted.bed -b "$anno_exons" -sorted -wao | sort | uniq > ${sample}_allinfo_exon_intersect_annotation
    done

    cd ../../  # Move back to the parent directory
done

declare -A pairs
pairs=(
  ["P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M8VIS_prom_probe_synaptic_pc_combine_05184e74_a1f68796_20240217_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
  ["P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M9VIS_prom_probe_synaptic_pc_combine_ba275943_d2c25ecb_20240625_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/mouse_probe/synaptic"
  ["P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"]="/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/check_probelmatic_genes/real_geneID_P56M9VIS_prom_probe_junction_0716_accept_allinfo_missed_cutoff10_fineGenes"
  ["P56M2HPC_minion_probe_all_junction_20240125_AllInfoCutoff1"]="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/Mouse.junction.probe.geneID.txt"
  ["P56M2HPC_minion_low_expression_junction_20240202_AllInfoCutoff1_1518FineReadList"]="/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/P56M2HPC_Missed_JunctionProbeGenes_fineReadOnly"
)

# Loop through each pair
for folder in "${!pairs[@]}"; do
    cd "$folder/exon_coverage/"   # Change directory

    targets=${pairs[$folder]}
    echo "folder = $folder"
    echo "target list = $targets"

    script=scripts/scisorseq_scripts/exon_coverage.py

    samples=("accept_bcFilter" "control_bcFilter")

    for sample in "${samples[@]}"; do
        input=${sample}_allinfo_exon_intersect_annotation
        python $script --mode "allinfo_to_annotation_exonblock" --input $input --out_prefix ${sample} --targets $targets
    done

    cd ../../  # Move back to the parent directory
done
```

#### scatter plot 

use `control_bcFilter_exonPerOntargetGene.tsv` `accept_bcFilter_exonPerOntargetGene.tsv`
use `control_bcFilter_readPerOntargetExon.tsv` `accept_bcFilter_readPerOntargetExon.tsv`

reads per target exon, R scatter plot

```bash
source("scripts/general/R_plot_functions.R")

plotDir = glue("{plotRoot}/exon_coverage_umicorrect")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}
base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/"

##### main figures
exp_list <- c("P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1")
names <- c("SYNnND VIS-1")

source("scripts/general/R_plot_functions.R")

plot_list <- list()
for (i in 1:length(exp_list)) {
        p <- exon_coverage_scatter_plot(exp_list[i],names[i]) # function defined in source script
        plot_list[[i]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 1)
prefix = glue("{plotDir}/exon_coverage_main")
save_png(combined_plot, prefix, 11, 5 * length(exp_list))
save_pdf(combined_plot, prefix, 11, 5 * length(exp_list))
```



## Analysis of differential exon inclusion 
Differential exon inclusion between neuronal (ExciteL23, ExciteL4, ExciteL5, ExciteL6, ExciteNP, InhibNeuron) and non‑neuronal (Astrocytes, Microglia, OPCs, DivOPCs, COPs, MFOLs, MOLs) cell types was assessed using a previously described statistical framework [19], which is implemented in ScISOr-ATAC as the “casesVcontrols” function [20]. 

Specifically, prior to assessing statistical significance with Fisher’s exact test, a 2 × 2 contingency table was required to pass a chi-squared criterion for data sufficiency (expected count ≥ 5 in the majority of table cells). Only exons passing this filter were tested, and the resulting P values from all tested exons were corrected for multiple comparisons using the Benjamini–Yekutieli method. 

We compared the enriched and control pools at three sequential levels: all exons, exons passing the chi‑squared criterion, and significant exons (FDR < 0.05). A one‑sided paired t‑test was applied to the paired observations from the three datasets (SYNnND VIS‑1, VIS‑2, LA‑SYNnND VIS‑2) at each level (Figure 4E–G). The number of significant exons per dataset is shown in Figure 4D. Overlaps of significant exons between pools are visualized in Figure 4H. Exon‑usage plots in Figure 4I–J were generated using ScisorWiz [57]. 

### Neuron vs glia

scisorATAC parameter `min_reads`: minimum number of reads for sum of 2 allInfos for a given exon. default = 10

run scisorATAC with MinNumReads = 5

```bash
# SYNnND M8 Prom Junction Allinfo umicorrected 
/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz

exp="P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
sample="P56M8VIS_SYNnND_umicorrect"
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}

conda activate scisorATAC
mkdir scisorATAC_umicorrect && cd $_

###### Separate by cell type
accept_allinfo=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz
control_allinfo=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/AllInfo_IncompleteReads.filtered.corrected.gz
ls -hl $accept_allinfo $control_allinfo

mkdir accept_neuron_VS_nonNeuron_input && cd $_

zcat $accept_allinfo | awk '/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP|InhibNeuron/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}'| gzip -c > ${sample}_accept_AllInfo_neuron.gz
zcat $accept_allinfo | awk '/Astrocytes|Microglia|OPCs|DivOPCs|COPs|MFOLs|MOLs/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | gzip -c > ${sample}_accept_AllInfo_nonNeuron.gz
# 其他celltype：Endo, Ependymal, Fibroblasts, Macrophages, NIPCs, Pericytes

cd ..

mkdir control_neuron_VS_nonNeuron_input && cd $_
zcat $control_allinfo | awk '/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP|InhibNeuron/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}'| gzip -c > ${sample}_control_AllInfo_neuron.gz
zcat $control_allinfo | awk '/Astrocytes|Microglia|OPCs|DivOPCs|COPs|MFOLs|MOLs/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | gzip -c > ${sample}_control_AllInfo_nonNeuron.gz

cd ..

######## Files mouse
PathToChromosomeFile="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/mouse.all.chr.1-20.XY"
PathToAnnotation="/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/gencode.vM34.annotation.gtf.gz"
PathToCellTypeFile="/home/xil4009/store_tilgnerlab/SC_P56_M8_VIS/sc_mouse_celltype_file.txt"

NumThreads=12
ci_low=0.05
ci_high=0.95
MinNumReads=5
OL_fraction=0.8

##### Accept, neuron vs non-neuron
out_path="accept_neuron_VS_nonNeuron"
mkdir $out_path && cd $_

caseListPath=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/caseListPath_accept_neuron
controlListPath=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/controlListPath_accept_nonNeuron

touch $caseListPath $controlListPath
echo /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/${sample}_accept_AllInfo_neuron.gz > $caseListPath
echo /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron_input/${sample}_accept_AllInfo_nonNeuron.gz > $controlListPath

# conda activate scisorATAC
# no scratchLocal on new server 
echo "/tmp/" > tmpdirs; time bash /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/v1.1a_exonInclusion_CTspecific_case_control.sh $caseListPath $controlListPath tmpdirs $PathToChromosomeFile $NumThreads $PathToAnnotation /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/ 0.05 0.95 $MinNumReads zcat 0.8 $PathToCellTypeFile &> report

# R --vanilla --slave --args cases_vs_controls.counts.passed-chi-sq-crit.tab.gz two.sided BY < /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/v1.1a_fisherforNatan.r | gzip -c > cases_vs_controls.counts.pVal.FDR.LOR.tab.gz

cd ..

##### Control, neuron vs non-neuron
out_path="control_neuron_VS_nonNeuron"
mkdir $out_path && cd $_

caseListPath=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/caseListPath_control_neuron
controlListPath=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/controlListPath_control_nonNeuron

touch $caseListPath $controlListPath
echo /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/${sample}_control_AllInfo_neuron.gz > $caseListPath
echo /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron_input/${sample}_control_AllInfo_nonNeuron.gz > $controlListPath

echo "/tmp/" > tmpdirs; time bash /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/v1.1a_exonInclusion_CTspecific_case_control.sh $caseListPath $controlListPath tmpdirs $PathToChromosomeFile $NumThreads $PathToAnnotation /athena/tilgnerlab/scratch/nab4004/snisor-new-probes/cassette-exon-scripts/other-scripts/ 0.05 0.95 $MinNumReads zcat 0.8 $PathToCellTypeFile &> report
```


Count exons in scisorATAC output (reads_above_5, passed_chi, sigFDR)

```bash
probe=Mouse.junction.probe.geneID.txt

#############
## get dPSI of exons pass the min_reads cutoff
# 8 cols
**# col1 : exon, col2: inc-case, col3:exc-case, col4: inc-control, col5: exc-control, col6: sum-case, col7: sum-control, col8: dPSI**
zcat cases_vs_controls.counts.tab.gz | awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$2+$3"\t"$4+$5}' | awk -F"\t" '$6>0 && $7>0 {print $0"\t"($2/$6)-($4/$7)}' | sort -u > cases_vs_controls.counts.dPSI

## get exons of target genes
# add column geneID
cat cases_vs_controls.counts.dPSI | awk -F"\t" 'OFS="\t" {split($1,a,"_"); split(a[4],geneid,"."); print geneid[1],$1,$2,$3,$4,$5,$6,$7,$8 }' | sort | join -t $'\t' -1 1 -2 1 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9 $probe - | tr ' ' '\t' > cases_vs_controls.counts.dPSI.byjunctionGeneID

## get dPSI values for exons tested by chi-sqr
# 11 cols
# $1: exon, $2: inc-case, $3:exc-case, $4: inc-control, $5: exc-control, $6: pval, $7: FDR, $8: LOR, $9: sum-case, $10: sum-control, $11: dPSI
zcat cases_vs_controls.counts.pVal.FDR.LOR.tab.gz | awk -F"\t" '{print $0"\t"$2+$3"\t"$4+$5}' | awk -F"\t" '$9>0 && $10>0 {print $0"\t"($2/$9)-($4/$10)}' | sort -u > cases_vs_controls.counts.pVal.FDR.LOR.dPSI

## get exons of target genes
# add column geneID
cat cases_vs_controls.counts.pVal.FDR.LOR.dPSI | awk 'OFS="\t" {split($1,a,"_"); split(a[4],geneid,"."); print geneid[1],$0}' | sort | join -t $'\t' -1 1 -2 1 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12 $probe - | tr ' ' '\t' > cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID

# FDR is $8
cat cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID | awk '{if ($8<=0.05) {print $0}}' > cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR

cd ..

wc -l */cases_vs_controls.counts.dPSI.byjunctionGeneID # Total exons
wc -l */cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID  # Exons passing chi2 criterion
wc -l */cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR # Significant exons
```


### Fig 4E-G: Total tested, chi2-testable, significant exons

#### Get stats

```r
#exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
#exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
#exp = "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList" 

scisorATAC_stats <- function (enrich_file, control_file, title) {
  enrich_data <- read.table(enrich_file, header = FALSE)
  control_data <- read.table(control_file, header = FALSE)
  # Extract exon IDs
  enrich_ids <- enrich_data[, 2]
  control_ids <- control_data[, 2]
  # Identify exon IDs unique 
  enrich_uniq_ids <- setdiff(enrich_ids, control_ids)
  control_uniq_ids <- setdiff(control_ids, enrich_ids)
  overlap_ids <- intersect(enrich_ids, control_ids)

  message(exp)
  message(title)
  message("Enrich Control Overlap Enrich_unique Control_unique")
  message(paste(length(enrich_ids), length(control_ids),length(overlap_ids), length(enrich_uniq_ids), length(control_uniq_ids), sep=" "))
}

enrich = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.dPSI.byjunctionGeneID")
control = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.dPSI.byjunctionGeneID")
scisorATAC_stats(enrich,control, "Num exon tested in scisorATAC: ")

enrich = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID")
control = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID")
scisorATAC_stats(enrich,control, "Num exon passed chisq criterion (tested): ")

enrich = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
control = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
scisorATAC_stats(enrich,control,"Num exon tested and are sig: ")
```

#### Boxplot + paired t-test

```r
library(ggplot2)
library(ggsignif)
library(reshape2)
library(patchwork) 
library(glue)

source("scripts/general/R_plot_functions.R")

# Define dataframe 
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
    header = "Number of sigFDR exons \n(neuron vs non-neuron)"
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
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5, linewidth = 0.3,
                 position = position_dodge(width = 0.3)) +  
    geom_signif(
      comparisons = list(c("Control", "Enriched")),
      annotations = glue("p = {round(p_value, 3)}"),
      y_position = max(df_melted$value) * 1.1,
      textsize = 3.5,
      tip_length = 0.03,
      size = 0.4,
      vjust = -0.5 
    ) +
	  labs( x = NULL, y = header, fill = NULL, color = "Dataset" ) +
    scale_fill_manual(values = c("Control" = "skyblue", "Enriched" = "#f2a93b")) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    coord_cartesian(ylim = c(0, max(df_melted$value) * 1.2)) +
    theme_jennie(plot_title_size = 10, base_size = 10) +
    guides(fill = "none")  
}

# combine plots 
plot_total <- plot_single(data_list$total$df, data_list$total$header)
plot_chi2 <- plot_single(data_list$chi2_testable$df, data_list$chi2_testable$header)
plot_sigFDR <- plot_single(data_list$sigFDR$df, data_list$sigFDR$header)

combined_plot <- plot_total + plot_chi2 + plot_sigFDR

title="paired_t_test"
prefix = glue("{plotDir}/{title}_combined")
save_png(combined_plot, prefix, 14, 4)
save_pdf(combined_plot, prefix, 14, 4)
```


### Fig 4D: Bar plot significant exons  


```r
source("scripts/general/R_plot_functions.R")

plotDir = glue("{plotRoot}/**scisorATAC_sig_exons_umicorrect**/")
if (!file.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

######## only sigFDR
df_sigFDR <- data.frame( # LA run correct list 
  dataset = c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"),
  Enriched = c(307, 266, 46),
  Control = c(216, 154, 19)
)
df_sigFDR$dataset <- factor(df_sigFDR$dataset, levels=c("SYNnND VIS-1", "SYNnND VIS-2", "LA-SYNnND VIS-2"))

df_melted <- reshape2::melt(df_sigFDR, id.vars = c("dataset"))

colors <- c("Enriched" = "#f2a93b", "Control" = "skyblue")
bar_width <- 0.6

p <- ggplot(df_melted, aes(x = dataset, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
  labs(title = "Number of sigFDR exons \n(neuron vs non-neuron)",
        x = NULL, y = "Exon number", fill = "Pool") +
  geom_text(aes(label = sprintf("%.0f", value)), position = position_dodge(width = bar_width), vjust = -0.5, size = 4, family = "ArialMT") +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0, max(df_melted$value) * 1.1)) +
  scale_y_continuous(expand = c(0, 0)) +
    theme_jennie(plot_title_size=12, base_size=12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

title="neuron_vs_nonNeuron_sigFDR_exons"
prefix = glue("{plotDir}/{title}")
save_png(p, prefix, 6, 5)
save_pdf(p, prefix, 6, 5)
```


### Fig 4H: Venn diagram common/uniq significant exons


```r
#install.packages("VennDiagram")  # Install the package if not already installed
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # disable log

source("scripts/general/R_plot_functions.R")
plotDir = glue("{plotRoot}/scisorATAC_sig_exons_umicorrect/")

### Define
exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
prefix = glue("{plotDir}/venn_SYNnND_VIS1_umi.png")
# myCol <- c("dodgerblue", "goldenrod1")
# 1e90ff  ffbf00
title = "sigFDR exons in enriched & control\n(SYNnND VIS-1)"

exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
prefix = glue("{plotDir}/venn_SYNnND_VIS2_umi.png")
myCol <- c("#efba7a", "#f1776b")
title = "sigFDR exons in enriched & control\n(SYNnND VIS-2)"

exp = "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList"
prefix = glue("{plotDir}/venn_LA-SYNnND_VIS2_umi.png")
myCol <- c("#efba7a", "#f1776b")
title = "sigFDR exons in enriched & control\n(LA-SYNnND VIS-2)"

### Main script
accept = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
control = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")

accept_ids <- as.vector(read.table(accept, header = FALSE)[, 2])
control_ids <- as.vector(read.table(control, header = FALSE)[, 2])

exon_list <- list(
  "Enriched" = accept_ids,
  "Control" = control_ids
)

cat(length(exon_list$Enriched), length(exon_list$Control))

myCol <- c("#F2A93B", "#87CEEB")
venn.diagram(
  x = exon_list,
  filename = prefix,
  output = TRUE,
  imagetype = "png",
  height = 1000 ,
  width = 1000 ,
  resolution = 300,

  # Circles
  lty = 'blank',
  fill = myCol,
  # alpha = 0.5,

  # Numbers
  cex = 0.8,
  fontfamily = "Arial",
  fontface = "plain",

    # Set names
    cat.cex = 0.8,
  cat.fontfamily = "Arial",
  cat.fontface = "plain",
  cat.dist = c(0.07, 0.08), # Modify
  cat.pos = c(-50, 50), # Modify

  margin = 0.2,
  main = title,
  main.cex = 0.8,               # Title size
  main.fontface = "plain",     # Title font face
  main.fontfamily = "Arial"    # Title font family
)
```


## Fig 4I-J: ScisorWIZ 
For all accepted unique exons, rank by dPSI, pick genes for plotting

```r
exp = "P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212_AllInfoCutoff1"
exp = "P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
exp = "P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList" 

accept = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/accept_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")
control = glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/scisorATAC_umicorrect/control_neuron_VS_nonNeuron/cases_vs_controls.counts.pVal.FDR.LOR.dPSI.byjunctionGeneID.sigFDR")

accept_data <- read.table(accept, header = FALSE)
control_data <- read.table(control, header = FALSE)

# Extract exon IDs
accept_ids <- accept_data[, 2]
control_ids <- control_data[, 2]

# Identify exon IDs unique to the accept file
accept_uniq_ids <- setdiff(accept_ids, control_ids)
length(accept_uniq_ids)

# Filter the accept data to keep only the unique exon IDs
accept_uniq_ids_data <- accept_data[accept_data[, 2] %in% accept_uniq_ids, ]

# Rank the unique exon IDs by the value in column 12 (from largest to smallest)
accept_uniq_ids_data <- accept_uniq_ids_data[order(-abs(accept_uniq_ids_data[, 12])), ]

head(accept_uniq_ids_data)

write.table(accept_uniq_ids_data, glue("/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/{exp}/**scisorWiz_umicorrect**/ranked_unique_accept_exons_{exp}.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
```
    

prep allinfo 11 cols

```bash
conda activate scisorATAC
# mkdir **scisorWiz_umicorrect** && cd $_

# M9 Prom LA-SYNnND 
exp="**P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList**" 
sample="P56M9VIS_LA-SYNnND_umicorrect"

# M9 Prom SYNnND 
exp="P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1"
sample="P56M9VIS_SYNnND_umicorrect"

accept_allinfo=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/accept_bcFilter/mmalign/LongReadInfo/**AllInfo_IncompleteReads.filtered.corrected.gz**

###### Rename cell type
ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP -> ExcitatoryNeurons
InhibNeuron -> InhibitoryNeurons
Astrocytes
Microglia
OPCs|DivOPCs -> OPCs
COPs|MFOLs|MOLs -> Oligodendrocytes

zcat $accept_allinfo | awk '{
  gsub(/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP/, "ExcitatoryNeurons");
  gsub(/InhibNeuron/, "InhibitoryNeurons");
  gsub(/OPCs|DivOPCs/, "OPCs");
  gsub(/COPs|MFOLs|MOLs/, "Oligodendrocytes");
  print
}' | gzip -c > ${sample}_accept_AllInfo_rename.gz

zcat ${sample}_accept_AllInfo_rename.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | gzip -c > ${sample}_accept_AllInfo_rename_11cols.gz

control_allinfo=/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/${exp}/scisorseqr/control_bcFilter/mmalign/LongReadInfo/**AllInfo_IncompleteReads.filtered.corrected.gz**

zcat $control_allinfo | awk '{
  gsub(/ExciteL23|ExciteL4|ExciteL5|ExciteL6|ExciteNP/, "ExcitatoryNeurons");
  gsub(/InhibNeuron/, "InhibitoryNeurons");
  gsub(/OPCs|DivOPCs/, "OPCs");
  gsub(/COPs|MFOLs|MOLs/, "Oligodendrocytes");
  print
}' | gzip -c > ${sample}_control_AllInfo_rename.gz

zcat ${sample}_control_AllInfo_rename.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tNoTSS\tNoPolyA\t"$7"\t"$8"\t"$9}' | gzip -c > ${sample}_control_AllInfo_rename_11cols.gz

###### make celltype2color file
```
    
keep only useful read per gene

```bash
# look for cooridnates covering the exon's coord ($9 in 11col allinfo)

# M9 LA
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList/**scisorWiz_umicorrect**

### accept 
zcat P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000027601.14 -v exon=chr3_19262655_19262691 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols_Mtfr1.gz
zcat P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000057789.15 -v exon=chr17_27240637_27240656 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols_Bak1.gz

### control 
zcat P56M9VIS_LA-SYNnND_umicorrect_control_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000027601.14 -v exon=chr3_19262655_19262691 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_LA-SYNnND_umicorrect_control_AllInfo_rename_11cols_Mtfr1.gz

zcat P56M9VIS_LA-SYNnND_umicorrect_control_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000057789.15 -v exon=chr17_27240637_27240656 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_LA-SYNnND_umicorrect_control_AllInfo_rename_11cols_Bak1.gz

# M9 SYNnND 
cd /home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1/**scisorWiz_umicorrect**

# accept 
zcat P56M9VIS_SYNnND_umicorrect_accept_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000037685.16 -v exon=chr5_67913993_67914037 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_SYNnND_umicorrect_accept_AllInfo_rename_11cols_Atp8a1.gz

zcat P56M9VIS_SYNnND_umicorrect_accept_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000031256.12 -v exon=chrX_132973164_132973223 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_SYNnND_umicorrect_accept_AllInfo_rename_11cols_Cstf2.gz

### control 
zcat P56M9VIS_SYNnND_umicorrect_control_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000037685.16 -v exon=chr5_67913993_67914037 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_SYNnND_umicorrect_control_AllInfo_rename_11cols_Atp8a1.gz

zcat P56M9VIS_SYNnND_umicorrect_control_AllInfo_rename_11cols.gz | awk -v gene=ENSMUSG00000031256.12 -v exon=chrX_132973164_132973223 'BEGIN{split(exon,givenExon,"_");}{if($2!=gene){next;} n=split($9,a,";%;"); split(a[2],b,"_"); if(b[2]>givenExon[3]){next;} split(a[n],c,"_"); if(c[3]<givenExon[2]){next;}   print; }' | gzip -c > P56M9VIS_SYNnND_umicorrect_control_AllInfo_rename_11cols_Cstf2.gz
```
    

plot 

```r
conda activate /home/caf4010/miniconda3/envs/R_new

library(ScisorWiz)
library(glue)

###### M9 cutoff1
scisorwiz= "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M9VIS_prom_probe_junction_pc_combine_0185f42f_651fc02d_20240716_AllInfoCutoff1/scisorWiz_umicorrect"
gene.query.list <- c("Atp8a1", "Cstf2")  # , "Trrap"

######## M9 LA-SYNnND 
scisorwiz = "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M9VIS_prom_junction_LACutoff10_fineGene_bf857f38_20240927_1480FinReadList/scisorWiz_umicorrect"
gene.query.list <- c("Mtfr1","Bak1") # "Tbc1d23","Stx3",

########
cTypeFile <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/P56M8VIS_prom_probe_junction_pc_chunk08_3d630332_63h_20240212/scisorWiz/ct2color.txt"
annotation.file <- "/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/gencode.vM34.annotation.gtf.gz"

for (i in 1:length(gene.query.list))
{
gene.query <- gene.query.list[i]

# allInfoFile <- glue("{scisorwiz}/P56M9VIS_LA-SYNnND_umicorrect_accept_AllInfo_rename_11cols_{gene.query}.gz")
allInfoFile <- glue("{scisorwiz}/P56M9VIS_SYNnND_umicorrect_accept_AllInfo_rename_11cols_{gene.query}.gz")

outpath <- glue("{scisorwiz}/accept/")
output.dir1 <- paste0(outpath,gene.query)
ScisorWiz_AllInfo(gencodeAnno = annotation.file,
                  AllInfoInput = allInfoFile,
                  cellTypeFile = cTypeFile,
                  gene = gene.query, cluster = 1, ci = .05,
                  mismatchCutoff = .05,
                  outputDir = output.dir1,
                  interactive = FALSE)

# allInfoFile <- glue("{scisorwiz}/P56M9VIS_LA-SYNnND_umicorrect_control_AllInfo_rename_11cols_{gene.query}.gz")
allInfoFile <- glue("{scisorwiz}/P56M9VIS_SYNnND_umicorrect_control_AllInfo_rename_11cols_{gene.query}.gz")

outpath <- glue("{scisorwiz}/control/")
output.dir1 <- paste0(outpath,gene.query)
ScisorWiz_AllInfo(gencodeAnno = annotation.file,
                  AllInfoInput = allInfoFile,
                  cellTypeFile = cTypeFile,
                  gene = gene.query, cluster = 1, ci = .05,
                  mismatchCutoff = .05,
                  outputDir = output.dir1,
                  interactive = FALSE)
}
```

