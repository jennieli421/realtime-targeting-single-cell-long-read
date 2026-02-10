#! /bin/bash -l

#SBATCH --partition=scu-cpu   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=isoquant
#SBATCH --time=120:00:00   # HH/MM/SS 
#SBATCH --mem=256G   # memory requested, units available: K,M,G,T

# stats of isoquant output for a sample

# Read sample name and target from the command line
sample=$1
target_geneID2txID=$2

tbl=../${sample}/${sample}/${sample}.read_assignments.tsv 
#read2bc=../${sample}/${sample}/${sample}.barcodes.tsv
read2bc=../${sample}/${sample}/${sample}.barcoded_reads_0.tsv


if [ ! -f "$tbl" ] || [ ! -f "$read2bc" ]; then
   echo "Cannot find '$tbl' or '$read2bc'. Exiting."
   exit 1
fi

## seperate into geneID and pc_txID seperately
cat $target_geneID2txID | awk -F"\t" '{print $1}' | sort -u > target_geneID
cat $target_geneID2txID | awk -F"\t" '{print $2}' | sort -u > target_txID

target_gene=target_geneID
target_tx=target_txID

echo "sample: $sample"
echo "target geneID: $(wc -l $target_gene)"
echo "target transcriptID: $(wc -l $target_tx)"


#########

# Create a file for the output
output_file="isoquant_stats_summary_${sample}.txt"


echo -e "sample \t metric \t $sample" > $output_file # first line

# 1. get barcoded reads
# map readID to barcode
cat $read2bc | awk -F"\t" '$2!="*" && $2!="barcode" {print $1"\t"$2}' | sort -u > readID2barcode

barcoded_reads=$(wc -l < readID2barcode)
echo -e "$sample \t barcoded \t $barcoded_reads" >> $output_file


# read_id	chr	strand	isoform_id	gene_id	assignment_type	assignment_events	exons	additional_info

# 2. filter read_assignment.tsv for geneAssigned reads
# remove .version
cat $tbl | awk -F "\t" 'NR>3 { if ($6!="noninformative" && $6!="intergenic") {split($4,tid,"."); split($5,gid,"."); print $1"\t"$2"\t"$3"\t"tid[1]"\t"gid[1]"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} }' | sort -u > read_assignments.geneAssigned

# map readID to geneID 
cat read_assignments.geneAssigned | awk -F"\t" '{print $1 "\t" $5}' | sort -u > geneAssigned_readID2geneID

# select readID freq==1 
cat geneAssigned_readID2geneID | cut -f1 | sort | uniq -c | awk '{$1=$1; print}' | awk -F" " '{if ($1==1) print $2}' > uniqGeneAssigned_readID
# select readID freq>=2 
cat geneAssigned_readID2geneID | cut -f1 | sort | uniq -c | awk '{$1=$1; print}' | awk -F" " '{if ($1>1) print $2}' > multipleGeneAssigned_readID


echo -e "$sample \t assigned_to_any_gene \t $(cat geneAssigned_readID2geneID | cut -f1 | sort -u | wc -l)" >> $output_file
echo -e "$sample \t assigned_to_multiple_genes \t $(wc -l < multipleGeneAssigned_readID)" >> $output_file
echo -e "$sample \t assigned_to_unique_gene \t $(wc -l < uniqGeneAssigned_readID)" >> $output_file



# 3. filter read_assignment.geneAssigned for uniqGeneAssigned reads (tab delimited)
join -j1 -o1.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9 uniqGeneAssigned_readID read_assignments.geneAssigned | tr ' ' '\t' > read_assignments.uniqGeneAssigned
tbl_uniqGene=read_assignments.uniqGeneAssigned

# 4. get spliced   # | sort -u 
cat $tbl_uniqGene | awk -F"\t" '{ split($8, blocks, ","); if (length(blocks) >= 2) { print $0 } }' > read_assignments.uniqGeneAssigned_spliced 
tbl_spliced=read_assignments.uniqGeneAssigned_spliced 

# get barcoded 
join -j1 -o1.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,1.2 readID2barcode $tbl_uniqGene | tr ' '  '\t' >  read_assignments.uniqGeneAssigned_barcoded
tbl_barcoded=read_assignments.uniqGeneAssigned_barcoded

# get spliced & barcoded 
join -j1 -o1.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,1.2 readID2barcode $tbl_spliced | tr ' '  '\t' >  read_assignments.uniqGeneAssigned_spliced_barcoded
tbl_spliced_barcoded=read_assignments.uniqGeneAssigned_spliced_barcoded

uniq_spliced_reads=$(cut -f1 $tbl_spliced | sort -u | wc -l)
uniq_barcoded_reads=$(cut -f1 $tbl_barcoded | sort -u | wc -l)
uniq_spliced_barcoded_reads=$(cut -f1 $tbl_spliced_barcoded | sort -u | wc -l)

echo -e "$sample \t uniqGeneAssigned_spliced \t $uniq_spliced_reads" >> $output_file
echo -e "$sample \t uniqGeneAssigned_barcoded \t $uniq_barcoded_reads" >> $output_file
echo -e "$sample \t uniqGeneAssigned_spliced_barcoded \t $uniq_spliced_barcoded_reads" >> $output_file
echo -e "" >> $output_file

#########

# get reads on-target by geneID
# awk 'NR==FNR { geneid[$1]; next } { if ($5 in geneid) print $0 }' $target_gene $tbl_barcoded | sort -u > read_assignments.uniqGeneAssigned_barcoded_targetedByGene
# awk 'NR==FNR { geneid[$1]; next } { if ($5 in geneid) print $0 }' $target_gene $tbl_spliced_barcoded | sort -u > read_assignments.uniqGeneAssigned_barcoded_targetedByGene

# sort by col5 'geneID'

sort -k5 $tbl_uniqGene | join -1 1 -2 5 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_gene - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_targetedByGene
sort -k5 $tbl_spliced | join -1 1 -2 5 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_gene - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_spliced_targetedByGene

sort -k5 $tbl_barcoded | join -1 1 -2 5 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_gene - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_barcoded_targetedByGene
sort -k5 $tbl_spliced_barcoded | join -1 1 -2 5 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_gene - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_spliced_barcoded_targetedByGene
tbl_uniqGene_targetByGene=read_assignments.uniqGeneAssigned_targetedByGene
tbl_spliced_targetByGene=read_assignments.uniqGeneAssigned_spliced_targetedByGene
tbl_barcoded_targetByGene=read_assignments.uniqGeneAssigned_barcoded_targetedByGene
tbl_spliced_barcoded_targetByGene=read_assignments.uniqGeneAssigned_spliced_barcoded_targetedByGene

# get reads  on-target by pc_txID
# sort by col4 'txID'   
sort -k4 $tbl_uniqGene | join -1 1 -2 4 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_tx - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_targetedByPCTX
sort -k4 $tbl_spliced | join -1 1 -2 4 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_tx - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_spliced_targetedByPCTX

sort -k4 $tbl_barcoded | join -1 1 -2 4 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_tx - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_barcoded_targetedByPCTX
sort -k4 $tbl_spliced_barcoded | join -1 1 -2 4 -o2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 $target_tx - | tr ' ' '\t' | sort -u > read_assignments.uniqGeneAssigned_spliced_barcoded_targetedByPCTX
tbl_uniqGene_targetByPCTX=read_assignments.uniqGeneAssigned_targetedByPCTX
tbl_spliced_targetByPCTX=read_assignments.uniqGeneAssigned_spliced_targetedByPCTX
tbl_barcoded_targetByPCTX=read_assignments.uniqGeneAssigned_barcoded_targetedByPCTX
tbl_spliced_barcoded_targetByPCTX=read_assignments.uniqGeneAssigned_spliced_barcoded_targetedByPCTX


uniqGene_targetByGene=$(cut -f1 $tbl_uniqGene_targetByGene | sort -u | wc -l)
spliced_targetByGene=$(cut -f1 $tbl_spliced_targetByGene | sort -u | wc -l)
barcoded_targetByGene=$(cut -f1 $tbl_barcoded_targetByGene | sort -u | wc -l)
spliced_barcoded_targetByGene=$(cut -f1 $tbl_spliced_barcoded_targetByGene | sort -u | wc -l)

uniqGene_targetByPCTX=$(cut -f1 $tbl_uniqGene_targetByPCTX | sort -u | wc -l)
spliced_targetByPCTX=$(cut -f1 $tbl_spliced_targetByPCTX | sort -u | wc -l)
barcoded_targetByPCTX=$(cut -f1 $tbl_barcoded_targetByPCTX | sort -u | wc -l)
spliced_barcoded_targetByPCTX=$(cut -f1 $tbl_spliced_barcoded_targetByPCTX | sort -u | wc -l)


echo -e "$sample \t uniqGeneAssign_targetByGene \t $uniqGene_targetByGene" >> $output_file
echo -e "$sample \t spliced_targetByGene \t $spliced_targetByGene" >> $output_file
echo -e "$sample \t barcoded_targetByGene \t $barcoded_targetByGene" >> $output_file
echo -e "$sample \t spliced_barcoded_targetByGene \t $spliced_barcoded_targetByGene" >> $output_file

echo -e "$sample \t uniqGeneAssign_targetByPCTX \t $uniqGene_targetByPCTX" >> $output_file
echo -e "$sample \t spliced_targetByPCTX \t $spliced_targetByPCTX" >> $output_file
echo -e "$sample \t barcoded_targetByPCTX \t $barcoded_targetByPCTX" >> $output_file
echo -e "$sample \t spliced_barcoded_targetByPCTX \t $spliced_barcoded_targetByPCTX" >> $output_file
echo -e "" >> $output_file

genes_targeted_by_barcoded=$(cut -f5 $tbl_barcoded_targetByGene | sort -u | wc -l)
genes_targeted_by_spliced_barcoded=$(cut -f5 $tbl_spliced_barcoded_targetByGene | sort -u | wc -l)
tx_targeted_by_barcoded=$(cut -f4 $tbl_barcoded_targetByPCTX | sort -u | wc -l)
tx_targeted_by_spliced_barcoded=$(cut -f4 $tbl_spliced_barcoded_targetByPCTX | sort -u | wc -l)

echo -e "$sample \t genes_targeted_by_barcoded \t $genes_targeted_by_barcoded/$(wc -l < $target_gene)"  >> $output_file
echo -e "$sample \t genes_targeted_by_spliced_barcoded \t $genes_targeted_by_spliced_barcoded/$(wc -l < $target_gene)"  >> $output_file
echo -e "$sample \t tx_targeted_by_barcoded \t $tx_targeted_by_barcoded/$(wc -l < $target_tx)"  >> $output_file
echo -e "$sample \t tx_targeted_by_spliced_barcoded \t $tx_targeted_by_spliced_barcoded/$(wc -l < $target_tx)"  >> $output_file


#### Plot: frequency of targeted pctxID or geneeID. if no read, 0. #####

# for all targeted geneID, count n_read per Gene ($5)
cat $tbl_barcoded | awk -F"\t" '{print $5 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_barcoded_readsPerGene 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_barcoded_readsPerGene $target_gene > uniqGeneAssigned_barcoded_readsPerGene_targetGenes

cat $tbl_spliced_barcoded | awk -F"\t" '{print $5 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_spliced_barcoded_readsPerGene 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_spliced_barcoded_readsPerGene $target_gene > uniqGeneAssigned_spliced_barcoded_readsPerGene_targetGenes

cat $tbl_uniqGene | awk -F"\t" '{print $5 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_readsPerGene 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_readsPerGene $target_gene > uniqGeneAssigned_readsPerGene_targetGenes

cat $tbl_spliced | awk -F"\t" '{print $5 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_spliced_readsPerGene 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_spliced_readsPerGene $target_gene > uniqGeneAssigned_spliced_readsPerGene_targetGenes

# for all targeted tx, count n_read per Transcript $4
cat $tbl_barcoded | awk -F"\t" '{print $4 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_barcoded_readsPerTranscript 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_barcoded_readsPerTranscript $target_tx > uniqGeneAssigned_barcoded_readsPerTranscript_targetTX

cat $tbl_spliced_barcoded | awk -F"\t" '{print $4 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_spliced_barcoded_readsPerTranscript
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_spliced_barcoded_readsPerTranscript $target_tx > uniqGeneAssigned_spliced_barcoded_readsPerTranscript_targetTX

cat $tbl_uniqGene | awk -F"\t" '{print $4 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_readsPerTranscript 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_readsPerTranscript $target_gene > uniqGeneAssigned_readsPerTranscript_targetTX

cat $tbl_spliced | awk -F"\t" '{print $4 "\t" $1}' | sort -u | awk '{print $1}' | sort | uniq -c | awk '{print $2 "\t" $1}' > uniqGeneAssigned_spliced_readsPerTranscript 
awk 'NR==FNR { count[$1] = $2; next } { print $1 "\t" (count[$1] ? count[$1] : 0) }' uniqGeneAssigned_spliced_readsPerTranscript $target_gene > uniqGeneAssigned_spliced_readsPerTranscript_targetTX


# # count genes with 0 reads 
genes_missed_by_barcoded=$(cut -f2 uniqGeneAssigned_barcoded_readsPerGene_targetGenes | grep ^0 | wc -l)
genes_missed_by_spliced_barcoded=$(cut -f2 uniqGeneAssigned_spliced_barcoded_readsPerGene_targetGenes | grep ^0 | wc -l)
tx_missed_by_barcoded=$(cut -f2 uniqGeneAssigned_barcoded_readsPerTranscript_targetTX | grep ^0 | wc -l)
tx_missed_by_spliced_barcoded=$(cut -f2 uniqGeneAssigned_spliced_barcoded_readsPerTranscript_targetTX | grep ^0 | wc -l)

echo -e "$sample \t genes_missed_by_barcoded \t $genes_missed_by_barcoded" >> $output_file
echo -e "$sample \t genes_missed_by_spliced_barcoded \t $genes_missed_by_spliced_barcoded" >> $output_file
echo -e "$sample \t tx_missed_by_barcoded \t $tx_missed_by_barcoded" >> $output_file
echo -e "$sample \t tx_missed_by_spliced_barcoded \t $tx_missed_by_spliced_barcoded" >> $output_file


#### On-target/spliced & barcoded reads 

echo -e "$sample \t targetByGene_over_spliced_barcoded_reads \t $spliced_barcoded_targetByGene/$uniq_spliced_barcoded_reads = $(echo "scale=2; $spliced_barcoded_targetByGene / $uniq_spliced_barcoded_reads" | bc)" >> $output_file
echo -e "$sample \t targetByPCTX_over_spliced_barcoded_reads \t $spliced_barcoded_targetByPCTX/$uniq_spliced_barcoded_reads = $(echo "scale=2; $spliced_barcoded_targetByPCTX/$uniq_spliced_barcoded_reads" | bc)" >> $output_file

#### On-target / barcoded reads 

echo -e "$sample \t targetByGene_over_barcoded_reads \t $barcoded_targetByGene/$uniq_barcoded_reads = $(echo "scale=2; $barcoded_targetByGene/$uniq_barcoded_reads" | bc)" >> $output_file
echo -e "$sample \t targetByPCTX_over_barcoded_reads \t $barcoded_targetByPCTX/$uniq_barcoded_reads = $(echo "scale=2; $barcoded_targetByPCTX/$uniq_barcoded_reads" | bc)" >> $output_file

echo -e "" >> $output_file

# reads per exon 
spliced_barcoded_readsPerExon=uniqGeneAssigned_spliced_barcoded_readsPerExon

# python /home/xil4009/store_tilgnerlab/my_scripts/isoquant_scripts/isoquant_exon_stats.py "$tbl_barcoded" "$spliced_barcoded_readsPerExon"
python /Users/jennie2000/Desktop/analysis_realtime_targeting/my_scripts_20250311/isoquant_scripts/isoquant_exon_stats.py "$tbl_barcoded" "$spliced_barcoded_readsPerExon"

echo -e "$sample \t exon_targeted_by_spliced_barcoded \t $(cut -f1 $spliced_barcoded_readsPerExon | sort -u | wc -l)" >> $output_file

# python calculate_reads_per_exon.py "$tbl_barcoded"
# python calculate_reads_per_exon.py "$tbl_spliced_barcoded"
# python calculate_reads_per_exon.py "$tbl_spliced_barcoded_targetByGene"
# python calculate_reads_per_exon.py "$spliced_barcoded_targetByPCTX"


echo "Results have been saved to $output_file"


