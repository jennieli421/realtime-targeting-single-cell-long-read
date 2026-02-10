#! /bin/bash -l

#SBATCH --partition=panda   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=coverage_stats
#SBATCH --time=8:00:00   # HH/MM/SS 
#SBATCH --mem=50G   # memory requested, units available: K,M,G,T

# stats of scisorseq output for a sample

# Define usage function
usage() {
  echo "Usage: $0 <sampleID> <probe> "
  exit 1
} 

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  usage
fi


# Parse command line arguments
sample=$1
probe=$2

echo "stats for sample: $sample"
echo "using probe $probe" 

root=$PWD

####
# cutoff=5
####

outdir="${root}/analysis_target_gene_coverage/${sample}"
# outdir="${root}/analysis_target_gene_coverage/${sample}_cutoff${cutoff}"

mkdir -p $outdir

cd $outdir


###### mapped & spliced ######

LRProcessingOutput="${root}/scisorseqr/${sample}/mmalign/LRProcessingOutput"

# Generate transcript to readID list from spliced reads file in scisorseqr output folder ‘LRProcessingOutput’
zcat ${LRProcessingOutput}/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz | awk -F"\t" 'split($2,a,".") {print a[1]"\t"$1}' | sort | uniq > map.best.transcriptWise.gene2reads.txt

# On-target genes and spliced-mapped-barcoded reads 
sort $probe | join -j1 -a1 -o1.1,2.1,2.2 - map.best.transcriptWise.gene2reads.txt | column -t | awk -F" " '$1==$2 {print $2"\t"$3}' | sort | uniq > Ontarget_reads_Genes_map.best.transcriptWise.txt


# Reads number per on-target gene
cat Ontarget_reads_Genes_map.best.transcriptWise.txt | cut -f1 | sort | uniq -c | tr -s " " | awk -v sample="$sample" -F" " '{print sample "\t" $1 "\t" $2}' > "${sample}_Gene2Reads_mapBestTranscriptWise_freq.txt"

# On-target genes (sort)
sort $probe | join -j1 -a1 -o1.1,2.1,2.2 - map.best.transcriptWise.gene2reads.txt | column -t | awk -F" " '$1==$2 {print $2}' | sort | uniq > Targeted_Genes_map.best.transcriptWise.txt

# Missed genes (sort)
sort $probe | join -j1 -a1 -o1.1,2.1,2.2 - map.best.transcriptWise.gene2reads.txt | column -t | awk -F" " '$1!=$2 {print $1}' | sort | uniq > Missed_Genes_map.best.transcriptWise.txt


###### mapped & spliced & allinfo ######
allinfo="${root}/scisorseqr/${sample}/mmalign/LongReadInfo/AllInfo.gz"
# allinfo="${root}/scisorseqr/${sample}/mmalign/LongReadInfo_cutoff${cutoff}/AllInfo.gz"

echo "AllInfo: $allinfo"

zcat $allinfo | awk -F"\t" 'split($2,a,".") {print a[1]"\t"$1}' | sort | uniq > AllInfo.gene2reads.txt
join -j1 -a1 -o1.1,2.1,2.2 $probe AllInfo.gene2reads.txt | column -t | awk -F" " '$1==$2 {print $2"\t"$3}' | sort | uniq > AllInfo_Ontarget_Gene2reads.txt
join -j1 -a1 -o1.1,2.1,2.2 $probe AllInfo.gene2reads.txt | column -t | awk -F" " '$1!=$2 {print $1}' | sort | uniq > AllInfo_Missed_Genes.txt
join -j1 -a1 -o1.1,2.1,2.2 $probe AllInfo.gene2reads.txt | column -t | awk -F" " '$1==$2 {print $1}' | sort | uniq > AllInfo_Ontarget_Genes.txt


# Create a file for the output
output_file="scisorseq_stats_summary_${sample}.txt"
# output_file="scisorseq_stats_summary_${sample}_cutoff${cutoff}.txt"

echo -e "metric \t $sample" > $output_file # header line
# echo -e "metric \t ${sample}_cutoff${cutoff}" > $output_file # header line


echo -e "spliced \t $(wc -l < map.best.transcriptWise.gene2reads.txt)" >> $output_file

# on_target_spliced_mapped_barcoded_reads
echo -e "spliced_targetByGene \t $(wc -l < Ontarget_reads_Genes_map.best.transcriptWise.txt)" >> $output_file

# on-target genes 
echo -e "genes_targeted_by_spliced \t $(wc -l < Targeted_Genes_map.best.transcriptWise.txt)" >> $output_file

# missed genes
echo -e "genes_missed_by_spliced \t $(wc -l < Missed_Genes_map.best.transcriptWise.txt)" >> $output_file

# AllInfo reads 
echo -e "allinfo \t $(wc -l < AllInfo.gene2reads.txt)" >> $output_file

# AllInfo on_target_spliced_mapped_barcoded_reads
echo -e "allinfo_targetByGene \t $(wc -l < AllInfo_Ontarget_Gene2reads.txt)" >> $output_file

# AllInfo on-target genes 
echo -e "genes_targeted_by_allinfo \t $(wc -l < AllInfo_Ontarget_Genes.txt)" >> $output_file

# AllInfo missed genes
echo -e "genes_missed_by_allinfo \t $(wc -l < AllInfo_Missed_Genes.txt)" >> $output_file



echo "=========== coverage_stats terminated `date` ==========="

exit


# ###### mapped ######
# Mapped="${root}/scisorseqr/${sample}/mmalign/Mapped"
# # skip header, print ontarget readID and geneID 
# cat ${Mapped}/${sample}_bed_readID_to_annotation_geneID.tsv | awk -F"\t" 'NR>1 {print $10"\t"$5}' | sort -u > Mapped_Gene2reads.txt
# # skip header, print ontarget readID and geneID 
# cat ${Mapped}/${sample}_Ontarget_bed_readID_to_annotation_geneID.tsv | awk -F"\t" 'NR>1 {print $10"\t"$5}' | sort -u > Mapped_Ontarget_Gene2reads.txt
# # skip header, print reads per ontarget gene 
# cat ${Mapped}/${sample}_readPerOntargetGene.tsv | tail -n +2  > Mapped_Ontarget_Genes_freq.txt



# # mapped reads 
# echo -e "mapped \t $(wc -l < Mapped_Gene2reads.txt)" >> $output_file
# # ontarget mapped reads
# echo -e "mapped_targetByGene \t $(wc -l < Mapped_Ontarget_Gene2reads.txt)" >> $output_file
# # mapped on-target genes 
# echo -e "genes_targeted_by_mapped \t $(wc -l < Mapped_Ontarget_Genes_freq.txt)" >> $output_file
