#! /bin/bash -l

#SBATCH --partition=panda   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=fq_by_condition
#SBATCH --time=8:00:00   # HH/MM/SS 
#SBATCH --mem=100G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=xil4009@med.cornell.edu
#SBATCH --mail-type=END 

# Define usage function
usage() {
  echo "Usage: $0 <fastq> <sample> "
  exit 1
} 

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  usage
fi

# Parse command line arguments
fastq=$1
sample=$2  # accept, control 

# echo "separating sample: $sample"
# if [ ! -f "$fastq" ]; then
#    echo "Cannot find '$fastq'. Exiting."
#    exit 1
# fi


conda activate readfish 

# Retrieve fastq files of each condition
if [ -d "fastq_by_condition/$sample" ]; then
  ls -hl "fastq_by_condition/$sample" 
  echo "Folder 'fastq_by_condition/$sample' already exists. exiting."
  exit 1
fi

mkdir -p fastq_by_condition/$sample 
echo $fastq

id_list=decision_readID/${sample}_decision_readID.txt
echo $id_list

# # handles fastq or fastq.gz
seqkit grep -f $id_list $fastq -o fastq_by_condition/${sample}/${sample}.fastq.gz



# Get reads length per condition 
outdir="analysis_read_length"
mkdir -p $outdir
zcat fastq_by_condition/${sample}/${sample}.fastq.gz | awk -v cond="$sample" '{if(NR%4==2) print cond "\t" length($1)}' > ${outdir}/${sample}_readLength.txt


echo "=========== fastq_by_condition terminated `date` ==========="
exit 




# condition=$(echo "$sample" | awk -F_ '{print $NF}')

# seqtk subseq $fastq decision_readID/${condition}_decision_readID.txt | gzip -c > fastq_by_condition/${sample}/${sample}.fastq.gz
# Process fastq file based on its compression
# if [[ $fastq == *.gz ]]; then
#     seqtk subseq <(gunzip -c $fastq) decision_readID/${condition}_decision_readID.txt | gzip -c > fastq_by_condition/${sample}/${sample}.fastq.gz
# else
#     seqtk subseq $fastq $id_list | gzip -c > fastq_by_condition/${sample}/${sample}.fastq.gz
# fi