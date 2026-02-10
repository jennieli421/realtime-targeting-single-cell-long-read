#! /bin/bash -l

#SBATCH --partition=panda   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=fq_id_length
#SBATCH --time=8:00:00   # HH/MM/SS 
#SBATCH --mem=100G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=xil4009@med.cornell.edu
#SBATCH --mail-type=END 

# Define usage function
usage() {
  echo "Usage: $0 <sample> "
  exit 1
} 

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
  usage
fi

# Parse command line arguments
sample=$1  # accept, control 

echo "read length analysis for sample: $sample"
sample_fq=fastq_by_condition/${sample}/${sample}.fastq.gz

# Retrieve sample_fq file
if [ ! -f "$sample_fq" ]; then
   echo "Cannot find '$sample_fq'. Exiting."
   exit 1
fi

conda activate readfish 

# Get reads length per sample 
outdir="analysis_read_length"
mkdir -p $outdir
# print sample, readid, read length 
zcat "$sample_fq" | awk -v cond="$sample" '{
    if(NR%4==1) {
        id=substr($1, 2)  # Remove the '@' sign
    } 
    if(NR%4==2) {
        print cond "\t" id "\t" length($1)
    }
}' > "${outdir}/${sample}_readID_readLength.txt"

# zcat fastq_by_condition/${sample}/${sample}.fastq.gz | awk -v cond="$sample" '{if(NR%4==1) id=$1; if(NR%4==2) print cond "\t" id "\t" length($1)}' > ${outdir}/${sample}_readID_readLength.txt
# seqkit fx2tab -nl $sample_fq | awk -v cond="$sample" '{print cond "\t" $1 "\t" $2}' > ${outdir}/${sample}_readID_readLength.txt
# zcat fastq_by_condition/${sample}/${sample}.fastq.gz | awk -v cond="$sample" '{if(NR%4==2) print cond "\t" length($1)}' > ${outdir}/${sample}_readLength.txt


echo "=========== get_fastq_readID_readLength terminated `date` ==========="
exit 