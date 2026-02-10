#! /bin/bash -l

#SBATCH --partition=scu-cpu   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=isoquant
#SBATCH --time=120:00:00   # HH/MM/SS 
#SBATCH --mem=256G   # memory requested, units available: K,M,G,T


# Define usage function
usage() {
  echo "Usage: $0 <sampleID> <bcFile> <reference> <genedb>"
  exit 1
} 

# Check for correct number of arguments
if [ "$#" -ne 4 ]; then
  usage
fi

# Parse command line arguments
sample=$1
bcFile=$2
reference=$3
genedb=$4
mode="tenX" 
# mode="double"

conda activate isoquant



# root="$PWD"
root="$PWD/$sample"
if [ ! -d "$root" ]; then
   mkdir -p "$root"
fi

# fqFile="${root}/../../fastq_by_condition/${sample}/*.fastq.gz"

# if [ ! -f "$fqFile" ]; then
#    echo "Cannot find '$fqFile'. Exiting."
#    exit 1
# fi
fqFile=$(find "${root}/../../fastq_by_condition/${sample}/" -maxdepth 1 -type f -name "*.fastq.gz")

if [ -z "$fqFile" ]; then
   echo "Cannot find any files matching '*.fastq.gz' in '${root}/../../fastq_by_condition/${sample}/'. Exiting."
   exit 1
fi

if [ $(echo "$fqFile" | wc -l) -gt 1 ]; then
   echo "More than one file found matching '*.fastq.gz' in '${root}/../../fastq_by_condition/${sample}/'. Exiting."
   exit 1
fi

echo "Crrent sample: ${sample}"
echo "fastq: $fqFile"
echo "barcode: $bcFile"
echo "reference: $reference"
echo "genedb: $genedb" 
echo "mode: $mode"

head -3 $bcFile



# genedb=/athena/tilgnerlab/scratch/caf4010/4_7_2023_Curio/isoquant/gencode.v35.annotation.db
# reference=/athena/tilgnerlab/store/hut2006/data/seq/genomes/H.sapiens/GRCh38/wholeGenomeUnzipped/GRCh38.fa

cd /home/xil4009/store_tilgnerlab/IsoQuant

./isoquant.py --reference $reference \
            --complete_genedb \
            --genedb $genedb \
            --fastq $fqFile \
            --barcode_whitelist $bcFile \
            --mode $mode -d nanopore -t 16 \
            -p $sample -o $root \
            --no_model_construction --count_exons # --no_secondary

# isoquant.py --resume -o <previous output folder>

echo "=========== run_isoquant terminated `date` ==========="

exit