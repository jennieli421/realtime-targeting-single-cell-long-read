#! /bin/bash -l

#SBATCH --partition=panda   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=Scisorseq
#SBATCH --time=96:00:00   # HH/MM/SS  
#SBATCH --mem=256G   # memory requested, units available: K,M,G,T  


# run scisorseqr with raw fastq (skip barcode filtering) : MMalign, MapAndFilter, InfoPerLongRead

# Define usage function
usage() {
  echo "Usage: $0 <sampleID>"
  exit 1
} 

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
  usage
fi

# Parse command line arguments
sample=$1

conda activate my_r_env
scripts="/athena/tilgnerlab/store/xil4009/my_scripts/scisorseq_scripts"

root="$PWD/${sample}_noBCFilter"
if [ ! -d "$root" ]; then
   mkdir -p "$root"
fi
fqFolder="${root}/../../fastq_by_condition/${sample}/"

if [ ! -d "$fqFolder" ]; then
   echo "Cannot find '$fqFolder'. Exiting."
   exit 1
fi

echo "Crrent sample: ${sample}_noBCFilter"

# skip barcode filtlering, use raw fastq.gz
filtered_fqFolder=$fqFolder

# mmaglin 
echo "=========== MMalign begin `date`==========="
# Make a output dir for this sample if doesn't exist and cd into it 
cd $root
mmalign="${root}/mmalign/"


if [ -d "$mmalign" ]; then
   echo "Folder 'mmalign' already exists. overwriting."
   rm -r $mmalign 
fi

mkdir -p $mmalign && cd $_

Rscript $scripts/MMalign.R $filtered_fqFolder

echo "=========== MMalign terminated `date` ==========="


# map & filter
echo "=========== MapAndFilter begin ==========="
# cd into the directory with MMoutput and Misc
mmalign="${root}/mmalign/"
cd $mmalign

Rscript ${scripts}/MapAndFilter.R

echo "=========== MapAndFilter terminated `date` ==========="

echo "Generating summary for MapAndFilter outputs" 

cd $mmalign/LRProcessingOutput

files=(
    "mapping.bestperRead.bam" # barcoded mapped 
    "mapping.bestperRead.noRiboRNA.bam" # barcoded mapped  
    "mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz" # mapped spliced
    "newIsoforms_vs_Anno_ignoreAnno/exonStretches.gz" # mapped spliced
)

echo "File Name    Number of Lines" > LRProcessing_${sample}_noBCFilter_summary 
for file in "${files[@]}"; do
   if [[ "$file" == *.bam ]]; then 
      echo "$file    $(samtools view -c $file)" >> LRProcessing_${sample}_noBCFilter_summary  ; 
   else 
      echo "$file    $(zcat $file | wc -l)" >> LRProcessing_${sample}_noBCFilter_summary  ; 
   fi 
done




# # info per long read 
# echo "=========== InfoPerLongRead begin ==========="
# # cd into the mmalign directory containing LRoutput
# cd $mmalign

# # FilteredDeconvBC_AllFiles.csv
# barcodeCSV="${root}/output/OutputFiltered/FilteredDeconvBC_AllFiles.csv"

# # Run R script 
# Rscript ${scripts}/InfoPerLongRead.R $barcodeCSV

# echo "=========== InfoPerLongRead terminated `date` ==========="
# echo "$sample"

exit