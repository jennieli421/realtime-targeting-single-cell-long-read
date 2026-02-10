#! /bin/bash -l

#SBATCH --partition=scu-cpu   # cluster-specific 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=Scisorseq
#SBATCH --time=96:00:00   # HH/MM/SS  
#SBATCH --mem=256G   # memory requested, units available: K,M,G,T  

# Full scisorseqr workflow: GetBarcodes, MMalign, MapAndFilter, InfoPerLongRead.
# barcode decovolution requires exact string matching 


# Define usage function
usage() {
  echo "Usage: $0 <sampleID> <bcFile>"
  exit 1
} 

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  usage
fi

# Parse command line arguments
sample=$1
bcFile=$2

# conda activate my_r_env
conda activate /home/weh4002/anaconda3/envs/seurat
scripts="/athena/tilgnerlab/store/xil4009/my_scripts/scisorseq_scripts"

root="$PWD/${sample}_bcFilter"
if [ ! -d "$root" ]; then
   mkdir -p "$root"
fi

fqFolder="${root}/../../fastq_by_condition/${sample}/"

if [ ! -d "$fqFolder" ]; then
   echo "Cannot find '$fqFolder'. Exiting."
   exit 1
fi

echo "Crrent sample: ${sample}_bcFilter"

# Get barcodes 
echo "=========== GetBarcodes begin `date` ==========="

head -3 $bcFile

cd $root
outputFolder="${root}/output/"

if [ -e "$outputFolder/OutputFiltered/Barcoded_AllFiles.fastq.gz" ]; then
   # 如果output存在，跳过
   echo "Proceed to MMalign. 'OutputFiltered/Barcoded_AllFiles.fastq.gz' already exists."
else
   Rscript ${scripts}/GetBarcodes.R $fqFolder $bcFile $outputFolder
fi

echo "=========== GetBarcodes terminated `date` ==========="


# mmaglin 
echo "=========== MMalign begin `date`==========="
# Make a output dir for this sample if doesn't exist and cd into it 
cd $root
mmalign="${root}/mmalign/"


if [ ! -d "$outputFolder" ]; then
   # 如果output不存在，终止
   echo "Folder '$outputFolder' doesn't exists. Exiting."
   exit 1
fi


if [ -e "$mmalign/MMoutput/Barcoded_AllFiles.bam" ]; then
   # 如果bam存在，跳过
   echo "Proceed to MapAndFilter. 'MMoutput/Barcoded_AllFiles.bam' already exists."
else
   mkdir -p $mmalign && cd $_
   filtered_fqFolder="${outputFolder}/OutputFiltered/"
   Rscript $scripts/MMalign.R $filtered_fqFolder
fi

echo "=========== MMalign terminated `date` ==========="


# map & filter
echo "=========== MapAndFilter begin ==========="
# cd into the directory with MMoutput and Misc
mmalign="${root}/mmalign/"
cd $mmalign

if [ -d "$mmalign/LRProcessingOutput" ]; then
   # 如果有旧的输出，删除LRProcessing
   echo "Folder 'LRProcessingOutput' already exists. deleting."
   rm -r $mmalign/LRProcessingOutput
   rm -r $mmalign/LongReadInfo
fi

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

echo "File Name    Number of Lines" > LRProcessing_${sample}_bcFilter_summary 
for file in "${files[@]}"; do
   if [[ "$file" == *.bam ]]; then 
      echo "$file    $(samtools view -c $file)" >> LRProcessing_${sample}_bcFilter_summary  ; 
   else 
      echo "$file    $(zcat $file | wc -l)" >> LRProcessing_${sample}_bcFilter_summary  ; 
   fi 
done



# info per long read 
echo "=========== InfoPerLongRead begin ==========="
# cd into the mmalign directory containing LRoutput
mmalign="${root}/mmalign/"
cd $mmalign

# FilteredDeconvBC_AllFiles.csv
barcodeCSV="${root}/output/OutputFiltered/FilteredDeconvBC_AllFiles.csv"

# Run R script 
Rscript ${scripts}/InfoPerLongRead.R $barcodeCSV

echo "=========== InfoPerLongRead terminated `date` ==========="
echo "$sample"

exit