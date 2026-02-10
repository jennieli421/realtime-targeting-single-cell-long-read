#! /bin/bash -l

# Define usage function
usage() {
  echo "Usage: $0 <probe> <x-hour>"
  exit 1
} 

# Parse command line arguments
probe=$1
x=$2  # hours time bin 


##### define sample and file pair 
declare -A pairs
pairs=(
  #["accept_noBCFilter"]="accept_readid_${x}htimebin_sorted"
  #["control_noBCFilter"]="control_readid_${x}htimebin_sorted"
  ["accept_bcFilter"]="accept_readid_${x}htimebin_sorted"
  ["control_bcFilter"]="control_readid_${x}htimebin_sorted"
)

# Loop through each pair
for sample in "${!pairs[@]}"; do
  echo "probe = $probe" 
  echo "x = $x" 
  id2timebin=${pairs[$sample]}
  echo "sample = $sample"
  echo "id2timebin = $id2timebin"
  
  LRProcessingOutput=../../scisorseqr/${sample}/mmalign/LRProcessingOutput
  if [ ! -d "$LRProcessingOutput" ]; then
    echo "Directory $LRProcessingOutput not found. Exiting loop."
    break
  fi

# Generate geneID to readID list from spliced reads file in scisorseqr output folder ‘LRProcessingOutput’
zcat ${LRProcessingOutput}/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz | awk -F"\t" '{split($1,readid,"."); split($2,geneid,"."); print readid[1]"\t"geneid[1]}' | sort -u | join -t $'\t' -1 1 -2 1 -o2.2,1.1,1.2,1.3 $id2timebin - | tr ' ' '\t' | sort -u > ${sample}_spliced_gene2reads_timebin.txt

# On-target, spliced / barcoded+spliced reads 
cat ${sample}_spliced_gene2reads_timebin.txt | join -j1 -a1 -o1.1,2.1,2.2,2.3,2.4 $probe - | column -t | awk -F" " '$1==$2 {print $2"\t"$3"\t"$4"\t"$5}' | sort -u > ${sample}_spliced_Ontarget_gene2reads_timebin.txt

# plot: enrich ratio over time 
cat ${sample}_spliced_Ontarget_gene2reads_timebin.txt | cut -f3,4 | sort | uniq -c |  awk -v sample="$sample" -F" " '{print $2"\t"$3"\t"$1"\t"sample"_splice"}' > ${sample}_spliced_Ontarget_readsPerTimebin.txt

done