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

    allinfo=../../scisorseqr/${sample}/mmalign/LongReadInfo/AllInfo.gz
    if [ ! -f "$allinfo" ]; then
    echo "File $allinfo not found. Exiting loop."
    break
    fi

    # # add timestamp to allinfo
    zcat $allinfo | sort -k1 | join -t $'\t' -1 1 -2 1 -o1.2,1.3,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9 $id2timebin - | tr ' ' '\t' > ${sample}_AllInfo_timebin

    # in allinfo_timebin, exon_chain is col $9 
    cat ${sample}_AllInfo_timebin | awk -F"\t" 'split($4,a,".") {print a[1]"\t"$3"\t"$9"\t"$1"\t"$2}' | sort -u > ${sample}_AllInfo_geneID_readID_exonChain_timebin.txt

    # get on-target genes 
    cat ${sample}_AllInfo_geneID_readID_exonChain_timebin.txt | join -j1 -a1 -o1.1,2.1,2.2,2.3,2.4,2.5 $probe - | column -t | awk -F" " '$1==$2 {print $2"\t"$3"\t"$4"\t"$5"\t"$6}' | sort -u > ${sample}_AllInfo_Ontarget_geneID_readID_exonChain_timebin.txt

    # plot: enrich ratio over time 
    cat ${sample}_AllInfo_Ontarget_geneID_readID_exonChain_timebin.txt | cut -f4,5 | sort | uniq -c |  awk -v sample="$sample" -F" " '{print $2"\t"$3"\t"$1"\t"sample"_splice_allinfo"}' > ${sample}_AllInfo_Ontarget_readsPerTimebin.txt 

    # plot: cumulative read count per gene over time 
    cat ${sample}_AllInfo_Ontarget_geneID_readID_exonChain_timebin.txt | cut -f1,4,5 | sort | uniq -c | awk -v sample="$sample" -F" " '{print $3"\t"$4"\t"$2"\t"$1"\t"sample"_splice_allinfo"}' > "${sample}_AllInfo_Ontarget_readsPerGene_Timebin.txt"
    # geneID  readID  exonChain  timebin  range 

    # get off-target genes 
    cat ${sample}_AllInfo_geneID_readID_exonChain_timebin.txt | join -a1 -o2.1,1.1,1.2,1.3,1.4,1.5 -e $'\t' - $probe | column -t | awk -F" " '$1!=$2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | sort -u | grep -v '^\s*$' > ${sample}_AllInfo_Offntarget_geneID_readID_exonChain_timebin.txt
    cat ${sample}_AllInfo_Offntarget_geneID_readID_exonChain_timebin.txt | cut -f1,4,5 | sort | uniq -c | awk -v sample="$sample" -F" " '{print $3"\t"$4"\t"$2"\t"$1"\t"sample"_splice_allinfo"}' > "${sample}_AllInfo_Offtarget_readsPerGene_Timebin.txt"


done