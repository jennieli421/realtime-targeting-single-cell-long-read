#!/bin/bash
# by Hagen Tilgner for mapping of 454/pacBio.MOLECULO reads : 5/2014

function checkForFile {
    file=$1;
    if [ ! -f $file ] 
    then
	echo "ERROR:"$file" is not a regular file ... exiting "	
	exit;
    fi
}

function checkForDirectory {
    dir=$1;
    if [ ! -d $dir ] 
    then
	echo "ERROR:"$dir" is not a regular directory ... exiting "	
	exit;
    fi
}


function getReadNumbersInParallel2 {

    echo "executing getReadNumbersInParallel2"
    
    infileGZ=$1;
    echo "infileGZ="$infileGZ
    if [ ! -f $infileGZ ] 
    then
	echo "ERROR:"$infileGZ" is not a regular file"
	exit;
    fi
    numberLines=`$unzipCommand $infileGZ | wc -l`
    echo $infileGZ" has "$numberLines" lines"
    
    outfile=$2;
    echo "outfile="$outfile
    if [ -e $outfile ]
    then
	echo "ERROR:"$outfile" exists in parallel"
	exit;
    fi

    n=$3;
    echo "n="$n
    tmp=$4;
    echo "tmp="$tmp
    chrs=$5;
    echo "chrs="$chrs
    collectCommand="cat "
    rmCommand="rm "
   
    exons2BeConsidered=$6
    echo "exons2BeConsidered="$exons2BeConsidered

    for c in `cat $chrs`; do  
	for s in - +; do
	    echo -e $c"\t"$s;
	    str="awk -f "${scriptDir}"/v1.1a.collectInclusionExclusionReads_leftAndRightSeparately_fromAllInforFile_sameGeneExc.awk -v chr="$c" -v strand="$s" -v fileGZ="$infileGZ" -v exons2BeConsidered="$exons2BeConsidered" -v unzipCommand="$unzipCommand" > "$tmp"/"$$.$c.$s
	    echo $str >> $tmp"/"$$"parallel.comp.guide";
	    collectCommand=$collectCommand" "$tmp"/"$$.$c.$s
	    rmCommand=$rmCommand" "$tmp"/"$$.$c.$s
	done
    done
    
    time python ${scriptDir}/v0.2.executeInParallel.py --commandFile $tmp"/"$$"parallel.comp.guide" --n $n
    echo $collectCommand;
    echo $rmCommand;
    `$collectCommand > $outfile`;
    gzip $outfile
    `$rmCommand`;
    rm $tmp"/"$$"parallel.comp.guide"

}



# when and where did we do this ? 
#############################
# when and where did we do this ? 
{ echo "#################################";
  echo "# RUNNING [$0]";
#  echo "##### on params [$1]";
  echo "# Current date:`date`";
  echo "# Current dir: `pwd`"; } 1>&2;

###############
# 0. the arguments
echo "+++++++++++++++++++++++++ 1. arguments";

echo "++++++++++++++++++ 1a. flexible";

ALLINFOGZ=$1
echo "ALLINFOGZ="$ALLINFOGZ
checkForFile $ALLINFOGZ

tmpDirFile=$2
echo "tmpDirFile="$tmpDirFile
checkForFile $tmpDirFile

chroms=${3};
echo "chroms="$chroms
checkForFile $chroms

numCPUsHighMemory=$4
echo "numCPUsHighMemory="$numCPUsHighMemory

annotationGZ=$5
echo "annotationGZ="$annotationGZ

scriptDir=${6}
echo "scriptDir="$scriptDir

minPSI=${7}
echo "minPSI="$minPSI

maxPSI=${8}
echo "maxPSI="$maxPSI

minReadNumber=$9
echo "minReadNumber="$minReadNumber

unzipCommand=${10}
echo "unzipCommand="$unzipCommand

minOLfrac=${11}
echo "minOLfrac="$minOLfrac

cellTypeFile=${12}
echo "cellTypeFile="$cellTypeFile


echo "++++++++++++++++++ 1b. deduced from arguments";
tmpdir1=`head -1 $tmpDirFile`
echo "tmpdir1="$tmpdir1
tmpdir1=$tmpdir1"tmp."$$
echo "tmpdir1="$tmpdir1
mkdir $tmpdir1
echo "tmpdir1="$tmpdir1



echo "+++++++++++++++++++++++++ 2. data organization";

echo "++++++++++++++++++ 2a. getting all the genes";

n=`$unzipCommand $ALLINFOGZ | wc -l`;
echo "# total number of reads: "$n;

echo "++++++++++++++++++ 2.b. getting inclusion and exclusion reads";

echo "+++++++++++++++ 2.b_1 getting all exons for all reads";

# $unzipCommand $ALLINFOGZ | awk '{if($2!~/@/ && $2!="none"){print $1"\t"$2"\t"$9;}}' | awk '{n=split($3,a,";%;"); if(n<3){print "ERROR:n="n" for read "$1 > "/dev/stderr";} for(i=3;i<=n-1;i++){print a[i]"\t"$2;}}' | sort | uniq -c | awk '{print $2"\t"$3"\t"$1;}' | gzip -c > all.internalExons.tab.gz

$unzipCommand $ALLINFOGZ | awk '{if($2!~/@/ && $2!="none"){print $1"\t"$2"\t"$7;}}' | awk '{n=split($3,a,";%;"); if(n<3){print "ERROR:n="n" for read "$1 > "/dev/stderr";} for(i=3;i<=n-1;i++){print a[i]"\t"$2;}}' | sort | uniq -c | awk '{print $2"\t"$3"\t"$1;}' | gzip -c > all.internalExons.tab.gz


echo "+++++++++++++++ 2.b_2 getting inclusion and exclusion reads and remving intron retention";
getReadNumbersInParallel2 $ALLINFOGZ sampleBoth.inc.exc.tab $numCPUsHighMemory $tmpdir1 $chroms all.internalExons.tab.gz

$unzipCommand sampleBoth.inc.exc.tab.gz | awk -v anno=$annotationGZ -v unzipCommand=$unzipCommand 'BEGIN{comm=unzipCommand" "anno; while(comm|getline){if($3=="exon"){exon[$1"_"$4"_"$5"_"$7]=1;  for(i=$4;i<=$5;i++){exonic[$1"_"i"_"$7]=1;}     st[$1"_"$4"_"$7]=1; en[$1"_"$5"_"$7]=1;}}}{split($1,a,"_"); x=0; for(i=a[2];i<=a[3];i++){ s=a[1]"_"i"_"a[5]; if(!(s in exonic)){x++;}  }   if(!(a[1]"_"a[2]"_"a[3]"_"a[5] in exon) && x>=70 && (a[1]"_"a[2]"_"a[5] in st || a[1]"_"a[3]"_"a[5] in en)){next;} print $0;}' | gzip -c > sampleBoth.inc.exc.tab.NoRetention.gz

# echo "+++++++++++++++ 2.b_3 getting genes that have >= 2 exons with >=minPSI and <=maxPSI exon inclusion";
# $unzipCommand sampleBoth.inc.exc.tab.NoRetention.gz | awk -v minPSI=$minPSI -v maxPSI=$maxPSI -v minOLfrac=$minOLfrac -v minReadNumber=$minReadNumber '{if($2+$3+$5==0 || $2+$4+$5==0){next;} if($2+$3+$4+$5<minReadNumber){next;} psi=($2+$3+$4)/($2+$3+$4+$5); leftPSI=($2+$3)/($2+$3+$5); rightPSI=($2+$4)/($2+$4+$5); split($1,a,"_"); if(($2+$3+$4+$5)/$6<minOLfrac){next;} if(rightPSI>=minPSI && rightPSI<=maxPSI && leftPSI>=minPSI && leftPSI<=maxPSI && psi>=minPSI && psi<=maxPSI){print $0; }}' | gzip -c > all.altExons.tab.gz

# echo "++++++++++++++++++ 2.c. getting inclusion and exclusion reads for each cell type";
# echo "+++++++++++++++ 2.c_1 in separate output files";

# for CT in `cat $cellTypeFile`; do
#     $unzipCommand $ALLINFOGZ | awk -v CT=$CT '$3==CT' | gzip -c > tmp.gz
#     getReadNumbersInParallel2 tmp.gz sampleBoth.inc.exc.tab.$CT $numCPUsHighMemory $tmpdir1 $chroms all.internalExons.tab.gz
#     #gzip sampleBoth.inc.exc.tab.$CT
#     #$unzipCommand sampleBoth.inc.exc.tab.$CT.gz | awk -v exonsToUse=all.altExons.tab.gz -v unzipCommand=$unzipCommand 'BEGIN{comm=unzipCommand" "exonsToUse; while(comm|getline){use[$1]=1;}}{if($1 in use){print $0;}}' | gzip -c > all.altExons.inc.exc.tab.$CT.gz 
# done
# rm tmp.gz



# echo "+++++++++++++++ 2.c_2 in combing the data";
# awk -v cellTypeFile=$cellTypeFile -v unzipCommand=$unzipCommand -v ASexons=all.altExons.tab.gz -v CTspecificTrunk=sampleBoth.inc.exc.tab -v minReadNumber=$minReadNumber 'BEGIN{comm="cat "cellTypeFile; ctNumber=0; while(comm|getline){ctNumber++; CT[ctNumber]=$1;}  comm=unzipCommand" "ASexons; while(comm|getline){exonsToUse[$1]=1; for(i=1;i<=ctNumber;i++){inc[$1"\t"CT[i]]=0; exc[$1"\t"CT[i]]=0; skip[$1"\t"CT[i]]=0; psi[$1"\t"CT[i]]="NA";}} for(i=1;i<=ctNumber;i++){comm=unzipCommand" "CTspecificTrunk"."CT[i]".gz"; while(comm|getline){if(!($1 in exonsToUse)){continue;} if($2+$3+$4+$5<minReadNumber){continue;} psi[$1"\t"CT[i]]=($2+$3+$4)/($2+$3+$4+$5);} close(comm);}  commentLine="exonID"; for(i=1;i<=ctNumber;i++){commentLine=commentLine"\t"CT[i];} print "#"commentLine; for(ex in exonsToUse){str=ex; for(i=1;i<=ctNumber;i++){str=str" "psi[ex"\t"CT[i]]; }      print str;}}' | gzip -c > all.altExons.matrix.tab.gz

# echo "+++++++++++++++++++++++++ 3. cleaning up ";
# rmdir $tmpdir1; 
# exit;
