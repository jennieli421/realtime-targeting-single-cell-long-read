library(scisorseqr)

# MapAndFilter(numThreads=12, filterFullLength=TRUE, 
#              polyABed ='/athena/tilgnerlab/scratch/anj2026/2020_03_11_LinearPCR_SpikedSequins/04_03_Testing/cagePolyAScripts/atlas.clusters_chr.mm10.2-0.bed.gz',
#              cageBed ='/athena/tilgnerlab/scratch/anj2026/2020_03_11_LinearPCR_SpikedSequins/04_03_Testing/cagePolyAScripts/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz',
#              annoGZ='/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/gencode.vM34.annotation.gtf.gz',
#              seqDir = '/athena/tilgnerlab/store/weh4002/genome/mm39', 
#              cp_distance=50)


MapAndFilter(annoGZ='/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/gencode.vM34.annotation.gtf.gz',
             numThreads=20,
             seqDir = '/athena/tilgnerlab/store/weh4002/genome/mm39',
             filterFullLength=FALSE)





##################
# polyABed = "/athena/tilgnerlab/store/shardwick/genomes/human/atlas.clusters.2.0.GRCh38.96_chrNames.bed.gz"
# cageBed = "/athena/tilgnerlab/scratch/anj2026/2020_11_30_snisorseqAnalysis/2021_04_06_ONTdataProcessing_BothRuns/updatedPhastcons/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
# annotation = "/athena/tilgnerlab/store/hut2006/data/annotations/H.sapiens/GRCh38/gencode.v34.annotation.gtf.gz"
# referenceDir = "/athena/tilgnerlab/store/hut2006/data/seq/genomes/H.sapiens/GRCh38/chromFa/"

# # new order of args 
# MapAndFilter(numThreads=12, filterFullLength=TRUE, 
#             polyABed=polyABed, 
#             cageBed=cageBed, 
#             annoGZ=annotation, 
#             seqDir=referenceDir,  
#             genomeVersion="GRCh38")

# original order of args 
# MapAndFilter(outputDir="LRoutput", 
#             annoGZ=annotation, 
#             seqDir=referenceDir, 
#             cageBed=cageBed,
#             polyABed=polyABed,
#             filterFullLength=TRUE,
#             genomeVersion="GRCh38",
#             numThreads=12 )


cat("Done. \n")




# MapAndFilter(annoGZ='/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/revise_annotation/synnnd_junction_probes_annotation.gtf.gz',
#              numThreads=20,
#              seqDir = '/athena/tilgnerlab/store/weh4002/genome/mm39',
#              filterFullLength=FALSE)