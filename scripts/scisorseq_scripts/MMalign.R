library(scisorseqr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
filtered_fqFolder <- args[1]

cat(paste("filtered_fqFolder:", filtered_fqFolder, "\n"))


MMalign(fqFolder = filtered_fqFolder,
        mmProgPath = '/home/xil4009/mambaforge/envs/my_r_env/bin/minimap2',
        refGenome='/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/GRCm39.genome.fa', 
        numThreads = 16)

# MMalign(fqFolder = filtered_fqFolder,
#         mmProgPath = '/pbtech_mounts/homes064/weh4002/anaconda3/bin/minimap2',
#         refGenome='/athena/tilgnerlab/scratch/weh4002/Mouse_Ref/GENCODE_GRCm39_vM34/GRCm39.genome.fa', 
#         numThreads = 16)


# mmPath = "/home/caf4010/minimap2-2.21_x64-linux/minimap2"
# reference = "/athena/tilgnerlab/store/hut2006/data/seq/genomes/M.musculus/mm10/chromFa"
# nThreads = 12


# MMalign(fqFolder=filtered_fqFolder, 
#         mmProgPath=mmPath, 
#         refGenome=reference, 
#         numThreads=nThreads)

# human_reference = "/athena/tilgnerlab/scratch/anj2026/2020_09_04_SpeciesComparison_PBMCs/references/H.sapiens/hg38.fa.gz"