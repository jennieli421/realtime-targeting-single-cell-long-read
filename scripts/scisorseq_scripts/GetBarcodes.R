library(scisorseqr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

fqFolder <- args[1]
bcFile <- args[2]
outputFolder <- args[3]

cat(paste("fqFolder:", fqFolder, "\n"))
cat(paste("bcFile:", bcFile, "\n"))
cat(paste("outputFolder:", outputFolder, "\n"))

GetBarcodes(fqFolder=fqFolder, BCClustAssignFile=bcFile, outputFolder=outputFolder, numProcesses=12, filterReads=TRUE)

cat("Done. \n")
