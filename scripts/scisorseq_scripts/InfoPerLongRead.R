library(scisorseqr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
barcodeCSV <- args[1]

cutoff <- 1

cat(paste("barcodeCSV:", barcodeCSV, "\n"))
cat(paste("minTimesIsoObserve: ", cutoff,"\n"))

InfoPerLongRead(barcodeOutputFile=barcodeCSV, 
                mapAndFilterOut="LRProcessingOutput/", 
                minTimesIsoObserve=cutoff,   #3  # default 5 
                rmTmpFolder=FALSE)

cat("Done. \n")