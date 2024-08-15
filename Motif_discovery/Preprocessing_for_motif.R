"
Created on Thr Aug 15 2024

@author: Jeong-Woon, Park
"

# Load the package related to enrichment analysis.
packages = c("tidyverse", "data.table", "rtracklayer")
lapply(packages, library, character.only = TRUE)

# Handling peak and gene information related to the selected 10 genes.
gene_list = c("CDH3", "ERBB2", "TYMS", "GREB1", "OSR1", "MYBL2", "FAM83D", "ESR1", "FOXC1", "NAT1")
peak_info = fread("peak.annotation.csv", data.table = FALSE)
peak_info = peak_info[peak_info$SYMBOL %in% gene_list, ]
peak_info = subset(peak_info, select = c("seqnames", "start", "end", "name", "Score", "strand"))
colnames(peak_info) <- c('chrom', 'chromStart', 'chromEnd', 'name', "score", "strand")
dim(peak_info) # [1] 813  6

peak_info = GenomicRanges::GRanges(peak_info)
rtracklayer::export.bed(object = peak_info, "Step3_Downstream_analysis/Results/Top10_peak.bed")
