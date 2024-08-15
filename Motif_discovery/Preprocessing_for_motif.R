"
Created on Sta July 20 2024

@author: Jeong-Woon, Park
"

# Set working directory
setwd("~/OneDrive - inu.ac.kr/SSU_project/논문작성/저널논문_가제본/Code & Result")

# Load the package related to enrichment analysis.
packages = c("tidyverse", "data.table", "rtracklayer")
lapply(packages, library, character.only = TRUE)

### 10 genes based on their SHAP values.
gene_list = c("CDH3", "ERBB2", "TYMS", "GREB1", "OSR1", "MYBL2", "FAM83D", "ESR1", "FOXC1", "NAT1")

# Handling peak and gene information related to top 20 genes.
peak_info = fread("Step1_Data_collection_preprocessing/Results/peak.annotation.csv", data.table = FALSE)
peak_info = peak_info[peak_info$SYMBOL %in% gene_list, ]
peak_info = subset(peak_info, select = c("seqnames", "start", "end", "name", "Score", "strand"))
colnames(peak_info) <- c('chrom', 'chromStart', 'chromEnd', 'name', "score", "strand")
dim(peak_info) # [1] 813  6

peak_info = GenomicRanges::GRanges(peak_info)
rtracklayer::export.bed(object = peak_info, "Step3_Downstream_analysis/Results/Top10_peak.bed")

# Save peak information of top20 genes as bed file format.
fwrite(peak_info, "./Step3_Downstream_analysis/Results/813_promoter_peak.bed", sep = "\t")
