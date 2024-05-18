# Load the packages related to this paper.
packages = c("data.table", "tidyverse", "org.Hs.eg.db", 
             "TCGAbiolinks", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene")
lapply(packages, library, character.only = TRUE)
set.seed(2024)

# Convert log2(FPKM + Pseudocount) to log2(TPM + 1)
fpkm_to_tpm <- function(df, Pseudocount) {
  df <- 2^df - Pseudocount
  df <- apply(df, 2, function(x) x / sum(as.numeric(x)) * 10^6)
  df <- log2(df + 1) %>%
    as.data.frame()
}

# Extract matched samples (only -01A samples) between RNA-seq and ATAC-seq.
match_sample <- function(RNA_exp, ATAC_peak){
  common_sample <- intersect(colnames(df_RNA), colnames(df_ATAC))
  common_sample <- common_sample[which(str_detect(common_sample, "-01A"))]
  return(common_sample)
}

# Make metadata except for normal-like breast cancer.
process_BRCA_data <- function(df) {
  TCGA_BRCA <- TCGAquery_subtype("brca")
  
  PAM50 <- TCGA_BRCA[, c("patient", "BRCA_Subtype_PAM50")]
  
  patients <- colnames(df)
  inds <- which(PAM50$patient %in% substr(patients, 1, 12))
  metadata <- as.data.frame(PAM50[inds, ])
  rownames(metadata) <- NULL
  metadata$patient <- paste0(metadata$patient, "-01A")
  
  metadata <- metadata %>%
    filter(BRCA_Subtype_PAM50 != "Normal")
  
  return(metadata)
}

# Spearman correlation analysis bewteen RNA-seq and ATAC-seq at the individual sample level.
correlation <- function(RNA_exp, ATAC_peak, metadata, gene_list, cor_threshold = 0.5, p.value = 0.05){
  cor_result <- data.frame(Gene = character(), r = numeric(), p.value = numeric())
  
  for (gene in gene_list) {
    x <- RNA_exp[gene, metadata$patient]
    y <- ATAC_peak[gene, metadata$patient]
    
    test <- cor.test(x, y, method = "spearman")
    tmp_df <- data.frame(Gene = gene, r = test$estimate, p.value = test$p.value)
    cor_result <- rbind(cor_result, tmp_df)
  }
  
  rownames(cor_result) <- NULL
  cor_result <- na.omit(cor_result)
  cor_result$adj.p.value <- p.adjust(cor_result$p.value, method = "fdr")
  
  keep_genes <- abs(cor_result$r) > cor_threshold & cor_result$adj.p.value < p.value
  cor_result <- cor_result[keep_genes, ] 
  
  return(cor_result)
}





