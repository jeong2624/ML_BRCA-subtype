# Load the packages related to this paper.
packages = c("data.table", "tidyverse", "org.Hs.eg.db", "GEOquery", "genefu", "ggpubr",
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
  common_sample <- intersect(colnames(RNA_exp), colnames(ATAC_peak))
  common_sample <- common_sample[which(str_detect(common_sample, "-01A"))]
  return(common_sample)
}

# Define a function to retrieve and process except for normal-like breast cancer.
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

process_GEO_metadata <- function(filename) {
    # Get metadata from GEO file
    metadata <- getGEO(filename = filename)

    # Extract relevant information
    sample_name <- metadata$geo_accession
    sample_id <- metadata$title
    
    # Extract PAM50 subtype information
    if ("pam50 subtype:ch1" %in% colnames(metadata@phenoData@data)) {
        PAM50_subtype <- metadata@phenoData@data$`pam50 subtype:ch1`
    } else if ("pam50:ch1" %in% colnames(metadata@phenoData@data)) {
        PAM50_subtype <- metadata@phenoData@data$`pam50:ch1`
    } else {
        print("Check metadata.")
    }
    
    # Combine extracted information into a data frame
    PAM50 <- data.frame(sample_id = sample_id, geo_accession = sample_name, BRCA_Subtype_PAM50 = PAM50_subtype)
    
    # Filter out samples labeled as "Normal"
    PAM50 <- PAM50[PAM50$BRCA_Subtype_PAM50 != "Normal", ]
    
    # Filter out samples with "repl" in sample_id
    PAM50 <- PAM50[!grepl("repl", PAM50$sample_id, ignore.case = TRUE), ]
    
    # Rename columns for clarity
    colnames(PAM50) <- c("title", "geo_accession", "BRCA_Subtype_PAM50")
    
    return(PAM50)
}

# Spearman correlation analysis bewteen RNA-seq and ATAC-seq at the individual sample level.
correlation <- function(RNA_exp, ATAC_peak, metadata, gene_list, cor_threshold = 0.5, p_value = 0.05) {
  
  # Initialize an empty data frame to store results
  cor_result <- data.frame(Gene = character(), r = numeric(), p.value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each gene in gene_list
  for (gene in gene_list) {
    
    # Extract expression data (RNA-seq) for the current gene
    x <- t(RNA_exp[gene, metadata$patient])
    
    # Extract peak data (ATAC-seq) for the current gene
    y <- t(ATAC_peak[gene, metadata$patient])
    
    # Perform Spearman correlation test
    correlation_test <- cor.test(x, y, method = "spearman")
    
    # Store the results in a temporary data frame
    tmp_df <- data.frame(Gene = gene, 
                         r = correlation_test$estimate, 
                         p.value = correlation_test$p.value,
                         stringsAsFactors = FALSE)
    
    # Append the temporary data frame to the results data frame
    cor_result <- rbind(cor_result, tmp_df)
  }
  
  # Filter genes based on correlation coefficient and p-value thresholds
  keep_genes <- cor_result$r >= cor_threshold & cor_result$p.value < p_value
  cor_result <- cor_result[keep_genes, ]
  
  # Remove rows with NA values, if any
  cor_result <- na.omit(cor_result)
  
  # Reset row names
  rownames(cor_result) <- NULL
  
  # Return the final filtered results
  return(cor_result)
}

# Correlation plot
correlation_plot <- function(RNA_exp, ATAC_peak, metadata, gene_list) {
  
  # Initialize an empty list to store plots
  plots <- list()
  
  # Loop through each gene in gene_list
  for (gene in gene_list) {
    
    # Extract expression data (RNA-seq) for the current gene
    x <- as.numeric(RNA_exp[gene, metadata$patient])
    
    # Extract peak data (ATAC-seq) for the current gene
    y <- as.numeric(ATAC_peak[gene, metadata$patient])
    
    # Extract PAM50 subtype from metadata
    PAM50_subtype <- as.factor(metadata$BRCA_Subtype_PAM50)
    
    # Create dataframe between RNA-seq, ATAC-seq, and PAM50 subtype
    df <- data.frame(x = x, y = y, PAM50 = PAM50_subtype)
    
    # Remove NA values
    df <- df[complete.cases(df), ]
    
    # Correlation plot with color by PAM50 subtype
    plot <- ggscatter(data = df, x = "x", y = "y", color = "PAM50",
                      xlab = "Gene expression", ylab = "Chromatin accessibility") +
      geom_smooth(aes(group = 1), method = "lm", color = "blue") +
      stat_cor(aes(group = 1), method = "spearman", cor.coef.name = "R", digits = 3) +
      ggtitle(paste0("Correlation plot for ", gene)) + theme(legend.position = "right")
    
    # Add plot to the list of plots
    plots[[gene]] <- plot
  }
  
  # Return the list of plots
  return(plots)
}
