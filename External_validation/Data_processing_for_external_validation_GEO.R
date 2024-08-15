"
Created on Thr Aug 15 2024

@author: Jeong-Woon, Park
"

source("Utils.R")

# Define Gene variable, which derived from correlation analysis between gene expression and peak signal in 68 matched TCGA-BRCA patients
Gene <- fread("common_gene.csv", data.table = FALSE)$Gene

# Load GSE81538 gene expression data.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81538
GSE81538_exp = fread(file = "GSE81538_gene_expression_405_transformed.csv.gz", data.table = FALSE) %>%
  column_to_rownames("V1") %>%
  fpkm_to_tpm(Pseudocount = 0.1)

GSE81538_exp <- GSE81538_exp[Gene, ]
dim(GSE81538_exp) # [1] 813   405

# Load GSE81538 metadata.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81538
GSE81538_metadata = process_GEO_metadata('GSE81538_series_matrix.txt.gz')
GSE81538_exp = GSE81538_exp[, GSE81538_metadata$title]
dim(GSE81538_exp) # [1] 813  383

all(colnames(GSE81538_exp) == GSE81538_metadata$title)
colnames(GSE81538_exp) <- GSE81538_metadata$geo_accession
GSE81538_metadata <- GSE81538_metadata[, c("geo_accession", "BRCA_Subtype_PAM50")]
rownames(GSE81538_metadata) <- NULL

# Load GSE135298 gene expression data.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135298
GSE135298_exp = fread(file = "GSE135298_OSLO2_EMIT0_RNA-seq.txt.gz", data.table = FALSE) %>%
  column_to_rownames("V1") %>%
  fpkm_to_tpm(Pseudocount = 0.1)

GSE135298_exp <- GSE135298_exp[Gene, ]
dim(GSE135298_exp) # [1] 813    93

# Load the GSE135298 metadata.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135298
GSE135298_metadata = process_GEO_metadata('GSE135298_series_matrix.txt.gz')
GSE135298_exp = GSE135298_exp[, GSE135298_metadata$title]
dim(GSE135298_exp) # [1] 813    82

all(colnames(GSE135298_exp) == GSE135298_metadata$title)
colnames(GSE135298_exp) <- GSE135298_metadata$geo_accession
GSE135298_metadata <- GSE135298_metadata[, c("geo_accession", "BRCA_Subtype_PAM50")]
rownames(GSE135298_metadata) <- NULL

# Save gene expression profile and metadata.
fwrite(GSE81538_exp, "GSE81538_Feature.csv", sep = "\t", row.names = TRUE)
fwrite(GSE81538_metadata, "GSE81538_PAM50.csv", sep = "\t")
fwrite(GSE135298_exp, "GSE135298_Feature.csv", sep = "\t", row.names = TRUE)
fwrite(GSE135298_metadata, "GSE135298_PAM50.csv", sep = "\t")
