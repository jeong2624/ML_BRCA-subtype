"
Created on Thr Aug 15 2024

@author: Jeong-Woon, Park
"

source("Utils.R")

# Define Gene variable, which derived from correlation analysis between gene expression and peak signal in 68 matched TCGA-BRCA patients
Gene <- fread("common_gene.csv", data.table = FALSE)$Gene

# Load GSE96058 gene expression data, then transform log2(fpkm + 0.1) to log2(TPM + 1)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
GSE96058_exp = fread(file = "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", data.table = FALSE) %>%
  column_to_rownames("V1") %>%
  fpkm_to_tpm(Pseudocount = 0.1)

GSE96058_exp <- GSE96058_exp[Gene, ]
dim(GSE96058_exp) # [1] 813  3409

# GSE96058 metadata - Illumina Hiseq2000 platform (for building model and test)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
GPL11154_metadata = process_GEO_metadata("GSE96058-GPL11154_series_matrix.txt.gz")

GPL11154_exp = GSE96058_exp[, GPL11154_metadata$title]
dim(GPL11154_exp) # [1] 813  2770

all(colnames(GPL11154_exp) == GPL11154_metadata$title)
colnames(GPL11154_exp) <- GPL11154_metadata$geo_accession
GPL11154_metadata <- GPL11154_metadata[, c("geo_accession", "BRCA_Subtype_PAM50")]
rownames(GPL11154_metadata) <- NULL

# GSE96058 metadata - Illumina Nextseq500 platform
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
GPL18573_metadata = process_GEO_metadata('GSE96058-GPL18573_series_matrix.txt.gz')
GPL18573_exp = GSE96058_exp[, GPL18573_metadata$title]
dim(GPL18573_exp) # [1] 813  282

all(colnames(GPL18573_exp) == GPL18573_metadata$title)
colnames(GPL18573_exp) <- GPL18573_metadata$geo_accession
GPL18573_metadata <- GPL18573_metadata[, c("geo_accession", "BRCA_Subtype_PAM50")]
rownames(GPL18573_metadata) <- NULL

# Save gene expression profile and metadata.
fwrite(GPL18573_exp, "GPL18573_Feature.csv", sep = "\t", row.names = TRUE)
fwrite(GPL18573_metadata, "GPL18573_PAM50.csv", sep = "\t")
fwrite(GPL11154_exp, "GPL11154_Feature.csv", sep = "\t", row.names = TRUE)
fwrite(GPL11154_metadata, "GPL11154_PAM50.csv", sep = "\t")
