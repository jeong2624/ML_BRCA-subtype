"
Created on Thr July 18 2024

@author: Jeong-Woon, Park
"

# Set working directory
setwd("~/OneDrive - inu.ac.kr/SSU_project/논문작성/저널논문_가제본/Code & Result")
source("Step1_Data_collection_preprocessing/Code/Utils.R")
Gene <- fread("Step2_Model_construction/common_gene.csv", data.table = FALSE)$Gene

# Load GSE96058 gene expression data, then transform log2(fpkm + 0.1) to log2(TPM + 1)
GSE96058_exp = fread(file = "Step2_Model_construction/Rawdata/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", data.table = FALSE) %>%
  column_to_rownames("V1") %>%
  fpkm_to_tpm(Pseudocount = 0.1)

GSE96058_exp <- GSE96058_exp[Gene, ]
dim(GSE96058_exp) # [1] 813  3409

# GSE96058 metadata - Illumina Hiseq2000 platform (for building model and test)
GPL11154_metadata = process_GEO_metadata("Step2_Model_construction/Rawdata/GSE96058-GPL11154_series_matrix.txt.gz")

GPL11154_exp = GSE96058_exp[, GPL11154_metadata$title]
dim(GPL11154_exp) # [1] 813  2770

all(colnames(GPL11154_exp) == GPL11154_metadata$title)
colnames(GPL11154_exp) <- GPL11154_metadata$geo_accession
GPL11154_metadata <- GPL11154_metadata[, c("geo_accession", "BRCA_Subtype_PAM50")]
rownames(GPL11154_metadata) <- NULL

# GSE96058 metadata - Illumina Nextseq500 platform
GPL18573_metadata = process_GEO_metadata('Step2_Model_construction/Rawdata/GSE96058-GPL18573_series_matrix.txt.gz')
GPL18573_exp = GSE96058_exp[, GPL18573_metadata$title]
dim(GPL18573_exp) # [1] 813  282

all(colnames(GPL18573_exp) == GPL18573_metadata$title)
colnames(GPL18573_exp) <- GPL18573_metadata$geo_accession
GPL18573_metadata <- GPL18573_metadata[, c("geo_accession", "BRCA_Subtype_PAM50")]
rownames(GPL18573_metadata) <- NULL

# Save gene expression profile and metadata.
fwrite(GPL18573_exp, "Step2_Model_construction/GPL18573_Feature.csv", sep = "\t", row.names = TRUE)
fwrite(GPL18573_metadata, "Step2_Model_construction/GPL18573_PAM50.csv", sep = "\t")
fwrite(GPL11154_exp, "Step2_Model_construction/GPL11154_Feature.csv", sep = "\t", row.names = TRUE)
fwrite(GPL11154_metadata, "Step2_Model_construction/GPL11154_PAM50.csv", sep = "\t")
