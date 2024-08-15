# 작업 디렉토리 불러오기
setwd("~/OneDrive - inu.ac.kr/SSU_project/논문작성/저널논문_가제본/Code & Result/Step1_Data_collection_preprocessing")
source("Code/Utils.R")

TCGA_BRCA_exp <- fread("Rawdata/TCGA-BRCA.htseq_fpkm.tsv.gz", data.table = FALSE) %>% 
  column_to_rownames("Ensembl_ID") %>%
  fpkm_to_tpm(., Pseudocount = 1) %>%
  as.matrix()

TCGA_BRCA_signal <- fread("Rawdata/brca%2Fbrca_peak_Log2Counts_dedup", data.table = FALSE) %>%
  column_to_rownames("sample") %>%
  as.matrix()

common_sample <- match_sample(TCGA_BRCA_exp, TCGA_BRCA_signal)
match_exp <- TCGA_BRCA_exp[, common_sample]
match_signal <- TCGA_BRCA_signal[, common_sample]

TCGA_BRCA_metadata <- process_BRCA_data(match_exp)
match_exp <- match_exp[, TCGA_BRCA_metadata$patient]
match_signal <- match_signal[, TCGA_BRCA_metadata$patient]

gencode.v22 <- fread("Rawdata/gencode.v22.annotation.csv", data.table = FALSE) %>%
  filter(id %in% rownames(match_exp))

# Peak annotation preprocessing
df_peak_anno <- fread("Rawdata/peak.annotation.csv", data.table = FALSE) %>%
  filter(name %in% rownames(match_signal))

match_exp <- match_exp[gencode.v22$id, ]
rownames(match_exp) <- gencode.v22$gene

match_signal <- match_signal[df_peak_anno$name, ]
rownames(match_signal) <- df_peak_anno$SYMBOL

# Filtering genes with median expression <= 1
keep_gene <- rowMedians(match_exp) > 1
match_exp <- match_exp[keep_gene, ]

# Filtering peaks with median signal <= 1
keep_peak <- rowMedians(match_signal) > 1
match_signal <- match_signal[keep_peak, ]

overlap_genes <- intersect(rownames(match_exp), rownames(match_signal))

correlation_result <- correlation(match_exp, match_signal, TCGA_BRCA_metadata, 
                                  gene_list = overlap_genes,
                                  cor_threshold = 0.65, p_value = 0.01)

# Overlap genes bewteen correlated genes and PAM50 signature genes.
data("pam50")
PAM50_genes <- rownames(pam50$centroids.map)
overlap_pam50 <- intersect(correlation_result$Gene, PAM50_genes)
cor_plot <- correlation_plot(match_exp, match_signal, TCGA_BRCA_metadata, 
                             gene_list = overlap_pam50)
theme_update(text = element_text(size = 12))
pdf("~/Desktop/논문_초본+결과_정리/Integration_결과/PAM50_correlation.pdf")
cor_plot
dev.off()

TCGA_BRCA_external <- TCGA_BRCA_exp[, which(!colnames(TCGA_BRCA_exp) %in% common_sample)]
TCGA_BRCA_external <- TCGA_BRCA_external[gencode.v22$id, ]
rownames(TCGA_BRCA_external) <- gencode.v22$gene
metadata <- process_BRCA_data(TCGA_BRCA_external)
TCGA_BRCA_external <- TCGA_BRCA_external[, which(colnames(TCGA_BRCA_external) %in% metadata$patient)]
metadata <- metadata %>%
  filter(patient %in% colnames(TCGA_BRCA_external))

TCGA_BRCA_external <- TCGA_BRCA_external[correlation_result$Gene, ]
TCGA_BRCA_external <- TCGA_BRCA_external[, metadata$patient]

fwrite(as.data.frame(TCGA_BRCA_external), "Rawdata/TCGA-BRCA.htseq_953_tpm.csv", row.names = TRUE)
fwrite(metadata, "Rawdata/TCGA-BRCA.953_metadata.csv")
fwrite(data.frame(Gene = correlation_result$Gene), "../Step2_Model_construction/common_gene.csv")
