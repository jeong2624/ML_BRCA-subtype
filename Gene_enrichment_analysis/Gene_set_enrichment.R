"
Created on Thr Aug 15 2024

@author: Jeong-Woon, Park
"

# Load the package related to enrichment analysis.
packages = c("tidyverse", "data.table", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
lapply(packages, library, character.only = TRUE)
options(enrichplot.colours = c("red","blue"))

# Define the symbols of interest
symbol_genes <- c("CDH3", "ERBB2", "TYMS", "GREB1", "OSR1", 
                  "MYBL2", "FAM83D", "ESR1", "FOXC1", "NAT1")
gene_info <- fread("peak.annotation.csv")

# Subset gene_info to match SYMBOL values and create a named vector
names(symbol_genes) <- gene_info[SYMBOL %in% symbol_genes, geneId]

# GO enrichment analysis
GO_BP <- enrichGO(gene = names(symbol_genes),
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff = 1,
                  qvalueCutoff = 0.1,
                  readable = TRUE)

df_GO_BP <- GO_BP@result %>%
  as.data.frame() %>%
  dplyr::select(-p.adjust)

dotplot(GO_BP, x = "GeneRatio", showCategory = 10, 
        color = "qvalue", font.size = 12, label_format = 100,
        title = "GO: Biological process") + aes(shape = I(16)) + 
  aes(color = qvalue) + 
  set_enrichplot_color(type = "color", name = "qvalue")

# KEGG pathway enrichment analysis
KEEG <- enrichKEGG(gene = names(symbol_genes),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.1)

dotplot(KEEG, x = "GeneRatio", showCategory = 10, 
        color = "qvalue", font.size = 12, label_format = 100,
        title = "KEGG pathway") + aes(shape = I(16)) + 
  aes(color = qvalue) + 
  set_enrichplot_color(type = "color", name = "qvalue")
