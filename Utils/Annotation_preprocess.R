"
Created on Thr Aug 15 2024

@author: Jeong-Woon, Park

"

source("Code/Utils.R")

# Extract gene symbols from the GEO gene expression datasets.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81538
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135298
GSE96058 <- fread("GPL11154-GPL18573_gene.csv")$Gene
GSE81538 <- fread("GSE81538_gene.csv")$Gene
GSE135298 <- fread("Rawdata/GSE135298_gene.csv")$Gene

# Gene annotation preprocessing
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap
gencode.v22 <- fread("gencode.v22.annotation.gene.probeMap", data.table = FALSE) %>%
  filter(!chrom %in% c("chrX", "chrY", "chrM")) %>%
  mutate(gene_length = chromEnd - chromStart + 1) %>%
  arrange(desc(gene_length)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(chrom)

# Extract common genes.
common_gene <- Reduce(intersect, list(gencode.v22$gene, GSE135298, GSE81538, GSE96058))
common_gene <- unique(common_gene)

gencode.v22 <- gencode.v22 %>%
  filter(gene %in% common_gene)

# Peak annotation preprocessing
# https://api.gdc.cancer.gov/data/d28f95fc-af3b-497d-806b-eb0625d6c831
peak <- readPeakFile("brca_promoter.txt")
peak_anno <- annotatePeak(peak, tssRegion = c(-2000, 2000), 
                             TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                             annoDb = "org.Hs.eg.db")

plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno, title = NULL)

df_peak_anno <- as.data.frame(peak_anno@anno) %>%
  filter(!seqnames %in% c("chrX", "chrY", "chrM"),
         str_detect(annotation, "Promoter")) %>%
  arrange(desc(Score)) %>%
  distinct(SYMBOL,.keep_all = TRUE) %>%
  arrange(seqnames)

fwrite(gencode.v22, "gencode.v22.annotation.csv")
fwrite(df_peak_anno, "peak.annotation.csv")
