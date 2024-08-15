# 작업 디렉토리 불러오기
setwd("~/OneDrive - inu.ac.kr/SSU_project/논문작성/저널논문_가제본/Code & Result/Step1_Data_collection_preprocessing")
source("Code/Utils.R")

# 데이터 불러오기 및 전처리
GSE96058 <- fread("Rawdata/GPL11154-GPL18573_gene.csv")$Gene
GSE81538 <- fread("Rawdata/GSE81538_gene.csv")$Gene
GSE135298 <- fread("Rawdata/GSE135298_gene.csv")$Gene

gencode.v22 <- fread("Rawdata/gencode.v22.annotation.gene.probeMap", data.table = FALSE) %>%
  filter(!chrom %in% c("chrX", "chrY", "chrM")) %>%
  mutate(gene_length = chromEnd - chromStart + 1) %>%
  arrange(desc(gene_length)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(chrom)

# 공통 유전자 찾기
common_gene <- Reduce(intersect, list(gencode.v22$gene, GSE135298, GSE81538, GSE96058))
common_gene <- unique(common_gene)

# 공통 유전자만 선택하고 필요없는 열 제거하기
gencode.v22 <- gencode.v22 %>%
  filter(gene %in% common_gene)

# Peak annotation preprocessing
peak <- readPeakFile("Rawdata/brca_promoter.txt")
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

fwrite(gencode.v22, "Rawdata/gencode.v22.annotation.csv")
fwrite(df_peak_anno, "Rawdata/peak.annotation.csv")
