# Integrative analysis of ATAC-seq and RNA-seq through machine learning identified 10 signature genes for breast cancer subtypes

![그림1](https://github.com/user-attachments/assets/36d82036-780b-45d0-a23e-fc5daefff18c)

Breast cancer is categorized into four main intrinsic subtypes, and distinguishing between these subtypes is crucial for providing personalized treatment to patients. However, systematic analyses exploring the connections between gene expression and chromatin accessibility using bulk RNA-seq and ATAC-seq data, coupled with machine learning algorithms, are lacking. In this study, we developed a classification model based on the integrative analysis of RNA-seq transcriptome and ATAC-seq epigenetic information. We identified 10 signature genes associated with these intrinsic subtypes, which were predominantly linked to immune responses, hormone signaling, cancer progression, and cellular proliferation.

## Data preparation
In this study, four gene expression datasets—GDC TCGA-BRCA, GSE96058, GSE81538, and GSE135298—along with ATAC-seq peak data for GDC TCGA-BRCA were obtained from the UCSC Xena browser (https://xenabrowser.net/datapages/) and the Gene Expression Omnibus (GEO) database (https://www.ncbi.nlm.nih.gov/geo/).

## Data analysis
The R scripts for data integration, preprocessing, and integrative analysis between RNA-seq and ATAC-seq, as well as gene enrichment analysis and motif discovery, and the Python scripts in Jupyter Notebook for feature selection, training, and external validation, along with the input and output files for this analysis, have been uploaded to this repository.

Since some of the data files are limited to 25MB or less, the data used for the analysis has been uploaded as the gz file. Please decompress them before use. So use it after decompressing.

## Installation
All analyses were implemented using R version 4.1.1 and Python 3.9.12 on a MacBook Pro 14 (M1).

### Dependencies
```
scikit-learn==1.5.1
matplotlib==3.9.2
pandas==2.2.2
probatus==3.1.0
numpy==1.26.1
shap==0.43.0
imbalanced-learn==0.12.3
yellowbrick==1.5
umap-learn==0.5.6
```
