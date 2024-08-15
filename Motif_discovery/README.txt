"
Created on Sat July 20 2024

@author: Jeong-Woon, Park
"

*** bedtools code
bedtools getfasta -fo Top10_peak.fa -fi hg38.fa.masked -bed top10_peak.bed

*** STREAM code
streme -oc streme_out -p Top10_peak.fa

## MEME-Chip
# Reference : https://github.com/BioinfGuru/memeMotifs
meme-chip -oc meme-chip_out -dna -order 2 -meme-mod zoops -meme-nmotifs 10 -meme-minsites 6 -meme-maxsites 30 -meme-p 4 -db ../Rawdata/JASPAR2024_CORE_non-redundant_pfms_meme.txt Top10_peak.fa -filter-thresh 0.01
