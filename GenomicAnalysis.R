install.packages("BioManager")
install.packages("NMF")
install.packages("pheatmap")
BiocManager::install("GenomicDataCommons")
BiocManager::install("maftools")
BiocManager::install("TCGAbiolinks")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("mclust")

library(GenomicDataCommons)
library(maftools)
library(TCGAbiolinks)
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('NMF')
library('pheatmap')
library(mclust)

#Recover LUSC maf file 
lusc.maf<-GDCquery_Maf(tumor="LUSC", save.csv = FALSE, directory = "LUSC", pipelines = "mutect2")
sort(colnames(lusc.maf))
lusc.maf <- read.maf(maf = lusc.maf)

#LUSC
#Shows sample summary.
getSampleSummary(lusc.maf)
#Shows gene summary.
getGeneSummary(lusc.maf)
#shows clinical data associated with samples
getClinicalData(lusc.maf)
#Shows all fields in MAF
getFields(lusc.maf)
#Writes maf summary to an output file with basename lusc.
write.mafSummary(maf = lusc.maf, basename = 'lusc')

# Supplementary Figure 1
plotmafSummary(maf = lusc.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# Gene selection
#oncoplot for frequently mutated genes.(Looking at those who were frequently mutated in more than 15% of the cohort)
oncoplot(maf = lusc.maf, top = 35)

# Supplementary Figure 2
OncogenicPathways(maf = lusc.maf)

# Figure 6
dgi = drugInteractions(maf = luad.maf, fontSize = 0.75)


