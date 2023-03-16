#In this script we used the data from PRJNA776243 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA776243)
#gene count matrix has no variant surface antigens

library(DESeq2)
library(tidyverse)
library(magrittr)
library(pheatmap)
library(RColorBrewer)

#setting the working directory
dir <- "C:/Users/tifig/OneDrive - Instituto Politécnico do Porto/Ambiente de Trabalho/AcademicLife/R_project_training/Testing_Pf_data"
setwd(dir)
list.files(dir)

metadata <- read.csv(file = "coldata_PRJNA776243")
coldata <- select(metadata, c(2, 7, 10, 15, 26, 28, 29))
count_matrix <- read.table(gzfile("GSE186820_Gene_counts_matrix.txt.gz"))

coldata$condition <- as.factor(coldata$clinical_presentation)
levels(coldata$condition) #check order (it's incorrect, we want control to be uncomplicated malaria)
coldata$condition %<>% relevel("uncomplicated malaria") #use function relevel from magrittr

#Make the deseqdataset using the coldata and the count matrix
dds <- DESeqDataSetFromMatrix(count_matrix, coldata, design = ~ condition)
dds

#Filtering low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
summary(keep) #shows how many columns will be removed, in this case 43
dds.nf <- dds[keep,]

#Transformation of data
round(colSums(counts(dds)) / 1e6, 1) #see library size, numbers of reads widely varies 
rld <- rlog(dds.nf, blind = FALSE) #so its better to use rlog
head(assay(rld), 3)

#See similarity between samples
pcaData <- plotPCA(rld, intgroup = c("condition", "Sample.Description"), returnData = T)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Sample.Description, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with RLOG data")


