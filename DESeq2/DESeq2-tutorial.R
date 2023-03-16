# This workflow was taken from http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

# tximeta is usually used with salmon
# when using DESeq2 to perform DEG analysis it is necessary to always use count, never use normalized values such as RPKM, FPKM etc

library(airway) # has the data of the experiment
library(tximeta) # reading data and import to DESeq2
library(magrittr)
library(DESeq2) # RNA-seq statistical analysis
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggbeeswarm)
library(apeglm)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db) # the organism annotation package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”), using Entrez Gene IDs (“eg”) 
library(ReportingTools)
library(Gviz) # for plotting GRanges (Plotting fold changes in genomic space)

# Reading in data with tximeta ===================================================================================================================================================

# Identify the directory of the package airway since it contains the data to be used
dir <- system.file("extdata", package = "airway", mustWork = T)

# list all files present in the directory
list.files(dir)

# list files present in the file quants from the directory
list.files(file.path(dir, "quants"))

# sample_table.csv ussually needs to be created for your own project
# its where you introduce information about the samples
csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)
coldata

# this workflow works only with the first two samples 
# so you need to extract them from the table of information
coldata <- coldata[1:2,]
# and than had the name of the samples and path to the quantification (BAM) file
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")


# tximeta imports data at the transcript level
se <- tximeta(coldata) #se stands for summarized experiment

# summarize the transcript-level quantifications to the gene level 
gse <- summarizeToGene(se)

dim(se)
dim(gse)
head(row.names(se))
head(row.names(gse))
# compare the dimensions of se and gse and the name of the rows
# they passed to gene ids

# Reading data in DESeq2 from SummarizedExperiment =====================================================================================================================

# Now that we have gene associated with counts we have to use DESeq2
# To use DESeq2 we need to create a DESeqDataSet object from the gse and se objects

data(gse) # Loads all samples provided in the airway package
gse # check what is the SummarizedExperiment object

assayNames(gse)
# gse is composed of 3 matrices (counts, abundance and length)
# counts =  estimated fragment counts for each gene and sample
# abundance = estimated transcript abundances in TPM
# length =  effective gene lengths which include changes in length due to biases as well as due to transcript usage

head(assay(gse),3)
# function assay pulls out the first matrix (counts)

colSums(assay(gse))
# function shows the amount of reads mapped for each sample (it sums all the counts)

rowRanges(gse)
seqinfo(rowRanges(gse))
# shows info about the genes
# rowRanges is the GRanges of the genes from the left-most position of all the transcripts to the right-most position of all the transcripts

colData(gse)
# reflects the data.frame that was provided to the tximeta function for importing the quantification data
# it has info about sample name, donor iD and treatment conditions

# From this point on data could be analysed using several statistical packages for DEG analysis
# This protocol uses DESeq2

# In your case where you dont have replicates you should use EdgeR and be really hard on the filtering!!
# To prevent high fold changes due to low read counts

# DESeq2 uses the DESeqDataSet object which is built on top of the SummarisedExperiment object
# DESeqDataSet has an associated design formula. 
# Experimental design has to be specified at the begging of the analysis

# For that first examine the columns of colData(gse)
gse$donor
gse$condition

gse$cell <- gse$donor #renaming the previous columns 
gse$dex <- gse$condition

levels(gse$dex)
# See the levels, when renaming levels, the order MUST be preserved!!
levels(gse$dex) <- c("untrt", "trt")
# With this line of code we change the name of the dex levels
# Untreated -> untrt and Dexamethasone -> trt

# Check millions of fragments mapped for each sample (data presented in millions)
round(colSums(assay(gse)) / 1e6, 1)

# In R the first level of a factor should be the reference level (control or untreated samples)
# To do this we can use the [levels(...) <- ] function in the magrittr package

gse$dex %<>% relevel("untrt")
gse$dex

# The simplest design formula for differential expression would be [~ condition], 
# where condition is a column in colData(dds) that specifies which of two (or more groups) the samples belong to.
# For the airway experiment [~ cell + dex]
# This condition tests what the effect of dexamethasone (dex) controlling for the effect of different cell line (cell)

# Experimental design has to be specified at the begging of the analysis
dds <- DESeqDataSet(gse, design = ~ cell + dex) #dds stands for DESeqDataSet

# Exploratory analysis and visualization----------------------------------------------------------------------------

# To visually explore sample relationship use filtered data
# To perform statistical data  analysis use raw read counts

# First step on filtering data is removing any rows (genes) that do not have alot of counts 
# In this case we removed only genes that had 0 counts but it is possible to increase the number

nrow(dds)
keep <- rowSums(counts(dds)) > 1 #keep only rows whose sum is above 1 (removed 26680 rows)
dds <- dds[keep,]
nrow(dds)

# other type of filtering. Ex:to keep rows with at least 3 samples with a count of 10 or higher
# keep <- rowSums(counts(dds) >= 10) >= 3 (which would filter 41071 from the original dds)

# Transformation of counts

# For RNA-seq counts, the expected variance grows with the mean. 
# For example, if one performs PCA directly on a matrix of counts or normalized counts (e.g. correcting for differences in sequencing depth),
# the resulting plot typically depends mostly on the genes with highest counts because they show the largest absolute differences between samples.
# A simple and often used strategy to avoid this is to take the logarithm of the normalized count values plus a pseudocount of 1; 
# however, depending on the choice of pseudocount, now the genes with the very lowest counts will contribute a great deal of noise to the resulting plot, 
# because taking the logarithm of small counts actually inflates their variance.

# DESeq2 offers two transformations for count data that stabilize the variance across the mean: 
# the variance stabilizing transformation (VST) (runs faster)
# the regularized-logarithm transformation or rlog (use if the library size of samples vary widely)

vsd <- vst(dds, blind = FALSE) #blind = FALSE -> differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
rld <- rlog(dds, blind = FALSE)

head(assay(vsd), 3)
head(assay(rld), 3)
# Both vst and rlog return a DESeqTransform object which is based on the SummarizedExperiment class.
# The transformed values are no longer counts, and are stored in the assay slot.

# See the effect of transformation in the data (using only 2 samples)

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized = TRUE)[, 1:2] + 1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

df <- rename(df, x = c(1),
              y = c(2))

lvls <- c("log2(x + 1)", "vst", "rlog") #define levels
df$transformation <- factor(df$transformation, levels = lvls) #transformation becomes a factor with 3 levels instead of characters

ggplot(df, aes(x = "x", y = "y")) +
  geom_hex(bins = 80) +
  coord_fixed() +
  facet_grid(.~transformation)

# log2(x + 1) and rlog are similiar while VST has a upward shift for the smaller values.
# log2(x + 1) shows genes with low counts with a lot of variance (near the 0)
# vst and rlog however compress differences for low count genes 

# Test similiarty between samples (sample distance)
# Use the transformed data to ensure we have a roughly equal contribution from all genes

sampleDists <- dist(t(assay(vsd))) #t argument is used to transpose the matrix since the dist function expects the different samples to be rows and the different dimensions to be columns
sampleDists

# output the values into a heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$dex, vsd$cell, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_cols = sampleDists,
         cluster_rows = sampleDists,
         col = colors)

# Another way to test distance between samples is by performing a PCA
plotPCA(vsd, intgroup = c("dex", "cell")) #intgroup colored the point based on treatment and donor

# This was done using DESeq2, however it is also possible to do it usig ggplot2
pcaData <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

# The plot tells us that there is some distance between cells but treatment is much stronger variable
# Tutorial also has another type PCA to perform on not Normally distributed data

# MDS plot is another way to measure distance between samples
# MDS plot is very similar to PCA plots, and are made by using the multidimensional scaling (MDS) function
# Useful when we don’t have a matrix of data, but only a matrix of distances

mds <- as.data.frame(colData(vsd)) %>% 
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

# results were similar to the PCA

# Differential expression analysis ======================================================================================================================================
# experimental design was already specifeid when we did the DESeqDataSet
# differential expression pipeline on the raw counts with a single call to the function DESeq:

dds <- DESeq(dds)

# function calculates:
# the estimation of size factors (controlling for differences in the sequencing depth of the samples)
# the estimation of dispersion values for each gene
# fitting a generalized linear model.

# the following steps explain how to extract out results tables of interest from the new dds object

res <- results(dds)
res
# matrix could have been produced by the following command
# res <- results(dds, contrast = c("dex", "trt", "untrt"))

# results() extracts the estimated log2 fold changes and p values for the last variable in the design formula.
# If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level. meaning a comaprison between treatment for each cell cells 
# comparison is printed at the top of the output: dex trt vs untrt.

# res is a dataframe that carries information on the meaning of the columns
mcols(res, use.names = TRUE) #prints the meaning of the columns from the res dataframe

# baseMean: average of the normalized count values, divided by the size factors, taken over all samples in the DESeqDataSet.
# log2FoldChange: effect size estimate. How much the gene’s expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples
# lfcSE: standard error estimate for the log2 fold change estimate
# p-value: (avoid using) test null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability
# padj: (use instead of p-value) false discovery rate FDR (adjusted p value, uses a multiple comparison test) the one you should use

summary(res)
# summarises the results 
# FDR is below 10 % (padj < 0.1) and log2 fold change threshold is O (threshold that selects down and up regulated genes)
# we can be more strict about which set of genes are considered significant:

res.05 <- results(dds, alpha = 0.05) #by defining a lower FDR threshold (in other words probability of a type I error)
table(res.05$padj < 0.05) #5% FDR (filter genes with FDR lower than 0.05)

# and

resLFC1 <- results(dds, lfcThreshold = 1) # by increasing the log2 fold change threshold 
table(resLFC1$padj < 0.1) # lfcThreshold = 1 tests for gene counts that more than doubled or less that halved due to treatment with dex
                          # filter genes with FDR lower than 0.1

# compare res, resLFC1 and res.05 with summary()

# Results for a comparison of any two levels of a variable can be extracted using the contrast argument to results
results(dds, contrast = c("cell", "N061011", "N61311"))
# Specify three values: name of the variable (cell), level for the numerator (N061011) and denominaor (N61311)

# Multiple testing (prevents the increase of type 1 errors to occur, when conducting multiple statistical tests)
# In high-throughput biology, we are careful to not use the p values directly as evidence against the null, but to correct for multiple testing!!!
# DESeq2 uses the Benjamini-Hochberg (BH) adjustment

# Plotting results ======================================================================================================================================================

# Counts plot
# Visualize the counts for a specific gene with plotCounts
topGene <- rownames(res)[which.min(res$padj)] #print the gene with lowest padj
plotCounts(dds, gene = topGene, intgroup = c("dex"))

# This type of graph can also be done with ggplot2
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex", "cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() + geom_beeswarm(cex = 3)

ggplot(geneCounts, aes (x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line(linewidth = 1)

# through the graphs we can see that the gene being studied ha more counts in the treated samples (gene could be upregulated in the presence of the drug)

# MA-plot
# overview of the differences between measurements taken in two samples, by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values. Each point is a gene

# before doing the plot you need to shrink the log2 fold changes for the comparison of dex treated vs untreated samples
resultsNames(dds)
res <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")
plotMA(res, ylim = c(-5, 5))

# if shrinkage was not performed
res.noshr <- results(dds, name = "dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))

# blue points are significant data points

# It is possible to label individual points (genes) in the MA-plot
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topGene, pos = 2, col = "dodgerblue")
})

# Histogram of the p value
# Histogram is best formed by excluding genes with very small counts (othewise spikes appear in the histogram)

hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
# Histogram of p values for genes with mean normalized count larger than 1.

# Gene clustering
# similar to what was performed when we were evaluating distance between samples based on counts. in this case clustering is performed based on the highly variable genes
# we also work with the vst transformed data

topVarGene <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20) # the 20 genes with the highest variance across samples 

mat <- assay(vsd)[topVarGene, ]
mat <- mat - rowMeans(mat) # variable = the amount by which each gene deviates in a specific sample from the gene’s average across all samples (X - Mean)
anno <- as.data.frame(colData(vsd)[, c("cell", "dex")]) # provide a data.frame that instructs the pheatmap function how to label the columns.
pheatmap(mat, annotation_col = anno)

# Note that a set of genes in the heatmap are separating the N061011 cell line from the others, and there is another set of genes for which the dexamethasone treated samples have higher gene expression.

# Independent filtering

# MA-plot highlights a property of RNA-seq data.
# For weakly expressed genes, we have no chance of seeing differential expression. (due to low read counts)
# This can be observed by examining the ratio of small p values (say, less than 0.05) for genes binned by mean normalized count

qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6)) # create bins using the quantile function
bins <- cut(resLFC1$baseMean, qs) #bin the genes by base mean using cut
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2))) # rename the levels of the bins using the middle point
fractionSig <- tapply(resLFC1$pvalue, bins, function(p) # calculate the ratio of p values less than 0.05 for each bin
  mean(p < .05, na.rm = TRUE))

barplot(fractionSig, xlab = "mean normalized count",
                     ylab = "fraction of small p values")
# This plot demonstrates that genes with very low mean count have little or no power, and are best excluded from testing.

# By removing the low mean count genes we are improving the performance on the multiple testing adjustment.
# Therefore, we can find more genes to be significant among those that we keep, and so improved the power of our test. 
# This approach is known as independent filtering.
# DESeq2 automatically performs independent filtering

# Annotating and exploring results =====================================================================================================================================

# Annotating gene names to IDs
columns(org.Hs.eg.db) #list all available key types

ens.str <- substr(rownames(res), 1, 15) # extract gene ids
# had a column to res with gene Symbol
res$symbol <- mapIds(org.Hs.eg.db, 
                     keys = ens.str, # row names (geneIDs) provided as key
                     column = "SYMBOL", # which information you want to had 
                     keytype = "ENSEMBL", # specify keytype
                     multiVals = "first") # if there are multiple possible values for a single input value it uses the first value
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = ens.str,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)

# Exporting results
# save the results table in a csv file
# in this example we just export the 100 top genes

resOrderedDF <- as.data.frame(resOrdered)[1:100, ] # convert res into a data.frama and include only the top 100 genes
write.csv(resOrderedDF, file = "results_airway_tutorial.csv")

# another way of exporting data is using the ReportingTools package
# It generates dynamic HTML documents

htmlRep <- HTMLReport(shortName = "report_airway_tutorial", 
                      title = "My report",
                      reportDirectory = "./report_airway_tutorial")
publish(resOrderedDF, htmlRep)
finish(htmlRep)
browseURL("file:///C:/Users/tifig/OneDrive%20-%20Instituto%20Polit%C3%A9cnico%20do%20Porto/Ambiente%20de%20Trabalho/AcademicLife/R_project_training/report/report.html",
          browser = "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe")

# This would enough for the DESeq2 analysis, from here on out its just bonuses

# Plotting fold change in genomic space
# To do this we use the DESeqDataSet object and turn it into a GRange when shrinking the data
resGR <- lfcShrink(dds, coef = "dex_trt_vs_untrt", 
                   type = "apeglm", 
                   format = "GRanges") # output is GRanges (instead of a dataframe) which is used for the plot
resGR
ens.strGR <- substr(names(resGR), 1, 15)
resGR$symbol <- mapIds(org.Hs.eg.db,
                       keys = ens.strGR,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

window <- resGR[topGene] + 1e6 # window of 1 million base pairs upstream and downstream from the gene with the smallest p value
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window] # create a subset of our full results, for genes within the window
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol) # list of NA or duplicated genes
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)# add the gene symbol as a name if the symbol exists and is not duplicated in our subset

status <- factor(ifelse(resGRsub$padj < 0.05 & !is.na(resGRsub$padj), 
                        "sig","nostig")) # create a vector specifying if the genes in this subset had a low value of padj

options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack() # create an axis track specifying our location in the genome
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status) # create a track that will show the genes and their names
d <- DataTrack(resGRsub,
               data = "log2FoldChange",
               baseline = 0,
               type = "h",
               name = "log2 fold change",
               strand = "+") # create a data track that will draw vertical bars showing the moderated log fold change produced by DESeq2 (which are only large when the effect is well supported by the information in the counts)

# plot the results 
plotTracks(list(g, d, a), 
           groupAnnotation = "group",
           notsig = "grey",
           sig = "hotpink") #colored by significance

