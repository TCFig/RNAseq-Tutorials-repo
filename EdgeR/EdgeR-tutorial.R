# The tutorial was taken from https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html 
# tag = gene

library(edgeR)
library(baySeq) # Example data set used for data analysis

# Putting the data into the right format for edgeR =====================================================================================================================

# Data set is a matrix (mobData) of counts acquired for three thousand small RNA loci from a set of Arabidopsis grafting experiments

# First extract data from the package
dir <- system.file("data", package = "baySeq", mustWork = T)
load(file.path(dir, "mobData.RData"))
head(mobData)
mobData[1,] <- c(21, 52, 4, 4, 86, 68) #the first row was wrong so i had to change the values


help("mobData") # see more info about the dataframe

# create groups in the data
mobDataGroups <- c("MM", "MM", "WM", "WM", "WW", "WW")
# MM = triple mutant shoot grafted onto triple mutant root
# WM = wild-type shoot grafted onto triple mutant root
# WW = wild-type shoot grafted onto wild-type root
data("mobAnnotation")
?mobAnnotation
head(mobAnnotation)

# edgeR works on a table of integer read counts, with rows corresponding to genes and columns to independent libraries
# edgeR stores data in an object called DGEList, which can be manipulated like any list in R
# You can make this in R by specifying the counts and the groups in the function DGEList()

d <- DGEList(counts = mobData, group = factor(mobDataGroups))
d

# All information is now contained in the object d
# other functions like d <- estimateCommonDisp(d) will add one more element to the list instead of overwriting d

# Filtering the data =====================================================================================================================================================

# First get rid of genes which did not occur frequently enough
# Example, remove genes with less than 100 counts per million (calculated with cpm()) for at least 2 samples

d.full <- d # keep the old DGELIst in case something goes wrong

apply(d$counts, 2, sum) # full library size

keep <- rowSums(cpm(d) > 100) >= 2 # keep genes with at least 100 counts per million for at least 2 samples
d <- d[keep, ]

dim(d)
dim(d.full) 

# 3000 genes reduced to 787. 
# For the filtered tags, there is very little power to detect differential expression, so little information is lost by filtering.

# After filtering, it is a good idea to reset the library sizes
d$samples$lib.size <- colSums(d$counts)
d.full$samples #library sizes reduce after filtering
d$samples

# size factor from DESeq is different from "norm factor" in edgeR.
# In edgeR library size and additional normalization are two separated columns.
# In all the downstream code, the lib.size and norm.factors are multiplied together to act as the effective library size

# Normalizing the data =====================================================================================================================================================

# edgeR normalizes by total count
# edgeR is concerned with differential expression analysis rather than with the qauntification of expression levels
# It is concerned with relative changes in expression levels between conditions, but not directly with estimating absolute expression levels

# calcNormFactors() normalizes for RNA composition by finding a set of scaling factors for the library sizes taht minimize the log-fold changes between the samples for most genes
# The default method computes for TMM (trimmed mean of M-values between each pair of samples)

d <- calcNormFactors(d)
d$samples # changed norm.factors from 1 to a different value

# Data Exploration =====================================================================================================================================================

# Produce a plot to see sample relations based on multidimensional scaling
# MDS-plot

plotMDS(d, method = "logFC", col = as.numeric(d$samples$group)) # makes the plot
legend("bottomleft", as.character(unique(d$samples$group)), col = 1:3, pch = 20) # add a legend to the bottomleft with the group of each sample

# Estimating the Dispersion =====================================================================================================================================================

# In the negative binomial model for RNA-seq data, each gene is given a dispersion parameter, and correctly estimating these dispersion parameters is vital to detecting differential expression.

# To do a DGE analysis using the NB (quasi negative binomial) model first it is required to estimate the dispersion parameter for each tag 
# Dispersion = measure of degree of inter-library variation for that tag

# Estimate dispersion by assuming everything has the same common dispersion
# We could also use a generalized linear model to estimate dispersion

# Naive form of estimating dispersion
d1 <- estimateCommonDisp(d, verbose = T) # first method was used to estimate dispersion
names(d1)

# For differential gene expression analysis it is common to use the Bayes tagwise dispersions
# Note: common dispersion needs to be estimated before estimating tagwise

d1 <- estimateTagwiseDisp(d1)
names(d1)

plotBCV(d1) # plots the tagwise biological coefficient of variation against log2-CPM

# Looking at the graph we see that a single estimate for the coefficient of variation is a bad model 
# Since tagwise dispersion does not follow the model but instead increases as the cpm increases

# GLM estimates of dispersion

# Fitting a model in edgeR
# 1st: fit the common dispersion
# 2nd: fit a trended model (if you do not fit a trend, the default is to use the common dispersion as a trend)
# 3rd: fit the tagwise dispersion which is a function of this model

design.mat <- model.matrix(~ 0 + d$samples$group) # creates dummy variables
colnames(design.mat) <- levels(d$samples$group) # change the name of the columns to match the groups
d2 <- estimateGLMCommonDisp(d, design = design.mat)
d2 <- estimateGLMTrendedDisp(d2, design = design.mat, method = "power")
d2 <- estimateGLMTagwiseDisp(d2, design = design.mat)
plotBCV(d2)

# In this method we model the tagwise dispersion based on the model derived from the glm model choosen

# Comparing the models in DESeq and edgeR

# DESeq only uses the gamma glm as its model (edgeR does not have gamma glm model)

# Dispersion and Biological Variation

# The dispersion of a gene is a measure of a gene's variance
# The dispersion can be interpreted as the square of the coefficient of biological variation
# Ex: difference in counts between two biological replicates is 40% so the genes disperion is 0.4^2 = 0.16

# Differential Expression =====================================================================================================================================================

# the function exactTest() conducts tagwise tests using the exact negative binomial test
# topTags() shows the n most significant tags
# BH algorithm is by default used to control the false discovery rate (FDR)

# Differential expression when using the naive data d1
et12 <- exactTest(d1, pair = c(1,2)) #Compares groups 1 and 2
et13 <- exactTest(d1, pair = c(1,3)) #Compares groups 1 and 3
et23 <- exactTest(d1, pair = c(2,3)) #Compares groups 2 and 3

topTags(et12, n = 10)
topTags(et13, n = 10)
topTags(et23, n = 10)

# function decideTestsDGE() identifies which genes are differently expressed

de1 <- decideTestsDGE(et12, adjust.method = "BH", p.value = 0.05, lfc = 2) # log fold change threshold is by defaul [-1, 1]
summary(de1)

# function plotSmear() generates a plot of the tagwise log-fold-change against log-cpm
# DE tags are highlighted on the plot

de1tags12 <- rownames(d1)[as.logical(de1)]
plotSmear(et12, de.tags = de1tags12)
abline(h = c(-2, 2), col = "blue")


# Differential expression using the GLM data

fit <- glmFit(d2, design = design.mat) # necessary to proceed with genewise statistical tests


lrt12 <- glmLRT(fit, contrast = c(1,-1,0)) # 1: group 1, -1: compared with group 2, 0: group 3 is out
lrt13 <- glmLRT(fit, contrast = c(1,0,-1))
lrt23 <- glmLRT(fit, contrast = c(0,1,-1))

topTags(lrt12, n = 10)
topTags(lrt13, n = 10)
topTags(lrt23, n = 10)

de2 <- decideTestsDGE(lrt12, adjust.method = "BH", p.value = 0.05, lfc = 2)
de2tags12 <- row.names(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags = de2tags12)
abline(h = c(-2,2), col = "blue")
