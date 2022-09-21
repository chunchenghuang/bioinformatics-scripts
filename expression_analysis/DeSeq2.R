##Install the library packages
install.packages("htmltools")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("apeglm")
BiocManager::install("countToFPKM")
BiocManager::install("biomaRt")
install.packages("pheatmap")
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("apeglm")

##library packages
library(htmltools)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ggpubr)
library(pasilla)
library(countToFPKM)
library(biomaRt)
library(dplyr)

##http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

##import data

pheno_data = read.csv("D:/rna_seq/SWcellline_20191202/SWcelline_20191202_list.csv", head=TRUE)
mx = read.table("D:/rna_seq_denovo_IDS/MPS_genenames.mx", head=TRUE)
row.names(mx) <- as.character(unlist(mx[,1]))

##as.vector(zebrafinch_brain_pre_3_mx[,1]) also works, unlist makes it to vectors
##make your sameples and gene names into rownames for pheno_data and your matrix
mx<-mx[,-1]
row.names(pheno_data)<-pheno_data[,1]
pheno_data[1]<-NULL
#or
#pheno_data[,-1 ,drop=FALSE]

all(rownames(pheno_data) %in% colnames(mx))
##%in% checks whether or not the object is contained in the other object. Whereas == is a logical operator that checks for identity properties.

all(rownames(pheno_data) == colnames(mx))
##should be true, the sample names need to be in the same order and number

dds <-DESeqDataSetFromMatrix(mx, pheno_data, ~ Treatment)

##Differential Expression Analysis
dds<-DESeq(dds)
head(dds)

##ouput summary
res <- results(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Treatment_Treated_vs_Control", type="apeglm")
plotMA(res, alpha=0.05, ylim=c(-2,2))
plotMA(resLFC, alpha=0.05, ylim=c(-2,2))

##subset pvalue < 0.05
res0.05<-subset(res, res$pvalue<0.05)

##order the table from small pvale to large pvale
res0.05 <- res0.05[order(res0.05$pvalue), ]
write.csv(res0.05, "sig0.05_order.csv")

##pvalue<0.01 genes
res0.01<-subset(res, res$pvalue<0.01)
res0.01<-res0.01[order(res0.01$pvalue),]
write.csv(res0.01, "sig_results_0.01.csv")

##heatmap
# this gives log2(n + 1)
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
df[2]<-NULL
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

