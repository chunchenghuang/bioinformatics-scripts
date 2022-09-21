#install.packages("dendextend")
#install.packages("pheatmap")
#install.packages("tidyverse")

library(dplyr)
library(dendextend)
library(pheatmap)

setwd("D:/heatmap/20200324/")

df <- read.csv("D:/heatmap/20200324/TPM_color.csv", header = TRUE, check.names=FALSE, row.names = 1)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

df_norm <- t(apply(df, 1, cal_z_score))

my_hclust_gene <- hclust(dist(df_norm), method = "complete")

##I found out my original data has too low expression level in general (maybe because I have half of the data? or it��s just brain��s expression level is generally low?), therefore was unable to do gene clustering. So, I used z-score normalized data instead of original data

pdf("heatmap_dendrogram.pdf")

as.dendrogram(my_hclust_gene) %>%
  
  plot(horiz = TRUE)

dev.off()

##%>% = | in linux

##my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 8)

##used cutree_rows so did not use gene cluster anymore

##columns by individual in brain region's order

##reorder the column so we align desired brain regions together

my_sample_col <- as.data.frame(pheno_data)
row.names(my_sample_col) <- my_sample_col[,1]
my_sample_col <- my_sample_col[,-1]
my_sample_col[,3] <- data.frame(Individual = c(1,2,3,3,1,1,1,2,2,2,3,3,4,4,4,4))
colnames(my_sample_col)[1]<-"Brain Region"

pdf("heatmap.pdf")

pheatmap(df_norm, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = TRUE)

dev.off()
