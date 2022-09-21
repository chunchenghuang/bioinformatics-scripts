#install.packages("tidyverse")
#install.packages("svglite")
library(tidyverse)
library(ggplot2)
library(svglite)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

setwd("D:/heatmap/20201014")

#import .csv file
gene_list<-as.character(read.csv("D:/heatmap/20201014/GeneList.txt", header= TRUE)$Gene)
df<-read.table("D:/heatmap/20201014/Variant.txt", header = TRUE,sep="\t")

#dplyr
df_all<-df%>%
  group_by(Sample)%>%
  expand(Gene)%>%
  as.data.frame()
df2 <- merge(df, df_all, by=c("Sample","Gene"), all=TRUE)
df2$Consequence <- as.factor(df2$Consequence)

#check level order which would be directly changed to number by as.numeric 
Consequence<-as.factor(levels(df2$Consequence))
df2$Consequence <- as.numeric(df2$Consequence)

#fill in NA with 0s
df2$Consequence[is.na(df2$Consequence)] <- 0

#turn sampels and groups into factor
df2$Sample <- factor(df2$Sample)

#subset desired gene_list
gene_list<-c(gene_list)
df2_subset <- df2[df2$Gene %in% gene_list, ]
df2_subset$Gene <- factor(df2_subset$Gene, levels = c(rev(gene_list)))
#df2_subset<-df2

#calculate counts of mutated genes for bar chart
df2_count <- df2_subset %>%
  group_by(Gene) %>%
  summarise(Total_Samples= n(), Mutated_Samples = sum(Consequence != "0")) %>%
  arrange(desc(Mutated_Samples)) %>%
  mutate(Percentage=Mutated_Samples/Total_Samples) %>%
  as.data.frame()

#rearrange levels in Gene
df2_subset$Gene <- factor(df2_subset$Gene, levels = rev(as.vector(df2_count$Gene)))
df2_count$Gene <- factor(df2_count$Gene, levels = rev(as.vector(df2_count$Gene)))

#counts for bar chart
df2_count=df2_count[which(df2_count$Percentage!=0),]
df2_subset$Gene <- factor(df2_subset$Gene, levels = rev(as.vector(df2_count$Gene)))
df2_subset<-df2_subset[!is.na(df2_subset$Gene),]
df2_count$Gene <- factor(df2_count$Gene, levels = rev(as.vector(df2_count$Gene)))

#ggplot2
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

#Group_label = c( "1" ="Group 1", "2" ="Group 2")
gene_variants<- ggplot(df2_subset, aes(Sample, Gene)) +
  geom_tile(aes(fill=factor(df2_subset$Consequence)),colour="white") +
#facet_grid(.~Group, scales="free_x", space="free_x", labeller = labeller(Group=Group_label))+
  scale_x_discrete(expand=c(0,0))+ 
  scale_y_discrete(expand=c(0,0))+ 
  scale_fill_manual("Consequence",
                    values=c("lightgrey",getPalette(length(Consequence))),
                    breaks=c(as.numeric(Consequence)),
                    labels=c(as.character(Consequence)))+
  theme(legend.position = "bottom",
        legend.key.width=unit(3, "mm"),legend.key.height = unit(3, "mm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
        #axis.text.x = element_text(angle = 90, hjust = 1 , vjust = 0.5))+
  guides(fill=guide_legend(nrow=2))

df2_subset_bar<-df2_subset[which(df2_subset$Consequence != "0"),]
variant_counts <- ggplot(df2_subset_bar, aes(factor(Gene)))+
  geom_bar(aes(fill=factor(df2_subset_bar$Consequence)),colour="white")+
  scale_x_discrete(labels=rev(c(percent(df2_count$Percentage, accuracy = 1))),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  #scale_y_continuous(position = "right")+
  scale_fill_manual(values=c("lightgrey",getPalette(length(Consequence))))+
  coord_flip()+
  theme(legend.position = "none",
        aspect.ratio=2.5,
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

final <- ggarrange(gene_variants, variant_counts, ncol = 2, common.legend=TRUE, legend = "bottom", widths=c(2, 1), heights=c(2, 1))

ggsave(file="D:/heatmap/20201014/heatmap_20201014.tiff", plot=final, device="tiff", width=10, height=8)
ggsave(file="D:/heatmap/20201014/heatmap_20201014.svg", plot=final, device="svg", width=10, height=8)
