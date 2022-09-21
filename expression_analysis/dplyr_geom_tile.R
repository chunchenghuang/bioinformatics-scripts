#install.packages("tidyverse")
#install.packages("svglite")
#install.packages("gridExtra")
#install.packages("ggpubr")
library(tidyverse)
library(ggplot2)
library(svglite)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

#import .csv file
gene_list<-as.character(read.table("D:/gene_mutation_heatmap/gene_list.txt",sep="\t", header= FALSE)$V1)
df<-read.table("D:/gene_mutation_heatmap/Pan_cancer_variaint_list.txt",sep="\t", header = TRUE)

#dplyr
df_all<-df%>%
  group_by(Group, Sample)%>%
  expand(Gene)%>%
  as.data.frame()
df2 <- merge(df, df_all, by=c("Group","Sample","Gene"), all=TRUE)
df2$Consequence <- as.factor(df2$Consequence)
df2<-unique(df2)

#check level order which would be directly changed to number by as.numeric 
Consequence<-as.factor(levels(df2$Consequence))
df2$Consequence <- as.numeric(df2$Consequence)

#fill in NA with 0s
df2$Consequence[is.na(df2$Consequence)] <- 0

#turn sampels and groups into factor
df2$Sample <- factor(df2$Sample)
df2$Group <- factor(df2$Group)

#subset desired gene_list
gene_list<-c(gene_list)
df2_subset<- df2[df2$Gene %in% gene_list, ]


#remove duplicated variants, only one variant per sample per gene, leaving the most important one
df2_subset_max <- df2_subset[which(df2_subset$Consequence==max(df2_subset$Consequence)),]
df2_subset_merge <- merge(df2_subset, df2_subset_max, by=c("Group","Sample","Gene"), all=TRUE)
df2_subset_merge <- df2_subset_merge[-c(which(df2_subset_merge$Consequence.y==2)),]
df2_subset_merge$Consequence.y <- NULL
df2_subset <- df2_subset_merge %>%
  rename(Consequence = Consequence.x) %>%
  bind_rows(df2_subset_max) %>%
  arrange(Consequence) %>%
  as.data.frame()

#separate groups
df2_subset_1<-df2_subset[which(df2_subset$Group==1),]
df2_subset_2<-df2_subset[which(df2_subset$Group==2),]

#calculate counts of mutated genes in group 1
df2_count_1 <- df2_subset_1 %>%
  group_by(Gene) %>%
  summarise(Total_Samples= n(), Mutated_Samples = sum(Consequence != "0")) %>%
  arrange(desc(Mutated_Samples)) %>%
  mutate(Percentage=Mutated_Samples/Total_Samples) %>%
  as.data.frame()

#rearrange levels in Gene
df2_subset_1$Gene <- factor(df2_subset_1$Gene, levels = rev(as.vector(df2_count_1$Gene)))
df2_count_1$Gene <- factor(df2_count_1$Gene, levels = rev(as.vector(df2_count_1$Gene)))

#calculate counts of mutated genes in group 1
df2_count_2 <- df2_subset_2 %>%
  group_by(Gene) %>%
  summarise(Total_Samples= n(), Mutated_Samples = sum(Consequence != "0")) %>%
  arrange(desc(Mutated_Samples)) %>%
  mutate(Percentage=Mutated_Samples/Total_Samples) %>%
  as.data.frame()

#rearrange levels in Group 1
df2_count_1=df2_count_1[which(df2_count_1$Percentage!=0),]
df2_subset_1$Gene <- factor(df2_subset_1$Gene, levels = rev(as.vector(df2_count_1$Gene)))
df2_subset_1<-df2_subset_1[!is.na(df2_subset_1$Gene),]
df2_count_1$Gene <- factor(df2_count_1$Gene, levels = rev(as.vector(df2_count_1$Gene)))

#rearrange levels in Group 2 
df2_count_2=df2_count_2[which(df2_count_2$Percentage!=0),]
df2_subset_2$Gene <- factor(df2_subset_2$Gene, levels = rev(as.vector(df2_count_2$Gene)))
df2_subset_2<-df2_subset_2[!is.na(df2_subset_2$Gene),]
df2_count_2$Gene <- factor(df2_count_2$Gene, levels = rev(as.vector(df2_count_2$Gene)))

#subset consequence levels

#ggplot2

#colour
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

#Group label
#Group_label = c( "1" ="Group 1", "2" ="Group 2")

#heatmap
gene_variants_1<- ggplot(df2_subset_1, aes(Sample, Gene)) + 
  geom_tile(aes(fill=factor(df2_subset_1$Consequence)),colour="white") +
  #facet_grid(.~Group, scales="free_x", space="free_x", labeller = labeller(Group=Group_label))+
  scale_x_discrete(expand=c(0,0))+ 
  scale_y_discrete(expand=c(0,0))+ 
  scale_fill_manual("Consequence",
                    values=c("lightgrey",getPalette(length(Consequence))),
                    breaks=c(as.numeric(Consequence)),
                    labels=c(as.character(Consequence)))+
  theme(legend.position = "bottom",
        legend.key.width=unit(3, "mm"),legend.key.height = unit(3, "mm"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank())+
        axis.text.x = element_text(angle = 90, hjust = 1 , vjust = 0.5))+
  guides(fill=guide_legend(nrow=2))



gene_variants_2<- ggplot(df2_subset_2, aes(Sample, Gene)) + 
  geom_tile(aes(fill=factor(df2_subset_2$Consequence)),colour="white") +
  #facet_grid(.~Group, scales="free_x", space="free_x", labeller = labeller(Group=Group_label))+
  scale_x_discrete(expand=c(0,0))+ 
  scale_y_discrete(expand=c(0,0))+ 
  scale_fill_manual("Consequence",
                    values=c("lightgrey",getPalette(length(Consequence))),
                    breaks=c(as.numeric(Consequence)),
                    labels=c(as.character(Consequence)))+
  theme(legend.position = "bottom",
        legend.key.width=unit(3, "mm"),legend.key.height = unit(3, "mm"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank())+
        #axis.text.x = element_text(angle = 90, hjust = 1 , vjust = 0.5),
        axis.title.x=element_blank())+
  guides(fill=guide_legend(nrow=2))

#bar chart
df2_subset_bar_1<-df2_subset_1[which(df2_subset_1$Consequence != "0"),]
variant_counts_1 <- ggplot(df2_subset_bar_1, aes(factor(Gene)))+
  geom_bar(aes(fill=factor(df2_subset_bar_1$Consequence)),colour="white")+
  scale_x_discrete(labels=rev(c(percent(df2_count_1$Percentage, accuracy = 1))),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  #scale_y_continuous(position = "right")+
  scale_fill_manual(values=c(getPalette(length(Consequence))))+
  coord_flip()+
  theme(legend.position = "none",
        aspect.ratio=3.1,
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df2_subset_bar_2 <- df2_subset_2[which(df2_subset_2$Consequence != "0"),]
variant_counts_2 <- ggplot(df2_subset_bar_2, aes(factor(Gene)))+
  geom_bar(aes(fill=factor(df2_subset_bar_2$Consequence)),colour="white")+
  scale_x_discrete(labels=rev(c(percent(df2_count_2$Percentage, accuracy = 1))),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  #scale_y_continuous(position = "right")+
  scale_fill_manual(values=c(getPalette(length(Consequence))))+
  coord_flip()+
  theme(legend.position = "none",
        aspect.ratio=3.1,
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#arrange plots
final_1 <- ggarrange(gene_variants_1, variant_counts_1,ncol = 2, common.legend=TRUE, legend = "bottom", widths=c(2, 1), heights=c(2, 1))
final_2 <- ggarrange(gene_variants_2, variant_counts_2,ncol = 2, common.legend=TRUE, legend = "bottom", widths=c(2, 1), heights=c(2, 1))

final_2
final_1

ggsave(file="final_1.svg", plot=final_1, device="svg", width=10, height=8)
ggsave(file="final_2.svg", plot=final_2, device="svg", width=10, height=8)

