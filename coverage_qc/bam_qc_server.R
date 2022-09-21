#install
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Rsamtools")
#BiocManager::install("rtracklayer")
#BiocManager::install("seqTools")

#library
args <- commandArgs(trailingOnly = TRUE)

dir_path<-args[1]
stat_file<-args[2]
stat_region_file<-args[3]
qual_table<-args[4]
DP_file<-args[5]
depth_lower_than<-args[6]
depth_lower_than<-as.numeric(depth_lower_than)
bam_file<-args[7]
setwd(dir_path)

stat <- read.table(stat_file, sep = '\t', row.names = 1)
stat_region <- read.table(stat_region_file, sep = '\t', row.names = 1)

#Total reads
TR <- stat["reads paired:",]
TRR <- stat_region["reads paired:",]

#average insert size
MIS <-stat["insert size average:",]

#Total mapped reads
MR <- stat["reads mapped:",]

#Total mapped reads in target regions
MRR <- stat_region["reads mapped:",]

#Depth
df <- read.table(DP_file, sep = '\t')

#Target region length
BR <- stat_region["bases inside the target:",]

#summing quality

q_table <- read.table(qual_table, sep = '\t')
sum_base <- sum(q_table)
bq_20 <- sum(tail(q_table, n=nrow(q_table)-21))
bq_30 <- sum(tail(q_table, n=nrow(q_table)-31))
qual_20 <- (bq_20/sum_base)*100
qual_30 <- (bq_30/sum_base)*100


#Mean Depth 
MD <- mean(df$V3)

#Percentage of bases covers >= 20% of mean
df_0.2 <- df[which(df$V3>=MD*0.2),]
P20M <- (nrow(df_0.2)/BR)*100

#Percentage of bases covers >= 10x
df_10 <- df[which(df$V3>=10),]
P10 <- (nrow(df_10)/BR)*100

#Percentage of bases covers >= 30x
df_30 <- df[which(df$V3>=30),]
P30 <- (nrow(df_30)/BR)*100

#Percentage of bases covers >= 100x
df_100 <- df[which(df$V3>=100),]
P100 <- (nrow(df_100)/BR)*100

#Percentage of bases covers >= 200x
df_200 <- df[which(df$V3>=200),]
P200 <- (nrow(df_200)/BR)*100

#Percentage of bases covers >= 250x
df_250 <- df[which(df$V3>=250),]
P250 <- (nrow(df_250)/BR)*100

##Percentage of duplication
DP <- (stat_region["reads duplicated:",]/TRR)*100

##Output table
Criteria <- c("Sample ID", "Total reads", "Median insert size", "% of bases >= Q20 in target region", 
"% of bases >= Q30 in target region", "Total mapped reads", "Total mapped reads in target region",
"Target region length", "Mean read depth in target region",
"% of bases covered at >= mean*0.2", "% of bases covered at >= 10x",
"% of bases covered at >= 30x", "% of bases covered at >= 100x",
"% of bases covered at >= 200x", "% of bases covered at >= 250x",
"% of duplication")

file_name<-basename(bam_file)
QC <- c(file_name,TR,MIS,qual_20,qual_30,MR,MRR,BR,MD,P20M,P10,P30,P100,P200,P250,DP)

output<-data.frame(cbind(Criteria,QC))
rownames(output) <- output[,1]
output[,1]<-NULL
QC_csv<-paste(file_name,"_QC.csv",sep="")
write.csv(output,file=QC_csv)

#output the desired low depth region
df_low <- df[which(df$V3<depth_lower_than),]
w<-paste0("No region's depth is lower than ", depth_lower_than)
if(nrow(df_low)==0){warning(w);quit(save="no",status=0)}
colnames(df_low) <- c("chr","region","depth")
end_position<-c(which(diff(df_low[,2]) != 1),length(df_low[,2]))
start_position<-c(1,which(diff(df_low[,2]) != 1)+1)
region_pairs<-length(end_position)
mean<-c()
for(i in 1:region_pairs){
  region_depth<-df_low[start_position[i]:end_position[i],3]
  region_mean<-mean(region_depth)
  mean[i]<-region_mean
}
chr<-df_low[end_position,1]
low_coverage_regions<-data.frame(chr=chr,start=c(df_low[start_position,2]-1),end=c(df_low[end_position,2]),mean=mean)
sample_name<-head(unlist(strsplit(file_name, "_")),1)
#output_txt<-paste(sample_name,"_low_coverage_regions.txt",sep="")
output_csv<-paste(sample_name,"_low_coverage_regions.csv",sep="")
#write.table(low_coverage_regions,file=output_txt,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
write.csv(low_coverage_regions,file=output_csv,row.names=FALSE)
