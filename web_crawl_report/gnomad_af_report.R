library(httr)
library(jsonlite)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

dir<-args[1]
cdir<-args[2]

#resource preparation
txt_ls<-list.files(path=dir,pattern=".txt",full.names = TRUE)
tbbhg19af_file<-list.files(path=cdir,pattern="TaiwanBioBank_hg19.csv",full.names = TRUE)
w<-"No vcf or vcf.gz files detected!"
if(length(txt_ls)==0){warning(w);quit(save="no",status=0)}
tbbhg19af<-fread(tbbhg19af_file) #fread reads large files faster
wt<-"No taiwan biobank allele frequency file detected!"
if(length(tbbhg19af_file)==0){warning(wt);quit(save="no",status=0)}

#loop through each vcf file
for(f in c(1:length(txt_ls))){
  vcf_file<-txt_ls[f]
  vcf <- read.table(vcf_file, header=T,  sep='\t', fill=T, comment.char = "#")
  #loop through each variant in each vcf file
  gaf<-c()
  eaaf<-c()
  taf<-c()
  for(i in c(1:nrow(vcf))){
    cacc<-vcf[i,2]
    pos<-vcf[i,3]
    ref<-as.character(vcf[i,4])
    alt<-as.character(vcf[i,5])
    chromv<-paste("\"",cacc,"-",pos,"-",ref,"-",alt,"\"",sep="")
    search<-paste("https://gnomad.broadinstitute.org/api?query={variant(dataset:gnomad_r2_1%20variantId:",chromv,"){genome{ac%20an%20populations{id%20ac%20an}}}}",sep="")
    data<-fromJSON(search)
    if (is.null(data$data$variant$genome)){
      gaf[i]<-NA
      eaaf[i]<-NA}else{
        gaf[i]<-data$data$variant$genome$ac/data$data$variant$genome$an
        popdf<-data.frame(data$data$variant$genome$populations)
        eaaf[i]<-popdf[which(popdf$id=="EAS"),]$ac/popdf[which(popdf$id=="EAS"),]$an
      }
    tbk<-as.numeric(tbbhg19af[which(tbbhg19af$Chr==cacc & tbbhg19af$Coordinate==pos& tbbhg19af$Alt==alt),4])
    if(is.na(tbk)){
      taf[i]<-NA}else{
        taf[i]<-tbk
      }
  }
  vcf["Global_AF"]<-gaf
  vcf["EAS_AF"]<-eaaf
  vcf["TAIWAN_AF"]<-taf
  ##merge with header and write output
  vcf_output<-gsub(vcf_file, pattern=".txt$", replacement=".genomad.txt")
  write.table(vcf, vcf_output, quote = FALSE, row.names = FALSE, sep="\t")
}