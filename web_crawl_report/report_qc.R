library(XML)
library(httr)

args <- commandArgs(trailingOnly = TRUE)

dir<-args[1]
#dir<-"D:/TIP_FANGS_程式撰寫/QIAseq-DNA_132452_CDHS-32578Z-210"
cdir<-args[2]
#cdir<-"D:/TIP_FANGS_程式撰寫"

#final output file name
date<-paste(strsplit(as.character(Sys.Date()),"-")[[1]],collapse="")
final_output_csv<-paste(strsplit(basename(dir),"_")[[1]][2],date,"report.csv",sep="_")
final_output_ls<-list.files(path=dir,pattern=final_output_csv,full.names = TRUE)
final_output_csv<-paste(dir,final_output_csv,sep="/")
w<-"Final report already existed!"
if(length(final_output_ls)!=0){warning(w);quit(save="no",status=0)}


#resource preparation
txt_ls<-list.files(path=dir,pattern=".anno.txt",full.names = TRUE)
w<-"No vcf text files detected!"
if(length(txt_ls)==0){warning(w);quit(save="no",status=0)}

Disease.levels <- list(
  "CPT1" = c('CPT1A'),
  "CPT2" = c('SLC25A20', 'CPT2'),
  "GA2" = c('ETFA','ETFB','ETFDH'),
  "VLCAD" =c('ACADVL'))

gene.list <- c('CPT1A','SLC25A20', 'CPT2','ETFA','ETFB','ETFDH','ACADVL')

chromaccs.levels<-list(
  "NC_000001.10" = "chr1",
  "NC_000002.11" = "chr2",
  "NC_000003.11" = "chr3",
  "NC_000004.11" = "chr4",
  "NC_000005.9" = "chr5",
  "NC_000006.11" = "chr6",
  "NC_000007.13" = "chr7",
  "NC_000008.10" = "chr8",
  "NC_000009.11" = "chr9",
  "NC_000010.10" = "chr10",
  "NC_000011.9" = "chr11",
  "NC_000012.11" = "chr12",
  "NC_000013.10" = "chr13",
  "NC_000014.8" = "chr14",
  "NC_000015.9" = "chr15",
  "NC_000016.9" = "chr16",
  "NC_000017.10" = "chr17",
  "NC_000018.9" = "chr18",
  "NC_000019.9" = "chr19",
  "NC_000020.10" = "chr20",
  "NC_000021.8" = "chr21",
  "NC_000022.10" = "chr22",
  "NC_000023.10" = "chrX",
  "NC_000024.9" = "chrY")

UA <- "Mozilla/5.0 (Linux; Android 6.0.1; Nexus 5X Build/MMB29P) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2272.96 Mobile Safari/537.36 (compatible; Googlebot/2.1; +http://www.google.com/bot.html)"

nm_transcripts<-c()
np_transcripts<-c()
for(g in c(1:length(gene.list))){
  gene<-gene.list[g]
  nuccore_search<-paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=",gene,"%5BGene%20Name%5D+AND+Homo+sapiens%5BOrganism%5D+AND+Refseq_select%5Bfilter%5D",sep="")
  nuccore_doc <- GET(nuccore_search, user_agent(UA))
  nuccore_data <- xmlParse(content(nuccore_doc, "text"))
  id=xmlValue(xmlRoot(nuccore_data)[["IdList"]][["Id"]])
  Sys.sleep(2.5)
  transcript_search= paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=xml&is_variationid&id=",id,sep="")
  transcript_doc <- GET(transcript_search, user_agent(UA))
  transcript_data <- xmlParse(content(transcript_doc, "text"))
  nm_transcripts[g] <- xmlValue(getNodeSet(xmlRoot(transcript_data), "//GBSeq_accession-version"))
  npta <- xmlValue(getNodeSet(xmlRoot(transcript_data), "//GBQualifier_value"))
  np_transcripts[g] <-npta[grep("^NP", npta)]
  Sys.sleep(2.5)
}

dp_ls<-list.files(path=dir,pattern=".variant-calling-input.bedgraph",full.names = TRUE)

commonsnp_file<-list.files(path=cdir,pattern="commonSNPhg19EASandTaiwan.csv",full.names = TRUE)
w<-"Missing commonsnp file!"
if(length(commonsnp_file)==0){warning(w);quit(save="no",status=0)}

commonsnphg19 <- read.csv(commonsnp_file, header=F, fill=T)

output_list<-list()
sample_name<-c()
qc=c()
for(f in c(1:length(txt_ls))){
  ##Clinvar annotation
  vcf_file<-txt_ls[f]
  sample_name[f]<-strsplit(basename(txt_ls[f]),"[_.-]")[[1]][3]
  #read in vcf to parse
  vcf<-read.table(vcf_file,sep="\t",header=TRUE,stringsAsFactors = FALSE)
  #filter snp
  vcf <- vcf[-c(which(paste(vcf$CHROM,vcf$POS,vcf$ALT,sep='-') %in% commonsnphg19$V1)),]
  #filter VCF with variant's frequency > 0.25
  if(nrow(vcf)!=0){vcf <- vcf[which(vcf["VMF"]>=0.25),]}
  if(nrow(vcf)==0){
    output<-vcf[,c("Gene_Name","VMF")]
    output[1,]<-rep("-",ncol(output))
    output[,c("NCBI Transcript ID","ClinVar Accession","Clinical Significance","cDNA_change","AA_change","Disease")]<-"-"
    output["SampleID"]=sample_name[f]}else{
      vcf$chromaccs <- as.factor(vcf$CHROM)
      levels(vcf$chromaccs) <- chromaccs.levels
      rcvs=c()
      dess=c()
      NM_num=c()
      cchange=c()
      pchange=c()
      for(i in c(1:nrow(vcf))){
        cacc<-vcf[i,"chromaccs"]
        pos<-vcf[i,"POS"]
        alt<-vcf[i,"ALT"]
        ref<-vcf[i,"REF"]
        gname<-vcf[i,"Gene_Name"]
        chromv<-paste(cacc,"%3Ag.",pos,ref,">",alt,"%20",gname,sep="")
        id_search <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=",chromv,sep="")
        id_doc <- GET(id_search, user_agent(UA))
        id_data <- xmlParse(content(id_doc, "text"))
        id=xmlValue(xmlRoot(id_data)[["IdList"]][["Id"]])
        if(length(xmlToList(xmlRoot(id_data)[["IdList"]]))!=1){
          #gspdi<-"NC_000004.11:159616699:G:A"
          gspdi<-paste(cacc,":",pos-1,":",ref,":",alt,sep="")
          vrecorder<-paste("https://grch37.rest.ensembl.org/variant_recoder/human/",gspdi,"?content-type=text/xml",sep="")
          nm_doc <- GET(vrecorder, user_agent(UA))
          nm_data <- xmlParse(content(nm_doc, "text"))
          nm_hgvsc <- xmlValue(getNodeSet(xmlRoot(nm_data), "//hgvsc"))
          nm_subset<- grep("^NM", nm_hgvsc)
          nm_number_ls<-strsplit(nm_hgvsc[nm_subset],"[:]")
          np_hgvsp <- xmlValue(getNodeSet(xmlRoot(nm_data), "//hgvsp"))
          np_subset<- grep("^NP", np_hgvsp)
          np_number_ls<-strsplit(np_hgvsp[np_subset],"[:]")
          if (length(nm_number_ls)==0){
            NM_num[i]<-"-"
            cchange[i]<-"-"
          }else{
            for (y in 1:length(nm_number_ls)){
              if (nm_number_ls[[y]][1] %in% nm_transcripts == TRUE){
                NM_num[i]<-nm_number_ls[[y]][1]
                cchange[i]<-nm_number_ls[[y]][2]
                break
              }else{
                NM_num[i]<-"-"
                cchange[i]<-"-"
              }
            }
          }
          if (length(np_number_ls)==0){
            pchange[i]<-"-"
          }else{
            for (y in 1:length(np_number_ls)){
              if (np_number_ls[[y]][1] %in% np_transcripts == TRUE){
              pchange[i]<-np_number_ls[[y]][2]
              break
              }else{
                pchange[i]<-"-"
              }
            }
          }
          chromv<-paste(NM_num[i],"%3A",cchange[i],sep="")
          Sys.sleep(2.5)
          id_search <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=",chromv,sep="")
          id_doc <- GET(id_search, user_agent(UA))
          id_data <- xmlParse(content(id_doc, "text"))
          id=xmlValue(xmlRoot(id_data)[["IdList"]][["Id"]])
          if(length(xmlToList(xmlRoot(id_data)[["IdList"]]))==0){
            rcvs[i]<-"-"
            dess[i]<-"Uncertain Significance"
            next
          }
          }
        Sys.sleep(2.5)
        rcv_search= paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=clinvarset&is_variationid&id=",id,sep="")
        rcv_doc <- GET(rcv_search, user_agent(UA))
        rcv_data <- xmlParse(content(rcv_doc, "text"))
        nodes_rcv <- getNodeSet(xmlRoot(rcv_data), "//ClinVarAccession")
        accs <-sapply(nodes_rcv, xmlGetAttr, "Acc")
        rcv_subset<-grep("^[R].*", accs)
        nodes_sig <-getNodeSet(xmlRoot(rcv_data), "//ClinicalSignificance")
        nodes_sig <-nodes_sig[rcv_subset]
        record_dates<-sapply(nodes_sig, xmlGetAttr, "DateLastEvaluated")
        if(list(NULL) %in% record_dates){
          nodes_sig<-nodes_sig[-which(sapply(record_dates, is.null)==TRUE)]
        }
        node_order<-rev(order(as.Date(sapply(nodes_sig, xmlGetAttr, "DateLastEvaluated"))))
        nodes_sig<-nodes_sig[node_order]
        nodes_title <- getNodeSet(xmlRoot(rcv_data), "//Title")
        nodes_title <- nodes_title[node_order]
        title<-xmlValue(nodes_title[[1]])
        rcvs[i]<-accs[rcv_subset][node_order][1]
        dess[i]<-xmlValue(nodes_sig[[1]][["Description"]])
        NM_num[i]<-strsplit(title,"[\\: ()]")[[1]][1]
        cchange_item<-strsplit(title,"[\\: ()]")[[1]][4]
        pchange_item<-strsplit(title,"[\\: ()]")[[1]][6]
        if(cchange_item==""){
          cchange[i]<-"-"
        }else{cchange[i]<-cchange_item}
        if(pchange_item==""){
          pchange[i]<-"-"
        }else{pchange[i]<-pchange_item}
        Sys.sleep(2.5)
      }
      output<-vcf[,c("Gene_Name","VMF")]
      output["NCBI Transcript ID"]=NM_num
      output["ClinVar Accession"]=rcvs
      output["Clinical Significance"]=dess
      output["SampleID"]=rep(sample_name[f],nrow(vcf))
      output["cDNA_change"]=cchange
      output["AA_change"]=pchange
      output$Disease <- as.factor(output$Gene_Name)
      levels(output$Disease) <- Disease.levels
    }
    #dp file calculation
    dp_file<-dp_ls[f]
    dp<-read.table(dp_file,sep="\t",skip=1)
    dp_low <- dp[which(dp[,4]<20),c(1,3,4)]
    if(nrow(dp_low)==0){
      output["QC"]=rep("PASS", nrow(output))}else{
        output["QC"]=rep("FAILED", nrow(output));
      end_position<-c(which(diff(dp_low[,2]) != 1),length(dp_low[,2]))
      start_position<-c(1,which(diff(dp_low[,2]) != 1)+1)
      region_pairs<-length(end_position)
      mean<-c()
      for(i in 1:region_pairs){
        region_depth<-dp_low[start_position[i]:end_position[i],3]
        region_mean<-mean(region_depth)
        mean[i]<-region_mean
      }
      chr<-dp_low[end_position,1]
      low_coverage_regions<-data.frame(chr=chr,start=c(dp_low[start_position,2]-1),end=c(dp_low[end_position,2]),mean=mean)
      output_csv<-paste(dir,"/",sample_name[f],"_low_coverage_regions.csv",sep="")
      write.csv(low_coverage_regions,file=output_csv,row.names=FALSE)
      }
  output_list[[f]]<-output
}
final_output<-do.call(rbind,output_list)
names(final_output)[names(final_output) == 'VMF']<-"AF"
final_output<-final_output[c("SampleID","Gene_Name","cDNA_change","AA_change","AF","Disease","NCBI Transcript ID","ClinVar Accession","Clinical Significance","QC")]
write.csv(final_output,file=final_output_csv,row.names=FALSE)
