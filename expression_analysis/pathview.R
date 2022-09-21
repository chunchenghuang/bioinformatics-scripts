#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("pathview")

library(pathview)

#demo
data("gse16873.d")

#import data manually
#filename=system.file("extdata/gse16873.demo", package = "pathview")
#gse16873=read.delim(filename, row.names=1)

#subtract disease samples (DCIS) from control samples (HN) to get expression change
#gse16873.d=gse16873[,2*(1:6)]-gse16873[,2*(1:6)-1]

#load the sample related pathways from KEGG pathway analysis
data(demo.paths)

#full list of KEGG pathway ID/names if needed
#human pathway IDs (in the form of hsa+5 digits)
#hsa05321 
#"Inflammatory bowel disease (IBD)"

#KEGG species code

data(korg)
#head(korg)

#ktax.id  tax.id  kegg.code scientific.name           common.name                     entrez.gnodes kegg.geneid ncbi.geneid
#[1,] "T01001" "9606"  "hsa"     "Homo sapiens"            "human"                         "1"           "374659"    "374659"   
#[2,] "T01005" "9598"  "ptr"     "Pan troglodytes"         "chimpanzee"                    "1"           "474020"    "474020"   
#[3,] "T02283" "9597"  "pps"     "Pan paniscus"            "bonobo"                        "1"           "100989900" "100989900"
#[4,] "T02442" "9595"  "ggo"     "Gorilla gorilla gorilla" "western lowland gorilla"       "1"           "101125212" "101125212"
#[5,] "T01416" "9601"  "pon"     "Pongo abelii"            "Sumatran orangutan"            "1"           "100172878" "100172878"
#[6,] "T03265" "61853" "nle"     "Nomascus leucogenys"     "northern white-cheeked gibbon" "1"           "105739221" "105739221"
#ncbi.proteinid uniprot 
#[1,] "NP_001273380" "Q8N4P3"
#[2,] "XP_001140087" "Q1XHV8"
#[3,] "XP_003811308" NA      
#[4,] "XP_018886437" "G3QNH0"
#[5,] "NP_001125944" "Q5R9G0"
#[6,] "XP_012359712" "G1RK33"

data(bods)
data(gene.idtype.list)
#check your species id.type before executing pathview
#gene.idtype.list
#[1] "SYMBOL"       "GENENAME"     "ENSEMBL"      "ENSEMBLPROT"  "UNIGENE"      "UNIPROT"      "ACCNUM"       "ENSEMBLTRANS"
#[9] "REFSEQ"       "ENZYME"       "TAIR"         "PROSITE"      "ORF"    
#head(bods)
#package          species       kegg code id.type
#[1,] "org.Ag.eg.db"   "Anopheles"   "aga"     "eg"   
#[2,] "org.At.tair.db" "Arabidopsis" "ath"     "tair" 
#[3,] "org.Bt.eg.db"   "Bovine"      "bta"     "eg"   
#[4,] "org.Ce.eg.db"   "Worm"        "cel"     "eg"   
#[5,] "org.Cf.eg.db"   "Canine"      "cfa"     "eg"   
#[6,] "org.Dm.eg.db"   "Fly"         "dme"     "eg"

data(paths.hsa)
head(paths.hsa,3)

pv.out<-pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                 + species = "hsa", out.suffix = "gse16873", kegg.native = F,
                 + sign.pos = demo.paths$spos[i])