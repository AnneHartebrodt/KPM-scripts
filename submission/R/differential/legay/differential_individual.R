require(data.table)
require(ggplot2)
require(DESeq2)
require(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
library(stringr)
require(optparse)
require(cowplot)
require(gProfileR)
#require(ggrepel)

## Remove version numbers of Transcripts form files
## cd "/home/anne/Documents/Master/MA/data/dataHD/quants/"
## bash ../../../scripts/salmon/cut_name.txt 

# option_list = list(
#   make_option(c("-q", "--quant"), type="character", default=NULL, 
#               help="quant dir", metavar="character"),
#   make_option(c("-m", "--map"), type="character", default=NULL, 
#               help="identifier mapping", metavar="character"),
#   make_option(c("-r", "--run"), type="character", default=NULL, 
#               help="runtable", metavar="character"),
#   make_option(c("-o", "--out"), type="character", default=NULL, 
#               help="output file name [default= %default]", metavar="character"),
#   make_option(c("-t", "--kpmin"), type="character", default=NULL, 
#               help="kpm input directory", metavar="character"),
#   make_option(c("-n", "--network"), type="character", default=NULL, 
#               help="input network", metavar="character"),
#   make_option(c("-p", "--pathway"), type="character", default=NULL, 
#               help="input network", metavar="character"),
#   make_option(c("-a", "--map2"), type="character", default=NULL, 
#               help="input network", metavar="character")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# 
# runtable<-opt$run
# quants.dir<-opt$quant
# out.dir<-opt$out
# map.file<-opt$map
# network<-opt$network
# pathway_id<-opt$pathway
# kpm.in<-opt$kpmin
# maap<-opt$map2

runtable<-"~/Masterarbeit/data/salmon/SraRunTable(2).txt"
map.file<-"~/Masterarbeit/data/identifier_mappings/mart_export.txt"
quants.dir<-"/home/anne/Documents/Master/MA/data/dataHD/quants/"
out.dir<-"/home/anne/Documents/Master/MA/differential_out/dataHD/DESEQ_ind_new/"
network<-"~/Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_NCBI700.tsv"
pathway_id<-"hsa05016"
map.file<-"~/Masterarbeit/data/identifier_mappings/mart_export.txt"
dir.create(out.dir, recursive = T)

##Transcript mapping
txdb<-EnsDb.Hsapiens.v86
tx2gene<-transcripts(txdb, return.type = "DataFrame")
tx2gene<-data.table(TXNAME = tx2gene$tx_id, GENEID = tx2gene$gene_id)


#tx2gene<-merge(tx2gene, mapping, by.x="GENEID", by.y="Gene stable ID")
#tx2gene<-tx2gene[,c("TXNAME","NCBI gene ID")]
#colnames(tx2gene)<-c("TXNAME","GENEID" )

# get all the filenames from the directory
dirlist<- list.dirs(quants.dir, recursive = F)
names<-sapply(basename(dirlist), function(x) gsub("_quant", "", x))
runs<-sapply(1:length(names), function(x) file.path(dirlist[x],"quant.sf"))



##Read and create sample information
age<-function(a){
  if(a<=45){
    return (1)
  }
  else if(a>45 & a<=60){
    return (2)
  }
  else if(a>60 & a<=75){
    return (3)
  }
  else{
    return(4)
  }
}
diag<-function(a){
  if(a=="Huntington's Disease"){
    return ("HD")
  }
  else{
    return("NN")
  }
}

##create annotation, read file, recode age and diagnosis to usable values
c<-fread(runtable)
c<-c[, c("Run", "age_of_death", "age_of_onset", "diagnosis", "rin")]
c$age_class<-sapply(c$age_of_death, function(x) age(x))
samples<-data.table(sample=names, run = runs)
ex<-sapply(samples$run, function(x) file.exists(x))
samples<-samples[ex]
samples<-merge(samples, c, by.x="sample", by.y="Run")
samples<-as.data.frame(samples)
rownames(samples)<-samples$sample
samples$diagnosis<-sapply(samples$diagnosis, function(x) diag(x))
samples$diagnosis<-as.factor(samples$diagnosis)
samples$age_class<-as.factor(samples$age_class)


##import count data
txi.salmon <- tximport(samples$run, type = "salmon", tx2gene = tx2gene)
samples<-samples[,c(5,6,7)]

all_res<-list()
DE_matrix<-NULL

for(i in which(samples$diagnosis=="HD")){
  ind<-c(1:49,i)
samples.2<-samples[ind,]
txi.salmon.2<-list(txi.salmon$abundance[,ind], txi.salmon$counts[,ind], txi.salmon$length[,ind], txi.salmon$countsFromAbundance)
names(txi.salmon.2)<-names(txi.salmon)

#run DE analysis
ddsTxi<-DESeqDataSetFromTximport(txi.salmon.2, colData = samples.2, design = as.formula("~diagnosis"))
ddsTxi$diagnosis<-relevel(ddsTxi$diagnosis, ref="NN")
dds<-DESeq(ddsTxi)
resLFC <- lfcShrink(dds, coef=2)
res<-results(dds)
all_res[[i]]<-res
resOrdered <- res[order(res$padj),]
resOrdered$geneName<-row.names(resOrdered)
resOrdered<-as.data.table(resOrdered)

#Correct for additional source for multiple testing
siglevel<-0.01/20
resOrdered$sig<-as.integer(resOrdered$padj<siglevel)

if(is.null(DE_matrix)){
  re<-resOrdered[,.(geneName, sig)]
  colnames(re)<-c("geneName", i)
  DE_matrix<-re
}
else{
  re<-resOrdered[,.(geneName, sig)]
  colnames(re)<-c("geneName", i)
  DE_matrix<-merge(DE_matrix, re, by="geneName")
}
}                 
  
fwrite(DE_matrix, file=file.path(out.dir, "individual_pvals_diag.tsv"), sep="\t", row.names = F)
DE_matrix[is.na(DE_matrix)]<-0
fwrite(DE_matrix, file=file.path(out.dir, "individual_pvals_complete_diag.tsv"), sep="\t", row.names = F)

#save differential expression analyses for later use
saveRDS(all_res,file.path(out.dir,"all_res_diag.rds"))
dat<-fread("~/Masterarbeit/differential_out/dataHD/DESEQ_ind/individual_pvals_complete_diag.tsv")

#res<-readRDS(file.path(out.dir,"all_res.rds"))


#which(rowSums(dat[, 2:ncol(dat)])>1)
# siglevel<-0.05
# DE_matrix<-NULL
# for(res in 1:length(res_all)){
#   if(is.null(res_all[[res]])){
#     next()
#   }
#   else{
#     if(is.null(DE_matrix)){
#       print(res)
#       re<- as.data.table(res_all[[res]]@listData)
#       re$geneName<-row.names(res_all[[res]])
#       
#       re$sig<-as.integer(re$padj<siglevel)
#     re<-re[,.(geneName, sig)]
#     colnames(re)<-c("geneName",res )
#     DE_matrix<-re
#     }
#     else{
#       print(res)
#       re<- as.data.table(res_all[[res]]@listData)
#       re$geneName<-row.names(res_all[[res]])
#       
#       re$sig<-as.integer(re$padj<siglevel)
#       re<-re[,.(geneName, sig)]
#       colnames(re)<-c("geneName",res )
#     DE_matrix<-merge(DE_matrix, re, by="geneName")
#     }
#   }
# }

