require(data.table)
require(ggplot2)
require(DESeq2)
require(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
library(stringr)

##Transcript mapping
txdb<-EnsDb.Hsapiens.v86
tx2gene<-transcripts(txdb, return.type = "DataFrame")
tx2gene<-data.table(TXNAME = tx2gene$tx_id, GENEID = tx2gene$gene_id)


# get all the filenames from the directory
dirlist<- list.dirs("/home/anne/Documents/Master/MA/data/salmon/quants", recursive = F)
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

c<-fread("Masterarbeit/data/salmon/SraRunTable(2).txt")
c<-c[, c("Run", "age_of_death", "age_of_onset", "diagnosis", "rin")]
c$age_class<-sapply(c$age_of_death, function(x) age(x))

samples<-data.table(sample=names, run = runs)
ex<-sapply(samples$run, function(x) file.exists(x))
samples<-samples[ex]

samples<-merge(samples, c, by.x="sample", by.y="Run")
samples<-as.data.frame(samples)
rownames(samples)<-samples$sample
samples<-samples[,c(5,6,7)]
samples$diagnosis<-sapply(samples$diagnosis, function(x) diag(x))
samples$diagnosis<-as.factor(samples$diagnosis)
samples$age_class<-as.factor(samples$age_class)
##import count data
txi.salmon <- tximport(samples$run, type = "salmon", tx2gene = tx2gene)


ddsTxi<-DESeqDataSetFromTximport(txi.salmon, colData = samples, design = as.formula("~age_class+rin+diagnosis"))
dds<-DESeq(ddsTxi)

resLFC <- lfcShrink(dds, coef=2)
res<-results(dds)
resOrdered <- res[order(res$padj),]

write.table(resOrdered, "/home/anne/Documents/Master/MA/data/differential_out/HD/DESEQ/results.tsv", sep = "\t" )

# for(i in grep("HD", condition)){
#   cont<-c(grep("HD", condition, invert = T),i)
#   sample_sub<-samples[cont]
#   
#   txi.salmon <- tximport(sample_sub$run, type = "salmon", tx2gene = tx2gene)
#   
#   ddsTxi<-DESeqDataSetFromTximport(txi.salmon, colData = sample_sub, design = ~condition)
#   dds<-DESeq(ddsTxi)
#   
#   resLFC <- lfcShrink(dds, coef=2)
#   res<-results(dds)
#   #head(res)
#   resOrdered <- res[order(res$pvalue),]
#   #head(resOrdered)
#   
#   write.table(resOrdered, paste0("/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub",i,".tsv"), sep = "\t" )
#   
# }
