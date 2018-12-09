require(data.table)
require(ggplot2)
require(DESeq2)
require(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)


# get all the filenames from the directory
dirlist<- list.dirs("/home/anne/Documents/Master/MA/data/quants/fastq", recursive = F)
names<-sapply(basename(dirlist), function(x) gsub("_quant", "", x))
runs<-sapply(1:length(names), function(x) file.path(dirlist[x],"quant.sf"))

condition <-c(rep("control", 4), rep("HD",4))

samples<-data.table(sample=names, run = runs , condition = condition)
dir<-"~/Masterarbeit/data"
files <- file.path(dir, "salmon_out", samples$run)
names(files) <- sapply(1:length(condition), function(x) paste0(condition[x], x))

#sapply(samples$run, function(x) print(str_c("sed -i -E's/\.[0-9]\t/\t /g'", fixed=T)))

txdb<-EnsDb.Hsapiens.v86
tx2gene<-transcripts(txdb, return.type = "DataFrame")
tx2gene<-data.table(TXNAME = tx2gene$tx_id, GENEID = tx2gene$gene_id)
txi.salmon <- tximport(samples$run, type = "salmon", tx2gene = tx2gene)

ddsTxi<-DESeqDataSetFromTximport(txi.salmon, colData = samples, design = ~condition)
dds<-DESeq(ddsTxi)

resLFC <- lfcShrink(dds, coef=2)
res<-results(dds)
#head(res)
resOrdered <- res[order(res$pvalue),]
#head(resOrdered)

write.table(resOrdered, "/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out.tsv", sep = "\t" )

for(i in grep("HD", condition)){
  cont<-c(grep("HD", condition, invert = T),i)
  sample_sub<-samples[cont]
  
  txi.salmon <- tximport(sample_sub$run, type = "salmon", tx2gene = tx2gene)
  
  ddsTxi<-DESeqDataSetFromTximport(txi.salmon, colData = sample_sub, design = ~condition)
  dds<-DESeq(ddsTxi)
  
  resLFC <- lfcShrink(dds, coef=2)
  res<-results(dds)
  #head(res)
  resOrdered <- res[order(res$pvalue),]
  #head(resOrdered)
  
  write.table(resOrdered, paste0("/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub",i,".tsv"), sep = "\t" )
  
}



