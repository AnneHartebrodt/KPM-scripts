require(data.table)
require(ggplot2)
require(limma)
require(tximport)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
library(stringr)
require(edgeR)


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

for(i in which(samples$diagnosis=="HD")){
  #ind<-c(1:49,i)
  ind<-which(samples$age_class==samples$age_class[i])
  ind<-c(ind[ind<50],i)
  samples.2<-samples[ind,]
  txi.salmon.2<-list(txi.salmon$abundance[,ind], txi.salmon$counts[,ind], txi.salmon$length[,ind], txi.salmon$countsFromAbundance)
  names(txi.salmon.2)<-names(txi.salmon)
  t.test<-rowttests(txi.salmon.2$counts, as.factor(samples.2$diagnosis))
  t.test<-as.data.table(t.test)
  t.test$names<-rownames(txi.salmon.2$counts)
  t.test<-t.test[, c(3:4)]
  t.test$p.value<-p.adjust(t.test$p.value, method = "BH")
  colnames(t.test)<-c(i, "geneName")
  # txi.salmon.2<-list(txi.salmon$abundance[,ind], txi.salmon$counts[,ind], txi.salmon$length[,ind], txi.salmon$countsFromAbundance)
  # names(txi.salmon.2)<-names(txi.salmon)
  # 
  # dge <- DGEList(counts=txi.salmon.2$counts)
  # design<-model.matrix(~ factor(samples.2$diagnosis))
  # #design <- model.matrix(~ 0+factor(samples.2$diagnosis))
  # #colnames(design) <- c("HD", "NN")
  # #contrast.matrix <- makeContrasts(HD-NN, levels=design)
  # #design$NN<-as.integer(samples.2$diagnosis=="NN")
  # #design$HD<-as.integer(samples.2$diagnosis=="HD")
  # #design<-design[,2:3]
  # #keep<-txi.salmon.2$counts[,50]>10.00
  # #dge <- dge[keep,,keep.lib.sizes=FALSE]
  # dge <- calcNormFactors(dge, method = "TMM")
  # #logCPM <- cpm(dge, log=TRUE, prior.count=3)
  # logCPM <- voom(dge, design)
  # fit <- lmFit(logCPM, design)
  # fit <- eBayes(fit, trend=TRUE)
  # tt<-topTable(fit, coef=ncol(design), n=Inf)
  # tt$geneName<-row.names(tt)
  # tt<-tt[, c(5,7)]
  # colnames(tt)<-c(i, "geneName")
   all_res[[i]]<-t.test
}
DE_all<-all_res[[50]]
for(ind in 51:length(all_res)){
  DE_all<-merge(DE_all, all_res[[ind]], by="geneName")
}
DE_all<-as.data.table(DE_all)
write.table(DE_all, file.path(out.dir, "limma_results.tsv"))

DE_all_sig<-DE_all<0.05
DE_all_sig[DE_all_sig==TRUE]<-1
DE_all_sig[DE_all_sig==FALSE]<-0
DE_all_sig<-as.data.table(DE_all_sig)
DE_all_sig$geneName<-DE_all$geneName
which(rowSums(DE_all_sig[,2:ncol(DE_all_sig)])>10)
DE_all_sig[is.na(DE_all_sig)]<-0
write.table(DE_all_sig, file.path(out.dir, "limma_indicator_matrix.tsv"), sep="\t", quote = F, row.names = F)
saveRDS(all_res, file.path(out.dir, "limma_all.rds"))

count<-txi.salmon$counts/txi.salmon$length
count<-txi.salmon$abundance
write.table(count,file.path(out.dir, "counts.txt"), quote = F, col.names = F)



dge <- DGEList(counts=txi.salmon$counts)
l=as.vector(unlist(txi.salmon$counts))
ggplot(as.data.table(l), aes(x=l))+geom_histogram()+scale_y_log10()
design<-as.data.table(samples$diagnosis)
design <- model.matrix(~ factor(samples$diagnosis))
#colnames(design) <- c("HD", "NN")
#contrast.matrix <- makeContrasts(HD-NN, levels=design)
#design$HD<-as.integer(samples$diagnosis=="HD")
#design$NN<-as.integer(samples$diagnosis=="NN")
#keep<-rowSums(txi.salmon$counts[,50:69]>10.00)>18
#design<-design[,2:3]
keep <- filterByExpr(dge, design = design, min.count=10, min.total.count=10)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "upperquartile")
#logCPM <- cpm(dge, log=TRUE, prior.count=3)
logCPM <- cpm(dge, log=TRUE,prior.count = 3)
logCPM<-voom(txi.salmon$abundance, design = design, normalize.method = "quantile")
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
tt<-topTable(fit, coef=ncol(design), n=Inf)
tt$geneName<-row.names(tt)
tt<-as.data.table(tt)
tt[adj.P.Val<0.001]
tt[logFC>5]
tt[,lo:=exp(B)]
tt[,prob:=lo/(1+lo)]
tt[prob>0.99]
