require(data.table)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
require(graph)
require(ggplot2)

txdb<-EnsDb.Hsapiens.v86
tx2gene<-transcripts(txdb, return.type = "DataFrame")

network<-fread("~/Masterarbeit/data/networks/human_StringDB.txt")

genes<-fread("~/Masterarbeit/data/Kegg_Huntington_EntrezID.csv")
genes<-data.table(name=unique(unlist(genes)))
map<-fread("~/Masterarbeit/data/human_Entrez_ENSG.tsv")
map<-map[!is.na(map$`NCBI gene ID`)]

huntington_ENSG<-map[map$`NCBI gene ID` %in% genes$name]

other_ENSG<-map[!(map$`Gene stable ID` %in% huntington_ENSG$`Gene stable ID`)]

study<-c(rep("control", 4), rep("sample", 4))
background<-sapply(other_ENSG$`Gene stable ID`, function(x) runif(8))


rr<-as.numeric(sapply(runif(100000), function(x) -log(1-(1-exp(-80))*x)/80))
ru<-runif(1000000)
ra<-c(rr,ru)
ra <- data.frame(x=ra, a=p.adjust(ra,method = "BH" ))
ggplot(ra, aes(x))+geom_histogram()
ggplot(ra,aes(x=a))+geom_histogram()



limma<-fread("~/Documents/Master/MA/data/HDarray/processed/pvalue_diff_exp_raw.txt")
limma<-limma[,2:ncol(limma), with=F]
lim<-data.frame(pval=unlist(limma))
ggplot(lim, aes(x=pval))+geom_histogram()

DES<-DESeq
DESeq<-fread("/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub_merged_unad.tsv")
des<-data.frame(pa=unlist(DESeq[,2:ncol(DESeq), with=F]))
ggplot(des, aes(x=pa))+geom_histogram()
