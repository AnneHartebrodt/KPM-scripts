library(BioNet)
library(DLBCL)
require(data.table)

#create graphNEL of Network
dat<-fread("Masterarbeit/data/networks/StringDB/Homo_Sapiens_String.tsv")
dat$V1<-as.character(dat$V1)
dat$V3<-as.character(dat$V3)
dat<-dat[!(V1==V2)]
g<-graphNEL()
print(Sys.time())
g<-addNode(node = unique(c(dat$V1, dat$V3)), g)
g<-addEdge(as.vector(dat$V1), as.vector(dat$V3), g, 1)
print(Sys.time())

pvana<-fread("Masterarbeit/data/DESeq_out/DESeq_Pvals_Entrez.tsv")
pvalues<-fread("Masterarbeit/data/DESeq_out/DESeq_out.tsv")
pvalue<-as.numeric(pvalues$pvalue)
names(pvalue)<-pvana$V1
pvalue<-pvalue[!is.na(pvalue)]


g<-subGraph(names(pvalue), g)

fb <- fitBumModel(pvalue, plot = TRUE)
scores <- scoreNodes(g, fb, fdr = 0.001)

print(Sys.time())
module<-runHeinz(g, scores)
module <- runFastHeinz(g, scores)
print(Sys.time())

plotModule(module, scores = scores)
