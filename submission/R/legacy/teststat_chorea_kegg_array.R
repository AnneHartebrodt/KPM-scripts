

require(data.table)
require(ggplot2)

chorea_genes<-fread("/home/anne/Documents/Master/MA/data/Kegg_Huntington_EntrezID.csv")
chorea_genes<-as.data.table(unique(unlist(chorea_genes)))


pvals_raw<-fread("/home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/huntington-gene-expression-DOWN.txt", header=T, sep="\t")

chor_p<-pvals_raw[ID %in% chorea_genes$V1]
chor_p<-chor_p[,2:39]
mean<-chor_p[,.(Mean = rowMeans(.SD))]

ggplot(as.data.table(unlist(chor_p)), aes(x=V1))+geom_histogram()

test<--2*sum(log(mean))

pchisq(test, 108, lower.tail = T)


