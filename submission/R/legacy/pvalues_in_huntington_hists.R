require(data.table)
require(ggplot2)

pvals_raw<-fread("/home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/huntington-gene-expression-DOWN.txt", header=T, sep="\t")

pvals<-as.data.table(unlist(pvals_raw[,2:ncol(pvals_raw)]))
colnames(pvals)<-c("Pvalue")

gpl<-ggplot(pvals, aes(x=Pvalue))+geom_histogram()+
  labs(x="P-value")+ggtitle("Histogram of p-values in original data")
plot(gpl)

#pvals_raw[,Mean:=lapply(.SD,mean), by = ID]
pvals_raw[,median:=lapply(.SD,median), by = ID]
#pvals_clean<-pvals_raw[which(pvals_raw$Mean!=1.0)]
pvals_clean<-pvals_raw[which(pvals_raw$median!=1.0)]
write.table(pvals_clean[,1:(ncol(pvals_clean)-1)], "/home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/huntington_no1.csv", sep = "\t", row.names = F, quote = F)

pg<-ggplot(pvals_clean, aes(x=Mean))+geom_histogram()
pg<-ggplot(pvals_clean, aes(x=median))+geom_histogram()
plot(pg)

pvals_clean2<-data.table(pval = unlist(pvals_clean[,2:(ncol(pvals_clean)-1)]))
gpl2<-ggplot(pvals_clean2, aes(x=pval))+geom_histogram()+
  labs(x="P-value")+ggtitle("Histogram of p-values in original data with samples with mean \n pval = 1 removed")
plot(gpl2)

