require(data.table)
require(ggplot2)
require(limma)


#read data
breast.data<-fread("Masterarbeit/data/breast_cancer/GSE20685_breast_cancer.csv", sep=",")

#call used to create sample data.
#grep "Series_sample_id" GSE20685_series_matrix.txt > colnames.txt
#grep "Sample_characteristics_ch1" GSE20685_series_matrix.txt > design.txt 
colna<-fread("Masterarbeit/data/breast_cancer/colnames.txt", header = F)
names<-unlist(str_split(colna$V2, " "))
names<-names[names!=""]
design<-fread("Masterarbeit/data/breast_cancer/design.txt", header = F)
th<-unlist(sapply(colnames(design), function(x) grep("subtype", unlist(design[,c(x), with=F]), value=T)))
type<-data.table(names, th)
type$th<-sapply(type$th, function(x) gsub("subtype: type ", "", x))
#create annotation
annot<-breast.data[,.(sample, age.at.diagnosis, os_event, os_time)]

typ<-as.factor(type$th)
age<-as.factor(annot$age.at.diagnosis)

data.mat<-t(data.mat)


#transform data matrix colums = patients, rows= genes
data.mat<-breast.data[, c(5:ncol(breast.data)), with=F]
col.names<-rownames(data.mat)
data.mat<-as.data.table(t(data.mat))
colnames(data.mat)<-col.names


design<-model.matrix(~typ)
fit<-lmFit(data.mat, design)
fit <- eBayes(fit)
summary(decideTests(fit[,-1]))

topTable(fit,coef=4,n=20)

pr<-prcomp(data.mat, center = T, scale. = T)
g1<-ggplot(as.data.table(pr$rotation), aes(PC1, PC2))+geom_point()

