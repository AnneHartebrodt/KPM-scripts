require(data.table)
require(ggplot2)


#read data
breast.data<-fread("Masterarbeit/data/breast_cancer/GSE20685_breast_cancer.csv", sep=",")

#create annotation
annot<-breast.data[,.(sample, age.at.diagnosis, os_event, os_time)]

#transform data matrix colums = patients, rows= genes
data.mat<-breast.data[, c(5:ncol(breast.data)), with=F]
col.names<-rownames(data.mat)
data.mat<-as.data.table(t(data.mat))
colnames(data.mat)<-col.names

