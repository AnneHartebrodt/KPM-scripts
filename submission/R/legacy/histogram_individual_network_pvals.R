require(data.table)
require(ggplot2)

setwd("/home/anne/Documents/Master/MA/Testing/out/dist/")
files<-dir("/home/anne/Documents/Master/MA/Testing/out/dist/")

files<-files[grep("^92", files, perl = T)]

ff<-sapply(files, function(x) fread(x))
df<-as.data.frame(ff)
dd<-melt(df)
ggplot(dd[1:1000,], aes(value))+geom_histogram()+facet_wrap(facets = "variable")+
  ggtitle("Individual p-value distributions of nodes")

ff<-data.table(ff=unlist(ff))
ggplot(ff, aes(x=ff))+geom_histogram()+
ggtitle("Distribution of p-values for subnetworks of size 46")
