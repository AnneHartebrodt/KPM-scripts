require(data.table)
require(ggplot2)
dnames<-dir("/home/anne/Masterarbeit/Testing/toydata/")
grep(dnames, "^0.3")
dnames[startsWith(dnames, "0.3")][1:5]

goldstandard<-params$goldstandard
goldstandard<-"~/Masterarbeit/data/toydata/expression_type.tsv"

gold<-fread(goldstandard)
gold<-gold[V2=="differential"]
gold<-gold[,1]
colnames(gold)<-c("gene")


jaccard<-function(g){
  sub1<-f1[V1==g,2:4] 
  sub1<-sub1[order(sub1$V2)]
  su<-unique(c(sub1$V2, sub1$V4))
  intersect_1<-length(which(su %in% gold$gene))
  union_1<-length(unique(c(su, gold$gene)))
  jaccard<-intersect_1/union_1
  return(jaccard)
}



plots<-list()
dd<-fread(paste0("/home/anne/Masterarbeit/Testing/toydata/",dnames[1], "/teststats.txt"))  
colnames(dd)<-c("name", "value")
basedata<-paste0("/home/anne/Masterarbeit/data/toydata/toy_general_mean",dnames[startsWith(dnames, "0.3")][1], ".tsv")
ddd<-fread(basedata)
plots[[1]]<-ggplot(ddd, aes(x=V2))+geom_histogram()


f1 <-fread(file.path("/home/anne/Masterarbeit/Testing/toydata/",dnames[1],"file.graph"))
f1[,.N, by = V1]
jaccard_dists<-sapply(unique(f1$V1), function(g) jaccard(g))
jaccard_dists<-as.data.table(jaccard_dists)
jaccard_dists$name<-dnames[1]

jj<-jaccard_dists

for(i in 2:125){
d<-fread(paste0("/home/anne/Masterarbeit/Testing/toydata/",dnames[i], "/teststats.txt"))
dd<-merge(dd, d, by.x = "name", by.y = "V1")
colnames(dd)[i+1]<-paste0("value", i)

basedata<-paste0("/home/anne/Masterarbeit/data/toydata/toy_general_mean",dnames[i], ".tsv")
ddd<-fread(basedata)
plots[[i]]<-ggplot(ddd, aes(x=V2))+geom_histogram()

if(file.exists(paste0("/home/anne/Masterarbeit/Testing/toydata/",dnames[i],"/file.graph"))){
  #print(file.exists(paste0("/home/anne/Masterarbeit/Testing/toydata",dnames[i],"file.graph")))
f1 <-fread(file.path("/home/anne/Masterarbeit/Testing/toydata",dnames[i],"file.graph"))
f1[,.N, by = V1]

jaccard_dists<-sapply(unique(f1$V1), function(g) jaccard(g))
jaccard_dists<-as.data.table(jaccard_dists)
jaccard_dists$name<-dnames[i]


jj<-rbind(jj, jaccard_dists)
}
}

ggplot(jj, aes(y=jaccard_dists, x=name))+geom_boxplot()

dd<-melt(dd, id.vars = "name")
ggplot(dd, aes(x=name,y=value, color = variable))+geom_line()

ggplot(jaccard)
