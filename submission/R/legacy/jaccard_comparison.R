file.location<-"/home/anne/Masterarbeit/Test_pipeline/results/String0.0/"

res<-NULL
res<-data.table(V1=c(), V2=c(), V3=c(),V4=c(),V5=c())
for(d in dir(file.location)){
print(d)
 dat<-fread(paste0(file.location, d, "/", "jacccards.txt"))
 dat$method<-gsub("/home/anne/Masterarbeit/Test_pipeline/results/String0.0/StringDB_", "", d)
 #dat$method<-gsub("jacccards.txt", "", dat$method)
 res<-rbind(res, dat)
}

res$method<-res$method[order(res$method, decreasing = T)]
g1<-ggplot(res, aes(y=V1, x=method))+geom_boxplot(fill ="#3182bd")+
  scale_x_discrete(labels=c("StringDB_normDegSum"="Node Degree\nNormalized Sum",
                            "StringDB_normSum"="Network Size\nNormalized Sum", "StringDB_sum"="Sum"))+
  xlab("Score Aggregation Method")+
  ylab("Jaccard index")+
  ggtitle("Comparison of Different Score Aggregation Methods Using Artificial Data")+
  theme(plot.title = element_text(size=18))
plot(g1)

unique(dat$method)
