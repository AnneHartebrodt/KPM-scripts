require(data.table)
require(ggplot2)


## Load data: networks and mapping files
net <- fread("~/Masterarbeit/data/networks/Biog.sif")
map<-fread("~/Masterarbeit/data/human_Entrez_ENSG.tsv")
map<-map[!is.na(map$`NCBI gene ID`)]

m<-7
# sample genes, that form the subnetworks
subnets<-list()
for (i in 1:m){
seeds<-sample(net$V1, 1)
subnet<-net[V1 %in% seeds]
subnet<-subnet[1:5]
seeds<-net[V1 %in% subnet$V3]
subnet<-rbind(subnet, seeds)
subnet<-cbind(subnet, rep(i, nrow(subnet)))
subnets[[i]]<-subnet
}

# remove genes, that overlap between the networks
for(i in 1:(m-1)){
  for(j in (i+1):m){
    v1<-which(subnets[[i]]$V1 %in% subnets[[j]]$V1)
    v12<-which(subnets[[i]]$V1 %in% subnets[[j]]$V3)
    v2<-which(subnets[[i]]$V3 %in% subnets[[j]]$V3)
    subnets[[i]]<-subnets[[i]][-unique(c(v1,v12,v2))]
  }
}


subnet<-rbindlist(subnets)

differential<-unique(c(subnet$V1, subnet$V3))
nondifferential<-unique(map$`NCBI gene ID`[!(map$`NCBI gene ID` %in% differential)])

for(mn in c(0.6)){
  
for(sd1 in seq(from = 0.01, to= 0.1, by= 0.02)){
diff<-lapply(differential, function(x) abs(rnorm(8,mean = 0.00025, sd = sd1)))
diff<-as.data.table(diff)
diff<-as.data.table(t(diff))

for(sd2 in c(0.1)){
nodiff<-lapply(nondifferential, function(x) abs(rnorm(8, mean = mn, sd = sd2)))
nodiff<-as.data.table(nodiff)
nodiff<-as.data.table(t(nodiff))

#write files
toydata<-rbind(diff, nodiff)
toydata$gene_name<-c(differential, nondifferential)
toydata<-data.table(toydata$gene_name, toydata[,1:8])

type<-data.table(gene = toydata$V1, 
                 type = c(rep("differential", length(differential)), 
                          rep("nondifferential", length(nondifferential))))
filename_samplewise<-paste0("/home/anne/Documents/Master/MA/data/toydata/toy_samplewise_mean", mn, "_sd1",sd1, "_sd2", sd2,".tsv")
write.table(toydata, filename_samplewise, 
            sep = "\t",
            col.names = F,
            row.names =F)   
filename_general<-paste0("/home/anne/Documents/Master/MA/data/toydata/toy_general_mean", mn, "_sd1",sd1, "_sd2", sd2,".tsv")
write.table(toydata[, 1:2], filename_general, 
            sep = "\t", 
            col.names = F, 
            row.names =F)
diagnostic_call<-paste("-numProc=1", paste0("-matrix1=",
                        filename_samplewise),
        "-datasetsFile=/home/anne/Documents/Master/MA/Testing/datasets.txt",
"-summaryFile=harkan1.txt -combineOp=OR -pathwaysStatsFile=harkan2.txt", 
paste0("-resultsDir=/home/anne/Documents/Master/MA/Testing/toydata/", mn, "_sd1",sd1, "_sd2", sd2 ),
"-graphFile=/home/anne/Documents/Master/MA/data/networks/Biog.sif", 
"-geneStatsFile=harkan3.txt -K=5 -L1=2 -maxsolutions=5", 
"-perturbation_technique=edgerewire -algo=FDR -strategy=FDR -Umove_bens -comparator=LET", 
"-significance_level=0.05",  
paste0("-L1_pvalues=",filename_general), "-use_double -mfHeader", 
"-validation_file=/home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/COAD-VAL-ENTREZ.txt",  
"-L1_pvaluecutoff=0.05 -randomized_graph_file=/home/anne/Documents/Master/MA/data/networks/Biog_randomized.sif",  
"-aggregation_method=sum")
diagnostic_file<-filename_general<-paste0("/home/anne/Documents/Master/MA/Testing/toydata/call/call", mn, "_sd1",sd1, "_sd2", sd2,".txt")


fileConn<-file(diagnostic_file)
writeLines(diagnostic_call, fileConn)
close(fileConn)


}
}
}

write.table(type, paste0("~/Masterarbeit/data/toydata/expression_type.tsv"),
            sep = "\t", 
            col.names = F , 
            row.names = F)




diff<-data.table(diff=unlist(diff))
nodiff<-data.table(nodiff=unlist(nodiff))

gg<-ggplot(diff, aes(x=diff))+geom_histogram()+geom_histogram(data=nodiff, aes(x=nodiff))+
  ggtitle("Toydata p-values")+xlab("p-value")+ylab("Count")
plot(gg)
ggsave(plot = gg,"~/Masterarbeit/data/toydata/plots/toydata.pdf")


# Next attempt

diff<-lapply(differential, function(x) abs(rnorm(8,mean = 0.025, sd = 0.005)))
diff<-as.data.table(diff)
diff<-as.data.table(t(diff))

  nodiff<-lapply(nondifferential, function(x) abs(runif(8, 0,1)))
  nodiff<-as.data.table(nodiff)
  nodiff<-as.data.table(t(nodiff))
  
  #write files
  toydata<-rbind(diff, nodiff)
  toydata$gene_name<-c(differential, nondifferential)
  toydata<-data.table(toydata$gene_name, toydata[,1:8])
  
  type<-data.table(gene = toydata$V1, 
                   type = c(rep("differential", length(differential)), 
                            rep("nondifferential", length(nondifferential))))
  filename_samplewise<-paste0("/home/anne/Documents/Master/MA/data/toydata/toy_samplewise_unif.tsv")
  write.table(toydata, filename_samplewise, 
              sep = "\t",
              col.names = F,
              row.names =F)   
  filename_general<-paste0("/home/anne/Documents/Master/MA/data/toydata/toy_general_mean_unif.tsv")
  write.table(toydata[, 1:2], filename_general, 
              sep = "\t", 
              col.names = F, 
              row.names =F)
  write.table(type, paste0("~/Masterarbeit/data/toydata/expression_type_unif.tsv"),
              sep = "\t", 
              col.names = F , 
              row.names = F)
hist(toydata$V2)
