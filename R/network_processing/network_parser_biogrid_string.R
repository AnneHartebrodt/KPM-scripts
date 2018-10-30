require(data.table)
require(ggplot2)
require(cowplot)

stringDb<-fread("/home/anne/Masterarbeit/data/networks/StringDB/9606.protein.links.v10.5.txt")
# only human proteins
stringDb$protein1<-sapply(stringDb$protein1, function(x) gsub("9606.", "", x))
stringDb$protein2<-sapply(stringDb$protein2, function(x) gsub("9606.", "", x))
#remove self loops
stringDb<-stringDb[protein1!=protein2]
copy<-stringDb
#stringDb<-copy

mapping<-fread("/home/anne/Masterarbeit/data/networks/mart_export.txt")
mapping<-mapping[!is.na(`Gene stable ID`)]
mapping<-mapping[!is.na(`NCBI gene ID`)]
mapping<-mapping[!is.na(`Protein stable ID`)]

#stringDb<-merge(stringDb, unique(mapping[,.(`NCBI gene ID`, `Protein stable ID`)]), by.x="protein1", by.y="Protein stable ID")
#stringDb<-merge(stringDb, unique(mapping[,.(`NCBI gene ID`, `Protein stable ID`)]), by.x="protein2", by.y="Protein stable ID")

stringDb<-merge(stringDb, unique(mapping[,.(`Gene stable ID`, `Protein stable ID`)]), by.x="protein1", by.y="Protein stable ID")
stringDb<-merge(stringDb, unique(mapping[,.(`Gene stable ID`, `Protein stable ID`)]), by.x="protein2", by.y="Protein stable ID")
stringDb<-unique(stringDb)
# high confidence scores >0.7 -> 700 in file according to String FAQ
stringDb<-stringDb[combined_score>700]
stringDb$interaction_type <- "pp"
stringDb<-stringDb[, c(4,6,5)]

write.table(stringDb,"/home/anne/Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS.tsv", col.names = F, row.names = F, quote = F,
            sep = "\t")

stringDb_flat<-data.table(nodes=c(stringDb$`NCBI gene ID.x`, stringDb$`NCBI gene ID.y`))
count<-stringDb_flat[,.N, by=nodes]
nrow(count[N>1000])
#ommited nodes: 14
gg_string<-ggplot(count, aes(x=N, fill=I("navyblue")))+
  geom_histogram(bins = 1000)+
  scale_y_log10()+
  coord_cartesian(xlim=c(0,1100))+
  ggtitle("Node Degree Distribution in StringDB Network (10.5)")+
  xlab("Node Degree")+
  ylab("Log10(#nodes)")
plot(gg_string)
ggsave("/home/anne/Masterarbeit/data/networks/StringDB/plots/stringDb_degree_distribution.pdf", gg_string)

#Some network stats:
nr_nodes<-length(unique(count$nodes))

sub<-stringDb[,c(3,2,1)]
colnames(sub)<-colnames(stringDb)
mm<-length(which(duplicated(rbind(stringDb, sub))))
nr_edges<-nrow(stringDb)-mm/2
#14247
#329611

biogrid<-fread("/home/anne/Masterarbeit/data/networks/biogrid/BIOGRID-ORGANISM-Homo_sapiens-3.4.163.tab2.txt")
biogrid_entrez<-biogrid[,2:3]
#remove self loops for KPM
biogrid_entrez<-biogrid_entrez[`Entrez Gene Interactor A`!=`Entrez Gene Interactor B`]
biogrid_entrez$interaction_type <-"pp"
biogrid_entrez<-biogrid_entrez[,c(1,3,2)]
biogrid_entrez<-unique(biogrid_entrez)
write.table(biogrid_entrez, "/home/anne/Masterarbeit/data/networks/biogrid/Homo_sapiens_biogrid.tsv", quote = F,
            sep = "\t", row.names = F, col.names = F)


biogrid_flat<-data.table(nodes=c(biogrid_entrez$`Entrez Gene Interactor A`, biogrid_entrez$`Entrez Gene Interactor B`))
count<-biogrid_flat[,.N, by=nodes]
nrow(count[N>1000])
# ommited nodes: 14
gg_bio<-ggplot(count, aes(x=N,fill=I("navyblue")))+
  geom_histogram(bins = 1000)+
  scale_y_log10()+
  coord_cartesian(xlim=c(0,1000))+
  ggtitle("Node Degree Distribution in Biogrid Network (V3.4.163)")+
  xlab("Node Degree")+
  ylab("Log10(#nodes)")
plot(gg_bio)
ggsave("/home/anne/Masterarbeit/data/networks/biogrid/plots/biogrid_degree_distribution.pdf", gg_bio)

# biogrid stats
nr_nodes<-length(unique(biogrid_flat$nodes))

sub<-biogrid_entrez[,c(3,2,1)]
colnames(sub)<-colnames(biogrid_entrez)
mm<-length(which(duplicated(rbind(biogrid_entrez, sub))))
nr_edges<-nrow(biogrid_entrez)-mm/2
#22792
#310706
