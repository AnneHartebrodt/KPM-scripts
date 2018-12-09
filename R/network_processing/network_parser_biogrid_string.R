require(data.table)
require(ggplot2)
require(cowplot)

outdir<-"/home/anne/Masterarbeit/data/networks/StringDB/"

stringDb<-fread("/home/anne/Masterarbeit/data/networks/StringDB/9606.protein.links.v10.5.txt")
# only human proteins
stringDb$protein1<-sapply(stringDb$protein1, function(x) gsub("9606.", "", x))
stringDb$protein2<-sapply(stringDb$protein2, function(x) gsub("9606.", "", x))
#remove self loops
stringDb<-stringDb[protein1!=protein2]
copy<-stringDb
#stringDb<-copy

mapping<-fread("/home/anne/Masterarbeit/data/identifier_mappings/mart_export.txt")
mapping<-mapping[!is.na(`Gene stable ID`)]
mapping<-mapping[!is.na(`NCBI gene ID`)]
mapping<-mapping[!is.na(`Protein stable ID`)]

#stringDb<-merge(stringDb, unique(mapping[,.(`NCBI gene ID`, `Protein stable ID`)]), by.x="protein1", by.y="Protein stable ID")
#stringDb<-merge(stringDb, unique(mapping[,.(`NCBI gene ID`, `Protein stable ID`)]), by.x="protein2", by.y="Protein stable ID")

mapAndWrite<-function(stringDB, idtype, cutoff, filename){
  print(c(idtype, "Protein stable ID"))
stringDb<-merge(stringDb, unique(mapping[, c(idtype, "Protein stable ID"), with=F]), by.x="protein1", by.y="Protein stable ID")
stringDb<-merge(stringDb, unique(mapping[ , c(idtype, "Protein stable ID"), with=F]), by.x="protein2", by.y="Protein stable ID")
stringDb<-unique(stringDb)
# high confidence scores >0.7 -> 700 in file according to String FAQ
stringDb<-stringDb[combined_score>cutoff]
stringDb$interaction_type <- "pp"
stringDb<-stringDb[, c(4,6,5)]
write.table(stringDb, filename, col.names = F, row.names = F, quote = F,
            sep = "\t")

stringDb_flat<-data.table(nodes=unlist(stringDb[,1], stringDb[,3]))
count<-stringDb_flat[,.N, by=nodes]
print(nrow(count[N>1000]))
#ommited nodes: 14
gg_string<-ggplot(count, aes(x=N, fill=I("navyblue")))+
  geom_histogram(bins = 1000)+
  scale_y_log10()+
  coord_cartesian(xlim=c(0,1100))+
  ggtitle(paste0("Node Degree Distribution StringDB v.10.5 (cutoff ", cutoff,")"))+
  xlab("Node Degree")+
  ylab("Log10(#nodes)")
plot(gg_string)
ggsave(paste0(outdir, "/plots/stringDb_degree_distribution_",idtype, "_", cutoff, ".pdf"), gg_string)

#Some network stats:
nr_nodes<-length(unique(count$nodes))

sub<-stringDb[,c(3,2,1)]
colnames(sub)<-colnames(stringDb)
mm<-length(which(duplicated(rbind(stringDb, sub))))
nr_edges<-nrow(stringDb)-mm/2

cat(paste0(paste(idtype, cutoff, "nrNodes", nr_nodes, "nrEdges", nr_edges, sep = "\t"), "\n"), file=file.path(outdir, "logger.txt"), append = T)

}

names<-c("ENS", "NCBI")
ids<-c("Gene stable ID", "NCBI gene ID")

for(idtype in 1:length(ids)){
for(cutoff in c(700, 900)){
filename<-file.path(outdir, paste0("Homo_Sapiens_String_", names[idtype] , 
                 cutoff,".tsv"))
mapAndWrite(stringDB, ids[idtype], cutoff, filename)
}
}





biogrid<-fread("/home/anne/Masterarbeit/data/networks/biogrid/BIOGRID-ORGANISM-Homo_sapiens-3.4.163.tab2.txt")
biogrid_entrez<-biogrid[,2:3]
#remove self loops for KPM
biogrid_entrez<-biogrid_entrez[`Entrez Gene Interactor A`!=`Entrez Gene Interactor B`]
biogrid_entrez$interaction_type <-"pp"
biogrid_entrez<-biogrid_entrez[,c(1,3,2)]
biogrid_entrez<-unique(biogrid_entrez)
write.table(biogrid_entrez, "/home/anne/Masterarbeit/data/networks/biogrid/Homo_sapiens_biogrid.tsv", quote = F,
            sep = "\t", row.names = F, col.names = F)

#same for ENS ids.
# this inflates the network artificially, since the mapping is ambigous. Not sure if it makes sense
# to use as is.
biogrid_ENS<-merge(biogrid_entrez, unique(mapping[, c("Gene stable ID", "NCBI gene ID"), with=F]),
                   by.x="Entrez Gene Interactor A", by.y="NCBI gene ID")
biogrid_ENS<-merge(biogrid_ENS, unique(mapping[, c("Gene stable ID", "NCBI gene ID"), with=F]),
                   by.x="Entrez Gene Interactor B", by.y="NCBI gene ID", allow.cartesian = T)

write.table(biogrid_entrez, "/home/anne/Masterarbeit/data/networks/biogrid/Homo_sapiens_biogrid_ENS.tsv", quote = F,
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
