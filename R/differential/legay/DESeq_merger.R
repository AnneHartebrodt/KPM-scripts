require(data.table)

files<-list()

for(i in 5:8){
  df<-fread(paste0("/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub",i,".tsv"), header = F)
  colnames(df)<-c("geneName","baseMean", "log2FoldChange" ,"lfcSE","stat","pvalue" ,"padj" )
  files[[i-4]]<-df
  }
df<-

resultFrame<-files[[1]][,c("geneName", "padj")]
resultFrameUnad<-files[[1]][,c("geneName", "pvalue")]
for(i in 2:length(files)){
  resultFrame<-merge(resultFrame,files[[i]][,c("geneName", "padj")], by = "geneName", all=T)
  resultFrameUnad<-merge(resultFrameUnad,files[[i]][,c("geneName", "pvalue")], by = "geneName", all=T)
}
colnames(resultFrame)<-c("geneName",sapply(1:length(files), function(x) paste0("HD", x)))
colnames(resultFrameUnad)<-c("geneName",sapply(1:length(files), function(x) paste0("HD", x)))
                         

write.table(resultFrame, "/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub_merged.tsv", sep="\t", row.names = F)  
write.table(resultFrameUnad, "/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub_merged_unad.tsv", sep="\t", row.names = F)  

maps<-fread("/home/anne/Documents/Master/MA/data/human_Entrez_ENSG.tsv")
maps<-maps[!is.na(maps$`NCBI gene ID`)]

mapping<-maps[`Gene stable ID` %in% resultFrame$geneName]


resEntrez<-merge(resultFrame, mapping, by.x = "geneName", by.y = "Gene stable ID")
resUnadEntrez<-merge(resultFrameUnad, mapping, by.x = "geneName", by.y = "Gene stable ID")
resEntrez<-resEntrez[,c(6,2:5), with=F]
colnames(resEntrez)<-c("geneName",sapply(1:length(files), function(x) paste0("HD", x)))
resultUnadEntrez<-resUnadEntrez[,c(6,2:5), with=F]
colnames(resultUnadEntrez)<-c("geneName",sapply(1:length(files), function(x) paste0("HD", x)))

resEntrez[is.na(resEntrez)]<-1
write.table(resEntrez, "/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub_merged_entrez.tsv", sep="\t", row.names = F)  
write.table(resUnadEntrez, "/home/anne/Documents/Master/MA/data/DESeq_out/DESeq_out_sub_merged_unad_entrez.tsv", sep="\t", row.names = F)  

