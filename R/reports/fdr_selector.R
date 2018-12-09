require(data.table)
require(ggplot2)

fdr<-fread("Masterarbeit/batch_test/result_true/biogrid_clu/fdr.tsv")
anno<-fread("Masterarbeit/batch_test/result_true/biogrid_clu/annotation.tsv")
eval<-fread("Masterarbeit/evaluation.tsv")


fdr<-merge(fdr, anno, by.x="V4", by.y="orig")

list_params<-sapply(colnames(anno), function(x) unique(anno[,x, with=FALSE]))
names(list_params)<-colnames(anno)
list_params<-list_params[-12]


select<-function(fdr, list_params){
fdr<-fdr[Network %in% list_params$Network & 
           Aggr_meth %in% list_params$Aggr_meth &
           FDR %in% list_params$FDR &
           Perc_perm %in% list_params$Perc_perm &
           Seed %in% list_params$Seed &
           Perm_meth %in% list_params$Perm_meth &
           Sele_meth %in% list_params$Sele_meth &
           Perm_hi_de %in% list_params$Perm_hi_de &
           Deg_hi %in% list_params$Deg_hi &
           background %in% list_params$background &
           break1 %in% list_params$break1
         ]
return(fdr)
}
fdr2<-select(fdr, list_params)

ggplot(data = fdr, aes(x = V1, y = V2, color = V3)) + 
  geom_line() + 
  ggtitle("Score Thresholds for each network size") + 
  xlab("Network size") + ylab("FDR threshold")+
  facet_wrap(~V4, scales = "free")+
  theme(legend.position = "bottom")


makeTable<-function(annotation_table, network_jaccard){
  table_new<-network_jaccard[,.SD[which.min(V3)], by=c("run", "filename")]
  table_new$id<-"min_V3"
  table<-rbind(table_new)
  table_new<-network_jaccard[,.SD[which.max(nrnodes)], by=c("run", "filename")]
  table_new$id<-"max_nodes"
  table<-rbind(table, table_new)
  table_new<-network_jaccard[,.SD[which.max(nrinteractions)], by=c("run", "filename")]
  table_new$id<-"max_inter"
  table<-rbind(table, table_new)
  table_new<-network_jaccard[,.SD[which.min(norm)], by=c("run", "filename")]
  table_new$id<-"min_norm"
  table<-rbind(table, table_new)
  table_new<-network_jaccard[,.SD[which.min(nrinteractions)], by=c("run", "filename")]
  table_new$id<-"min_inter"
  table<-rbind(table, table_new)
  table_new<-network_jaccard[,.SD[which.min(nrnodes)], by=c("run", "filename")]
  table_new$id<-"min_nodes"
  table<-rbind(table, table_new)
  return(table)
}


tab<-makeTable(anno, eval)
g1<-ggplot(tab[1:1000], aes(y=jaccard_dists, x=run, fill=id))+geom_boxplot()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), unit="cm"))+
  guides(fill=F)+scale_color_brewer(aesthetics = "fill")
#plot(g1)

min<-tab[,.SD[which.min(jaccard_dists)], by=c("run", "id")]
ggplot(min, aes(x=run, y=jaccard_dists))+geom_point()+facet_wrap(~id)

med<-min[,.SD[which.max(jaccard_dists)], by=c("run")]
ggplot(med, aes(x=run, y=jaccard_dists))+geom_point()
med[jaccard_dists>0.5]
