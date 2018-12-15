require(data.table)
require(ggplot2)
require(viridis)
require(gridExtra)
require(colorspace)
require(cowplot)
source("~/Masterarbeit/R/compare_results/functions.R")

fdr<-fread("Masterarbeit/final2/StringDB//random//fdr_all.tsv")
annot<-fread("Masterarbeit/final2/StringDB//random//anno_all.tsv")

#random results
eval<-fread("Masterarbeit/final2/StringDB/random/eval_all.tsv")
jj<-makePlot(annotation_table =annot[(Aggr_meth=="meanLog"| window==15)], eval[run %in% annot[(Aggr_meth=="meanLog"| window==15)]$orig])
plot(jj)
ggsave(jj,file="Masterarbeit/final2/figures/random_jaccard_indices_StringDB.pdf", width = 30, height = 15, units = "cm")

eval[, percMean:=rowMeans(.SD, na.rm = T) , .SDcols = c("perc1", "perc2", "perc3")]
eval$percMean[eval$percMean>1]<-1
#selec<-eval[run %in% annot2[Perc_perm %in% c(50,60)]$orig]
perce<-makePlot2(annotation_table =annot, eval)
ggsave(perce,file="Masterarbeit/final2/figures/random_fraction_overlap_StringDB.pdf", width = 40, height = 30, units = "cm")



#jj<-makePlot(annotation_table =annot, eval)
# all networks extremly large
# a<-"Masterarbeit/final2/bio"
# min_size<-c()
# max_score<-c()
# f<-c()
# for(file in list.dirs(a)){
#   if(file.exists(file.path(file,"file.nodes"))){
# nodes<-fread(file.path( file,"file.nodes"))
# scores<-fread(file.path( file,"file.stat"))
# min_size<-c(min_size, min(nodes[,.N, by="V2"]$N))
# max_score<-c(max_score, max(scores$V3))
# f<-c(f,file)
# }
# 
# }
# dat<-data.table(f, min_size)
# dd<-dat[,.N, by="f"]

fdr$V1<-as.factor(fdr$V1)
sum<-fdr[V4%in% annot$orig[c(1,6,11)]]
sum$V4<-as.factor(sum$V4)
labels<-c("random-meanLog", "random-normSum","random-sum")
levels(sum$V4)<-labels
g1<-ggplot(sum, aes(x=V1,y=V2))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~V4, scales = "free_y" ,ncol = 1, )+
  scale_x_discrete(breaks=seq(0, 400, 20))+ylab("Threshold value")+xlab("Network Size")+
  theme(strip.text.x = element_text(size=20))
plot(g1)

fdr2<-fread("Masterarbeit/final2/StringDB/greedy///fdr_all.tsv")
annot2<-fread("Masterarbeit/final2/StringDB/greedy//anno_all.tsv")
fdr2$V1<-as.factor(fdr2$V1)
labels<-c("greedy-meanLog", "greedy-normSum", "greedy-sum")
sum2<-fdr2[V4%in% annot2$orig[c(4,8,145)]]
sum2$V4<-as.factor(sum2$V4)
levels(sum2$V4)<-labels
g2<-ggplot(sum2, aes(x=V1,y=V2))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~V4, scales = "free_y" ,ncol = 1)+
  scale_x_discrete(breaks=seq(0, 400, 20))+ylab("Threshold Value")+xlab("Network Size")+
  theme(strip.text.x = element_text(size=20))
#plot(g2)
p<-plot_grid(g1,g2, ncol = 2)

ggsave(p,file="Masterarbeit/final2/figures/compare_thresholds_stringDB.pdf", width = 40, height = 30, units = "cm")



tracker<-NULL
thresho<-NULL
d<-list.dirs("/home/anne/Masterarbeit/final2/StringDB/random/normDegSum/largest/String900_normDegSum_0.05_20_235_nodeswap_median_20_200_random_largest_15", recursive = F)
for(dir in d){
  tr<-fread(file=file.path(dir, "file.tracker"))
  tr$dir<-basename(dir)
  tracker<-rbind(tracker, tr)
  thresh<-fread(file=file.path(dir, "teststats.txt"))
  thresh$dir<-basename(dir)
  thresho<-rbind(thresho, thresh)
}

tracker$dir<-as.factor(tracker$dir)
thresho$dir<-as.factor(thresho$dir)
tracker$V2<-as.factor(tracker$V2)
thresho$V1<-as.factor(thresho$V1)
#ggplot(thresho, aes(x=as.numeric(V1), y=V2, fill=dir))+geom_line()
#ggplot(tracker[dir=="result_10_8"], aes(x=as.numeric(V2), y=V3, fill=V1))+geom_line()

mean<-thresho[,mean(V2), by=V1]
colnames(mean)<-c("index", "mean")
mean$type<-"threshold"
mean2<-tracker[,mean(V3), by=V2]
colnames(mean2)<-c("index", "mean")
mean2$type<-"tracker"
mean<-rbind(mean,mean2)

plo<-ggplot(mean, aes(x=index, y=mean, group=type))+geom_line()+
  ylab("Log10(Mean(Threhold))")+xlab("Network Size")+scale_y_log10()+
  ggtitle("Comparison FDR-Theshold - Growing Network Score\nAggregation=Degree Normalized Sum")+
  scale_x_discrete(breaks=seq(0, 400, 20))
ggsave(plo, file = "Masterarbeit/final2/figures/fdr_vs_actual_score_StringDB.pdf", width = 20, height = 15, units = "cm")

#the same for the actual solution
  track<-fread("Masterarbeit/real_data/results/run_final_15_005//file.tracker")
test<-fread("Masterarbeit/real_data/results/run_final_15_005/teststats.txt")

track<-fread("Masterarbeit/final2/biogrid/greedy/edgerewire/maximum/normDegSum/15/20/biogrid_normDegSum_0.05_50_235_edgerewire_median_20_200_greedy_maximum_distance_15/result_19_11/file.tracker")
test<-fread("Masterarbeit/final2/biogrid/greedy/edgerewire/maximum/normDegSum/15/20/biogrid_normDegSum_0.05_50_235_edgerewire_median_20_200_greedy_maximum_distance_15/result_19_11/teststats.txt")

track$V1<-as.factor(track$V1)
test$V1<-as.numeric(test$V1)
track$V2<-as.numeric(track$V2)
track$V3<-as.numeric(track$V3)
colnames(track)<-c("track", "size", "score")
test$t<-"threshold"
colnames(test)<-c("size", "score","track" )
test$size<-test$size+1
#test[size>15]$score<-test[size==15]$score
#test<-test[1:15]

track<-rbind(test,track)
track[,col:=ifelse(track=="threshold", "threshold", "graph")]
track[,w:=ifelse(track=="threshold", 0.1, 0.05)]
track<-track[order(col)]
pal<-c(brewer.pal(9, "Blues")[3], rep(brewer.pal(9, "Blues")[9],1))
thre<-ggplot(track, aes(x=size,y=score, group=track,color=col, size=w))+geom_line()+
  scale_color_manual(aesthetics = "color",  values = pal)+
  scale_x_continuous(breaks=seq(0, 400, 20))+theme(legend.title = element_blank(), legend.position = "bottom")+
  ggtitle("Comparison FDR-Theshold - Growing Network Score\n
          Permutation =Edge-rewire; Aggregation Method = NormDegSum")+
  xlab("Network Size")+ylab("Score")+scale_size(range=c(0.25, 0.5), guide=FALSE)
thre
#ggplot()+geom_line(data=test[V1<16], aes(x=V1, y=V2))+scale_y_log10()
ggsave(thre, file = "Masterarbeit/final2/figures/fdr_vs_actual_score_real_StringDB.pdf", width = 20, height = 15, units = "cm")



eval<-fread("Masterarbeit/final2/StringDB/greedy/eval_all.tsv")
a<-annot2[Perc_perm %in% c(50,60) & FDR==0.05]
selec<-eval[run %in% a$orig]
jj<-makePlot(annotation_table =annot2, eval)
plot(jj)
ggsave(jj,file="Masterarbeit/final2/figures/jaccard_indices_StringDB.pdf", width = 40, height = 30, units = "cm")

selec<-eval[run %in% a$orig]
jj<-makePlot(annotation_table =a, selec)
ggsave(jj,file="Masterarbeit/final2/figures/jaccard_indices_StringDB_reduced.pdf", width = 40, height = 30, units = "cm")


eval[, percMean:=rowMeans(.SD, na.rm = T) , .SDcols = c("perc1", "perc2", "perc3")]
eval$percMean[eval$percMean>1]<-1
perce<-makePlot2(annotation_table =annot2, eval)
plot(perce)
ggsave(perce,file="Masterarbeit/final2/figures/fraction_overlap_StringDB.pdf", width = 40, height = 30, units = "cm")
a<-annot2[Perc_perm %in% c(50,60) & FDR==0.05]
selec<-eval[run %in% a$orig]
p<-makePlot2(a,selec)
ggsave(p,file="Masterarbeit/final2/figures/fraction_overlap_StringDB_reduced.pdf", width = 40, height = 30, units = "cm")

