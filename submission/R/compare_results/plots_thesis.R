require(data.table)
#require(data)
source("~/Masterarbeit/R/compare_results/functions.R")

fdr<-fread("Masterarbeit/final2/biogrid/random/all/fdr_all.tsv")
annot<-fread("Masterarbeit/final2/biogrid/random/all/anno_all.tsv")
fdr$V1<-as.factor(fdr$V1)

sum<-fdr[V4%in% annot$orig[c(1,7,16)]]
sum$V4<-as.factor(sum$V4)
labels<-c("random-meanLog", "random-normSum", "random-sum")
levels(sum$V4)<-labels
g1<-ggplot(sum, aes(x=V1,y=V2))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~V4, scales = "free_y" ,ncol = 1)+
  scale_x_discrete(breaks=seq(0, 400, 20))+ylab("Threshold value")+xlab("Network Size")+
  theme(strip.text.x = element_text(size=20))
#plot(g1)

fdr2<-fread("Masterarbeit/final2/biogrid/greedy//fdr_all.tsv")
annot2<-fread("Masterarbeit/final2/biogrid/greedy//anno_all.tsv")
fdr2$V1<-as.factor(fdr2$V1)
labels<-c("greedy-meanLog", "greedy-normSum", "greedy-sum")
sum2<-fdr2[V4%in% annot2$orig[c(4,8,262)]]
sum2$V4<-as.factor(sum2$V4)
levels(sum2$V4)<-labels
g2<-ggplot(sum2, aes(x=V1,y=V2))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~V4, scales = "free_y" ,ncol = 1)+
  scale_x_discrete(breaks=seq(0, 400, 20))+ylab("Threshold Value")+xlab("Network Size")+
  theme(strip.text.x = element_text(size=20))
#plot(g2)
p<-plot_grid(g1,g2, ncol = 2)

ggsave(p,file="Masterarbeit/final2/figures/compare_thresholds.pdf", width = 40, height = 30, units = "cm")

annot2<-fread("Masterarbeit/final2/biogrid/greedy//anno_all.tsv")
eval<-fread("Masterarbeit/final2/biogrid/greedy/eval_all.tsv")
annot2[Perc_perm==60]
selec<-eval[run %in% annot2[Perc_perm %in% c(50,60)]$orig]
jj<-makePlot(annotation_table =annot2[Perc_perm %in% c(50,60)], selec)
plot(jj)
ggsave(jj,file="Masterarbeit/final2/figures/jaccard_indices.pdf", width = 40, height = 30, units = "cm")

eval[, percMean:=rowMeans(.SD, na.rm = T) , .SDcols = c("perc1", "perc2", "perc3")]
eval$percMean[eval$percMean>1]<-1
selec<-eval[run %in% annot2[Perc_perm %in% c(50,60)]$orig]
perce<-makePlot2(annotation_table =annot2[Perc_perm %in% c(50,60)], selec)
ggsave(perce,file="Masterarbeit/final2/figures/fraction_overlap.pdf", width = 40, height = 30, units = "cm")
