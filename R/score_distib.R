require(data.table)
require(ggplot2)
require(cowplot)


teststats<-fread("/home/anne/Masterarbeit/Test_pipeline/results/String_.01/StringDB_normSum/result_19_61_97/pgen.txt", fill = T)
teststats<-as.data.table(teststats)
su<-rowSums(teststats, na.rm = T)-teststats$V1
dat<-data.table(Nsize=teststats$V1, rs=su)
sub<-dat[Nsize==20]$rs
sub<-sort(sub)
val<-sub[floor(length(sub)*0.05)]

ggplot(dat[Nsize == 20], aes(x=rs))+
  geom_histogram(aes(y=..density..), bins = 50, fill="lightblue")+
  geom_density(color="midnightblue")+
  geom_vline(xintercept =val, color="darkblue")+
  xlab("Network Score")+
  ylab("Density")+
  ggtitle("Distribution of Scores for Randomized Networks of Size 20")

ggsave("/home/anne/Masterarbeit/Test_pipeline/results/StringDB/result_10_85_73.tsv/score_distrib20.pdf")


