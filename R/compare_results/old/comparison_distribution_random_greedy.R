require(data.table)
require(ggplot2)
require(viridis)
require(colorspace)
require(cowplot)


data<-fread("Masterarbeit/random/results/run_meanLog_greedy/pgen.txt", fill = T)
data2<-fread("Masterarbeit/random/results/run_meanLog_random/pgen.txt", fill = T)

scores<-data[,.(V1,V3)]
#scores$heuristic<-"greedy"
scores2<-data2[,.(V1,V3)]
#scores2$heuristic<-"random"
#scores<-rbind(scores, scores2)


col<-viridis(5)
col<-c(col, lighten(viridis(5), 0.4))
scores$V1<-as.factor(scores$V1)
scores2$V1<-as.factor(scores2$V1)
colnames(scores)<-c("Network_size", "Score")
colnames(scores2)<-c("Network_size", "Score")
greedy<-ggplot(scores[Network_size %in% c(10, 56,101, 201, 301)], aes(x=Score, group=Network_size))+
  geom_density(aes(color=Network_size))+
  scale_color_manual(values=col)+
  ggtitle("Comparison of the threshold distibutions random vs. greedy")+
  ylab("Density")+xlab("Score (greedy)")+theme(legend.position = "bottom")
plot(greedy)
random<-ggplot(scores2[Network_size %in% c(10, 56,101, 201, 301)], aes(x=Score, group=Network_size))+
  geom_density(aes(color=Network_size))+
  scale_color_manual(values=col)+
  ylab("Density")+xlab("Score (random)")+theme(legend.position = "bottom" )

g<-plot_grid(greedy, random, ncol = 1)
ggsave(g, file=file.path("Masterarbeit/random/comparison_distibution.pdf"), height = 20, width = 20, units = "cm")
plot(g)
sol<-fread("Masterarbeit/random/results/run_meanLog_greedy/file.stat")
siz<-fread("Masterarbeit/random/results/run_meanLog_greedy/file.nodes")

sol2<-fread("Masterarbeit/random/results/run_meanLog_random/file.stat")
siz2<-fread("Masterarbeit/random/results/run_meanLog_random/file.nodes")

r<-siz[, .N, by="V2"]
g<-siz2[, .N, by="V2"]

siz2[V2=="graph0"]
siz[V2=="graph0"]

n<-fread("/home/anne/Masterarbeit/final/sample_data/biogrid/expression_type_24_1_59_19.tsv")
nn<-n[V2=="differential"]
