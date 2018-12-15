test<-abs(rnorm(100, mean=0.01, sd=0.01))
test<-data.table(test)
test$mod<-"foreground"
test2<-data.table(test=rnorm(1000, mean= 0.5, sd=0.15))
test2<-test2[test>0]
test2$mod<-"background"
test<-rbind(test, test2)


p<-ggplot(test, aes(x=test, y=..count.., group=mod, col=mod, fi))+
  geom_density()+scale_color_manual(values=c("#375A8CFF","#440154FF"))+
  theme(legend.title = element_blank(), plot.title = element_text(size=20),
        axis.text = element_text(size=15))+
  xlab("P-value")+ylab("Density")+
  ggtitle("Simplistic Foreground and Background Distributions")
ggsave(p, file="/home/anne/Masterarbeit/final_presentation/pictures/toydata_dist.pdf")
