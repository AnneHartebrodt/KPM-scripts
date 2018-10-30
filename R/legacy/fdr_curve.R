require(data.table)
require(ggplot2)

fdr<-fread("~/Masterarbeit/Testing/out_new/thresholds_1000.txt")
fdr2<-fread("~/Masterarbeit/Testing/out_new/thresholds.txt")
fdr3<-fread("~/Masterarbeit/Testing/out_new/thresholds.txt")
fdr3$V1<-1:nrow(fdr3)
fdr4<-fread("~/Masterarbeit/Testing/out_new/1531213101071thresholds.txt")
fdr4$V1<-1:nrow(fdr4)
  fdr4<-fread("~/Masterarbeit/Testing/out_new/1531213546244thresholds.txt")
fdr4$V1<-1:nrow(fdr4)

ggplot(data=fdr3, aes(x = V1, y=V2,color="10000 samples"))+geom_point()+geom_line()+
  stat_smooth()+geom_point(data=fdr2, aes(x = V1, y=V2, color="10000 samples"))+
  geom_line(data=fdr2, aes(x = V1, y=V2, color="10000 samples"))+
  stat_smooth(data=fdr2, aes(x = V1, y=V2, color="10000 samples"))

ggplot(data=fdr3, aes(x = V1, y=V2, color="100"))+geom_point()+geom_line()+
  stat_smooth()+geom_point(data = fdr4, aes(x = V1, y=V2, color="100_1"))+ geom_line(data=fdr4, aes(x = V1, y=V2, color="100_1"))+
  stat_smooth(data=fdr4, aes(x = V1, y=V2, color="100_1"))

write.table(fdr2, "~/Masterarbeit/Testing/out_new/thresholds10000.txt")
