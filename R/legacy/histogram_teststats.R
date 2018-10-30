require(data.table)
require(ggplot2)
require(viridis)

teststats<-fread("/home/anne/Documents/Master/MA/Testing/out/dist/pdistributionT.txt", header= FALSE, sep = "\t")
teststats<-as.data.table(t(teststats))[2:ncol(teststats),]
teststats<-melt(teststats)
colnames(teststats)<-c("NetworkSize", "Pval")

gp<-ggplot(teststats[which(teststats$NetworkSize %in% c("V1", "V10", "V40"))], aes(x=Pval, color = NetworkSize))+
  geom_density(stat = "ecdf")+scale_colour_viridis(discrete = T)
p<-ggplot(teststats[which(teststats$NetworkSize %in% c("V1", "V10", "V40"))], aes(x=Pval, color = NetworkSize))+
  geom_density(adjust = 1/10)+scale_colour_viridis(discrete = T)

#+ theme(legend.position="none")
#plot(gp)
#plot(p)  

teststatsF<-fread("/home/anne/Documents/Master/MA/Testing/out/dist/pdistributionT.txt", header= FALSE, sep = "\t")
teststatsF<-as.data.table(t(teststatsF))[2:ncol(teststatsF),]
teststatsF<-melt(teststatsF)
colnames(teststatsF)<-c("NetworkSize", "Pval")

gpF<-ggplot(teststatsF[which(teststatsF$NetworkSize %in% c("V1", "V10", "V40"))], aes(x=Pval, color = NetworkSize))+
  geom_density(stat = "ecdf")+scale_colour_viridis(discrete = T)+ggtitle("includeBackgorund = F")
pF<-ggplot(teststatsF[which(teststatsF$NetworkSize %in% c("V1", "V10", "V40"))], aes(x=Pval, color = NetworkSize))+
  geom_density(adjust = 1/10)+scale_colour_viridis(discrete = T)+ggtitle("includeBackgorund = F")

#+ theme(legend.position="none")
#plot(gpF)
#plot(pF)  
gridExtra::grid.arrange(gp, p, gpF, pF)


teststats2<-fread("/home/anne/Documents/Master/MA/Testing/out/dist/distribution.txt", header= FALSE, sep = "\t")
teststats2<-as.data.table(t(teststats2))[2:ncol(teststats2),]
teststats2<-melt(teststats2)
colnames(teststats2)<-c("NetworkSize", "Pval")
p2<-ggplot(teststats2[which(teststats2$NetworkSize %in% c("V5")),], aes(x=Pval, color = NetworkSize))+
  geom_density(adjust = 1/10)+scale_colour_viridis(discrete = T)
plot(p2) 

chisq.test(na.omit(teststats2[which(teststats2$variable %in% c("V10")),]$value))

# revert file order first for easier reading
fi<-fread("/home/anne/Documents/Master/MA/Testing/out/dist/pvals2.txt", 
          sep = "\t", fill = T, header=F)
fi<-t(fi)
fi<-fi[!is.na(fi$value)]

ggplot(fi[fi$variable %in% c("V1", "V10", "V44")], aes(x=value))+geom_histogram()+facet_wrap("variable")
