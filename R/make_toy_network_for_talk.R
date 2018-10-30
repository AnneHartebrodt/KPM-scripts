require(data.table)

dat<-fread("Masterarbeit/talk/test.csv")

un<-unique(c(dat$V1, dat$V4))

un<-data.table(node=un)
dat$pp<-"pp"
dat<-dat[,c(1,5,4)]


write.table(file ="Masterarbeit/talk/test_CLEAN.csv", dat, sep="\t", quote = F, row.names = F)

#5527





un$p<-c(0.004, 0.01,0.003,0.003, 0.55,0.9,0.005, 0.08,0.03,0.005, 0.01,0.03,0.001, 0.01,0.1)
write.table(file ="Masterarbeit/talk/test_VALUES.csv", un, sep="\t", quote = F, row.names = F)
