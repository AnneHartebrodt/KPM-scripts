anno<-fread("Masterarbeit/final2/biogrid/random/all/anno_all.tsv")

anno$Perc_perm<-10
anno$Perm_meth<-"nodeswap"
anno$Perm_hi_de<-50
anno$Deg_hi<-200

fwrite(anno,"Masterarbeit/final2/biogrid/random/all/anno_all_clean.tsv", sep = "\t")
