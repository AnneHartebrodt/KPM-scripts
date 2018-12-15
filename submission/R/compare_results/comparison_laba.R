require(data.table)

laba<-fread("Masterarbeit/Labadorf/labadorf_all_genes.csv")
res<-fread("Masterarbeit/differential_out/dataHD/DESEQ_new/results.tsv")
res<-res[-which(is.na(padj))]
res<-res[order(padj)]
res$rank<-1:nrow(res)

laba$ENSG<-gsub("\\.[0-9]*", "", laba$ENSG)
laba<-laba[order(padj)]
laba$rank<-1:nrow(laba)

r<-merge(laba, res, by.x="ENSG", by.y="V1")
r[, rankproduct:=rank.x*rank.y]
r<-r[order(rankproduct)]
d<-r[,diff:=abs(rank.x-rank.y)][rank.x<100 | rank.y<100]
d[, .(rank.x, rank.y, rankproduct, diff)]
