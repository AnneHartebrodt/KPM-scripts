require(data.table)
require(ggplot2)
require(viridisLite)
require(heatmap.2)
require(gplots)

testNodes<-fread("Masterarbeit/Test_pipeline/sample_networks/StringDB/21.graph")
testNodes<-unique(c(testNodes$V2, testNodes$V4))

sig<-c(5914,23054,2033,5931,79685,85509, 3553)
gene_exp<-abs(rnorm(sig, 0 , 0.005))

not_sig<-testNodes[!testNodes %in% sig]
gene_exp<-c(gene_exp, runif(not_sig))

test_datea<-data.frame(gene_name=c(sig, not_sig), exp_val=gene_exp)

write.table(test_datea, "Masterarbeit/poster/node_values.txt", row.names = F, col.names = F, sep="\t")
gp<-ggplot(test_datea, aes(x=gene_exp))+geom_histogram(bins=30)
plot(gp)


heat<-fread("Masterarbeit/Test_pipeline/sample_data/StringDB/samplewise_49_61_34.tsv")

heat_red<-cbind(head(heat), tail(heat[,2:ncol(heat)])/15)
heat_red<-rbind(heat_red, cbind(tail(heat,20)[1:10]/15, head(heat, 20)[11:20, 2:ncol(heat)]))
heat_red<-as.data.table(heat_red)
colnames(heat_red)<-sapply(seq(1,ncol(heat_red)), function(x) paste0("P", x))
rownames(heat_red)<-sapply(seq(1,nrow(heat_red)), function(x) paste0("G", x))
pal <-RColorBrewer::brewer.pal(9, "Blues")
pal<-colorRampPalette(pal)(200)
heatmap.2(as.matrix(heat_red[,2:ncol(heat_red)]), col = pal, dendrogram = "none", trace = "none",
          key = F, labRow = rownames(heat_red))

heatmap.2(as.matrix(heat_red[,2:ncol(heat_red)]), col = pal, dendrogram = "none", trace = "none",
          key = F, labRow = F, labCol = F)
