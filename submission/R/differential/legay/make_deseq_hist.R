require(data.table)
require(ggplot2)
require(cowplot)


data<-fread("Masterarbeit/data/DESeq_out/DESeq_out.tsv")
g<-ggplot(data, aes(pvalue ))+geom_histogram(aes(col=I("#0080FF"),fill=I("#0080FF")), bins = 100)+
  ggtitle("P-value Distribution of a Real DE Analysis")+
  xlab("P-Value")+ylab("Count")#+geom_hline(yintercept = 190, color="darkblue")
plot(g)
ggsave(plot = g,filename = "Masterarbeit/thesis/figures/DEpvaluehist.pdf")

g2<-ggplot(data, aes(padj ))+geom_histogram(aes(col=I("#0080FF"),fill=I("#0080FF")), bins = 100)+
  ggtitle("P-value Distribution of a Real DE Analysis (Adjusted P-values)")+
  xlab("P-Value")+ylab("Count")
plot(g2)
ggsave(plot = g2,filename = "Masterarbeit/thesis/figures/DEAdjpvaluehist.pdf")
