---
title: "Graph overlap"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r, echo = F}
require(data.table)
require(ggplot2)
```

```{r}
gold<-fread("/home/anne/Documents/Master/MA/data/Kegg_Huntington_EntrezID.csv")
gold<-data.table(gene= unique(unlist(gold)))
write.table(gold, file = "/home/anne/Documents/Master/MA/data/Kegg_Huntington_Unique.csv", row.names = F)

  f1 <-fread("/home/anne/Documents/Master/MA/Testing/out_new/graph_out/1531384983334file.graph")
f1[,.N, by = V1]

jaccard<-function(g){
  sub1<-f1[V1==g,2:4] 
  sub1<-sub1[order(sub1$V2)]
  su<-unique(c(sub1$V2, sub1$V4))
  intersect_1<-length(which(su %in% gold$gene))
  union_1<-length(unique(c(su, gold$gene)))
  jaccard<-intersect_1/union_1
  return(jaccard)
}

jaccard_dists<-sapply(unique(f1$V1), function(g) jaccard(g))
which.max(jaccard_dists)
boxplot(jaccard_dists)

pval<-fread("/home/anne/Documents/Master/MA/Testing/out_new/graph_out/1531384983334file.stat")
pval[V2<0.05]


nodes<-fread("/home/anne/Documents/Master/MA/Testing/out_new/graph_out/1531384983334file.nodes")
nodes[V2=="graph165"]

ggplot(nodes[V2 %in% pval[V2<0.05]$V1], aes(x=V3))+geom_histogram()+facet_wrap("V2")
ggplot(nodes[!(V2 %in% pval[V2<0.05]$V1)], aes(x=V3))+geom_histogram()+facet_wrap("V2")

pchisq(475.7314, df = 708)
```


