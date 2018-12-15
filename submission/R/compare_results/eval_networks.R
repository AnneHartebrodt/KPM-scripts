require(data.table)
require(gProfileR)
require(topGO)
source("Masterarbeit/R/compare_results/go_enrich_functions.R")

deseq<-fread("Masterarbeit/differential_out/dataHD/DESEQ_final/results.tsv")
top<-deseq[padj<0.001]

final<-"Masterarbeit/real_data/results/run_final_15_005/"
nodes<-fread(file.path(final,"file.nodes"), header = F)
nodes[,.N, by=V2]

graph<-fread(file.path(final,"file.graph"), header = F)
stat<-fread(file.path(final,"file.stat"), header = F)
maap.file<-"~/Masterarbeit/data/identifier_mappings/mart_export_description.txt"
maap<-fread(maap.file)
map.file<-"~/Masterarbeit/data/identifier_mappings/mart_export.txt"
map<-fread(map.file)

best<-stat[which.max(nodes[,.N, by=V2]$N)]$V1

# Test whether any new gene came up in the analysis compared to DE
un<-unique(graph$V1)
for(gt in 1:(length(un)-1)){
g1<-unique(c(graph[V1==un[gt]]$V2, graph[V1==un[gt]]$V4))
g2<-unique(c(graph[V1==un[gt+1]]$V2, graph[V1==un[gt+1]]$V4))
print(c(length(intersect(g1,g2)), length(g1), length(g2)))
print(length(which(!(g1 %in% top$V1))))
}

nodiff<-nodes[!(V1 %in% top$V1)]
#nodiff<-nodiff[N>1]
nodiff<-merge(nodiff, maap, by.x="V1", by.y="Gene stable ID")
map<-unique(map[,c(1,3)])
nodiff<-nodiff[,c("Gene name", "Gene description")]
nodiff<-nodiff[,.N, by=c("Gene name", "Gene description")]
nodiff<-merge(nodiff, map, by.x="V1", by.y="Gene stable ID")
fwrite(nodiff, file=file.path("Masterarbeit/real_data/results/plots/", "nodiff_genes_pKPM.tsv"))


union_grapgh<-unique(graph[,.(V2,V4)])
#Test whether huntingtin came up in the analysis
nodes[V1=="ENSG00000197386"]

#test overlap with Kegg pathway
pathway_id<-"hsa05016"
kegg<-KEGGREST::keggGet(pathway_id)
kegglist<-unlist(kegg[[1]]$GENE)
ind1<-seq(1,386,2)
ind2<-seq(2,386,2)
kegglist<-data.table(gene_name=kegglist[ind1], desc=kegglist[ind2])
kegglist$gene_name<-as.integer(kegglist$gene_name)
kegglist<-merge(kegglist, map, by.x="gene_name", by.y="NCBI gene ID")

uni.nodes<-unique(nodes$V1)
g.best<-nodes[V2==best]
kegglist[`Gene stable ID` %in% g.best$V1]
kegg_ol<-kegglist[`Gene stable ID` %in% uni.nodes]
fwrite(kegg_ol, file=file.path("Masterarbeit/real_data/results/plots/kegg_overlap.tsv"), sep = "\t")

#gene set enrichment
select<-function (allScore) {
  return(allScore < 0.01)
  #return(allScore %in% bionet.nodes)
}
go<-go_list()
go_dir<-"/home/anne/Masterarbeit/evalutation_all_methods/go"
#best graph = largest graph
geneL<-deseq$padj
names(geneL)<-deseq$V1
best.nodes<-unique(c(graph[V1==best]$V2, graph[V1==best]$V4))
geneL[which(!( deseq$V1 %in% best.nodes))]<-1
sampleGOdata.best <- new("topGOdata", 
                    description = "Simple session", ontology = "BP",
                    allGenes = geneL, 
                    nodeSize = 10,
                    geneSelectionFun= select,
                    annotationFun=annFUN.gene2GO,
                    gene2GO = go)
result_pkpm.best<-runTest(sampleGOdata.best, algorithm = "parentchild", statistic = "fisher")
terms.best<-annot(result_pkpm.best)
fwrite(terms.best, file=file.path(go_dir, "pkpm_best.tsv"), col.names = F, sep="\t")

#union graph
geneL<-deseq$padj
names(geneL)<-deseq$V1
uni.nodes<-unique(c(union_grapgh$V2, union_grapgh$V4))
geneL[which(!( deseq$V1 %in% uni.nodes))]<-1
sampleGOdata <- new("topGOdata", 
                    description = "Simple session", ontology = "BP",
                    allGenes = geneL, 
                    nodeSize = 10,
                    geneSelectionFun= select,
                    annotationFun=annFUN.gene2GO,
                    gene2GO = go)
result_pkpm<-runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
terms<-annot(result_pkpm)
fwrite(terms, file=file.path(go_dir, "pkpm_union.tsv"), col.names = F, sep="\t")
i<-intersect(terms.best[score<0.01]$name, terms[score<0.01]$name)

#
select<-function (allScore) {
  return(allScore < 0.001)
  #return(allScore %in% bionet.nodes)
}

geneL<-deseq$padj
names(geneL)<-deseq$V1
sampleGOdata.diff <- new("topGOdata", 
                    description = "Simple session", ontology = "BP",
                    allGenes = geneL, 
                    nodeSize = 10,
                    geneSelectionFun= select,
                    annotationFun=annFUN.gene2GO,
                    gene2GO = go)
result_pkpm.diff<-runTest(sampleGOdata.diff, algorithm = "parentchild", statistic = "fisher")
terms.diff<-annot(result_pkpm.diff)
fwrite(terms.best, file=file.path(go_dir, "deseq_top.tsv"), col.names = F, sep="\t")


# nondiff does not work to few genes
# geneL<-deseq$padj
# names(geneL)<-deseq$V1
# geneL[which(!( deseq$V1 %in% nodiff$V1))]<-1
# sampleGOdata.nodiff <- new("topGOdata", 
#                     description = "Simple session", ontology = "BP",
#                     allGenes = geneL, 
#                     nodeSize = 5,
#                     geneSelectionFun= select,
#                     annotationFun=annFUN.gene2GO,
#                     gene2GO = go)
# result_pkpm.nodiff<-runTest(sampleGOdata.nodiff, algorithm = "parentchild", statistic = "fisher")
# terms.nodiff<-annot(result_pkpm.nodiff)
# 
# 
# top.no.mod<-deseq
# top.no.mod[V1 %in% uni.nodes]$padj<-1
# top.no.mod[!(V1 %in% top$V1)]$padj<-1