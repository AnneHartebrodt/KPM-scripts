require(data.table)
require(ggplot2)
require(cowplot)
require(VennDiagram)
source("Masterarbeit/R/compare_results/go_enrich_functions.R")
require(graph)
require(viridis)

final<-"Masterarbeit/real_data/results/run_final_15_005/"
bionet.file<-"~/Masterarbeit/bionet_out/dataHD/StringDB/module.sif"
bionet.700.file<-"~/Masterarbeit/bionet_out/Huntington700//module.sif"
probal.file<-file.path(final,"file.nodes")
probal.stat.file<-file.path(final, "file.stat")
web.file<-"~/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/5/Pathway-01-NODES-(1).txt"
diff.file<-"~/Masterarbeit/differential_out/dataHD/DESEQ_final/results.tsv"
map.file<-"~/Masterarbeit/data/identifier_mappings/mart_export.txt"
maap.file<-"~/Masterarbeit/data/identifier_mappings/mart_export_description.txt"
out.dir<-"/home/anne/Masterarbeit/evalutation_all_methods/"
go_dir<-file.path(out.dir, "go")

maap<-fread(maap.file)
bionet<-fread(bionet.file, header=F)
bionet.700<-fread(bionet.700.file, header=F)
probal<-fread(probal.file)
#probal.stat<-fread(probal.stat.file)
web<-fread(web.file)
diff<-fread(diff.file)
colnames(diff)<-c("geneName", colnames(diff)[2:7])
map<-fread(map.file)

# choosing largest as best
best_probal<-probal[,.N, by=V2][which.max(probal[,.N, by=V2]$N)]$V2

probal.1<-probal[V2==best_probal]
probal.1<-data.table(nodes=probal.1$V1)
bionet.nodes<-data.table(nodes=unique(c(bionet$V1, bionet$V3)))
bionet.700.nodes<-data.table(nodes=unique(c(bionet.700$V1, bionet.700$V3)))
diff.nodes<-data.table(nodes=diff[padj<0.001]$geneName)


KPM_exception_nodes<-list()
i=5
#for (i in c(0,2,4,6,8,10)){
web.nodes<-web#####create the venn diagram
colnames(web.nodes)<-c("nodes")  
area1<-nrow(probal.1)
  area2<-nrow(bionet.nodes)
  area3<-nrow(diff.nodes)
  area4<-nrow(web.nodes)
  area5<-nrow(bionet.700.nodes)
  
  n12<-length(intersect(probal.1$nodes, bionet.nodes$nodes))
  n13<-length(intersect(probal.1$nodes, diff.nodes$nodes))
  n14<-length(intersect(probal.1$nodes, web.nodes$nodes))
  n15<-length(intersect(probal.1$nodes, bionet.700.nodes$nodes))
  
  n23<-length(intersect(bionet.nodes$nodes, diff.nodes$nodes))
  n24<-length(intersect(bionet.nodes$nodes, web.nodes$nodes))
  n25<-length(intersect(bionet.nodes$nodes, bionet.700.nodes$nodes))
  
  n34<-length(intersect(diff.nodes$nodes, web.nodes$nodes))
  n35<-length(intersect(diff.nodes$nodes, bionet.700.nodes$nodes))
  
  n45<-length(intersect(web.nodes$nodes, bionet.700.nodes$nodes))
  
  n123<-length(intersect(intersect(probal.1$nodes, bionet.nodes$nodes), diff.nodes$nodes))
  n124<-length(intersect(intersect(probal.1$nodes, bionet.nodes$nodes), web.nodes$nodes))
  n125<-length(intersect(intersect(probal.1$nodes, bionet.nodes$nodes), bionet.700.nodes$nodes))
  n134<-length(intersect(intersect(probal.1$nodes, diff.nodes$nodes), web.nodes$nodes))
  n135<-length(intersect(intersect(probal.1$nodes, diff.nodes$nodes), bionet.700.nodes$nodes))
  n145<-length(intersect(intersect(probal.1$nodes, web.nodes$nodes), bionet.700.nodes$nodes))
  
  n235<-length(intersect(intersect(bionet.nodes$nodes, diff.nodes$nodes), bionet.700.nodes$nodes))
  n234<-length(intersect(intersect(bionet.nodes$nodes, diff.nodes$nodes), web.nodes$nodes))
  n245<-length(intersect(intersect(bionet.nodes$nodes, web.nodes$nodes ), bionet.700.nodes$nodes))
  
  n345<-length(intersect(intersect(diff.nodes$nodes, web.nodes$nodes),bionet.700.nodes$nodes))
  
  n1234<-length(intersect(intersect(intersect(bionet.nodes$nodes, web.nodes$nodes), diff.nodes$nodes), probal.1$nodes))
  n1235<-length(intersect(intersect(intersect(bionet.nodes$nodes, bionet.700.nodes$nodes), diff.nodes$nodes), probal.1$nodes))
  n1245<-length(intersect(intersect(intersect(probal.1$nodes, web.nodes$nodes), bionet.700.nodes$nodes), bionet.nodes$nodes))
  n1345<-length(intersect(intersect(intersect(probal.1$nodes, web.nodes$nodes), bionet.700.nodes$nodes), diff.nodes$nodes))
  
  n2345<-length(intersect(intersect(intersect(bionet.nodes$nodes, web.nodes$nodes), bionet.700.nodes$nodes), diff.nodes$nodes))
  
  n12345<-length(intersect(intersect(intersect(intersect(bionet.nodes$nodes, web.nodes$nodes), bionet.700.nodes$nodes), diff.nodes$nodes), probal.1$nodes))
  # pdf(file.path("~/Masterarbeit/evalutation_all_methods/", paste0("vennKPM0001",i,".pdf")))
  # 
  # ggplot()+geom_blank()+ggtitle("Overlap between the gene sets in the best solutions")
  # draw.quintuple.venn(area1, area2, area3, area4, area5, n12, n13, n14, n15,
  #                n23, n24, n25, n34, n35, n45, n123, n124, n125, n134,
  #                n135, n145, n234, n235, n245, n345, n1234, n1235,
  #                n1245, n1345, n2345, n12345, 
  #                category = c("pKPM", "Bionet 900", "DE Genes", paste0("KPM\nk=", i), "Bionet 700"),
  #                at.pos = c(0, 287.5, 215, 145, 70), cat.dist =rep(0.2, 5), 
  #                fill = c("#eff3ff","#bdd7e7", "#6baed6","#3182bd","#08519c")
  #                )
  # dev.off()
  pdf(file.path("~/Masterarbeit/evalutation_all_methods/", paste0("vennKPM10.pdf")))
  ggplot()+geom_blank()+ggtitle("Overlap between the gene sets in the best solutions")
  draw.quad.venn(area1, area2, area3, area4, n12, n13, n14,
                      n23, n24,  n34, n123, n124, n134,  n234, n1234,
                      category = c("pKPM", "Bionet 900", "DE Genes", paste0("KPM\nk=5")),
                      cat.pos = c(350, 10,0,0), cat.dist =c(-0.25, -0.24,-0.15,-0.15) , 
                      fill = c("#eff3ff","#bdd7e7","#6baed6","#3182bd" ), rotation.degree = 180
  )
  dev.off()
  KPM_exception_nodes[[((i/2)+1)]]<-setdiff( web.nodes$nodes, diff.nodes$nodes)

pdf(file.path("~/Masterarbeit/evalutation_all_methods/", "bionet_overlap.pdf"))
ggplot()+geom_blank()+ggtitle("Overlap between the Bionet with String700 and String900")
draw.pairwise.venn(length(bionet.700.nodes$nodes), length(bionet.nodes$nodes), 
                   length(intersect(bionet.700.nodes$nodes,bionet.nodes$nodes )),
                   fill = c("#6baed6","#3182bd"), category = c("Bionet 700", "Bionet 900"))
dev.off()


web.only<-maap[`Gene stable ID` %in% unique(setdiff(web$value, diff[padj<0.01]$geneName))]
network<-"Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS900.tsv"
String900<-fread(network, header=F)
sub<-String900
sub$V1<-as.character(sub$V1)
sub$V3<-as.character(sub$V3)
g<-graphNEL()
g<-addNode(node = unique(as.character(c(sub$V1, sub$V3))), g)
g<-addEdge(as.vector(sub$V1), as.vector(sub$V3), g, 1)
#analyse the bionet unique genes
bionet.uniq<-setdiff(bionet.nodes$nodes, diff.nodes$nodes)
bionet.uniq.p<-diff[geneName %in% bionet.uniq]
bionet.outlier<-bionet.uniq.p[padj>0.005]
names<-fread("Masterarbeit/data/identifier_mappings/mart_export_names.txt")
g1<-ggplot(bionet.uniq.p, aes(x=padj))+geom_histogram(bins=50)+xlim(c(0,0.005))+
ggtitle("'High Scoring' Nodes in the Bionet Module")+xlab("P-Value")+ylab("#P-values")
ggsave(g1, file=file.path(out.dir, "histogram_only_bionet.pdf"), width = 20, height = 15, units = "cm")

bionet.hist<-degree(g)[bionet.uniq]
bionet.hist.top<-degree(g)[top[V1 %in% bionet.nodes$nodes]$V1]
wilcox.test.1<- wilcox.test(bionet.hist, bionet.hist.top)

hist.bio.1<-ggplot(data.table(bionet.hist), aes(x=bionet.hist))+geom_histogram()+
  ggtitle("Node Degree Distribution\np-value > 0.001")+xlab("Node Degree")+ylab("#Nodes")
hist.bio.2<-ggplot(data.table(bionet.hist.top), aes(x=bionet.hist.top))+geom_histogram()+
  ggtitle("Node Degree Distribution\np-value < 0.001")+xlab("Node Degree")+ylab("#Nodes")
pg<-plot_grid(hist.bio.2, hist.bio.1)
ggsave(pg, file=file.path(out.dir, "histogram_node_degrees.pdf"), width = 20, height = 15, units = "cm")

bionet.outlier<-merge(bionet.outlier, names, by.y = "Gene stable ID", by.x = "geneName")
bionet.outlier$degree<-degree(g)[bionet.outlier$geneName]
write.table(bionet.outlier, file = file.path(out.dir, "bionet_high_p_vals.tsv"))
bionet.uniq.p<-bionet.uniq.p[order(padj)]


#GeneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))


deg.web<-degree(g)[web.nodes$nodes]
deg.bionet<-degree(g)[bionet.nodes$nodes]
deg.probal<-degree(g)[probal.1$nodes]
deg.diff<-degree(g)[diff.nodes$nodes]
union.nodes<-unique(probal$V1)
deg.union<-degree(g)[union.nodes]

rel.union<-data.table(table(deg.union)/length(deg.union))
rel.union$module<-"pKPM Union"
colnames(rel.union)<-c("degree", "density", "module")
rel.bionet<-data.table(table(deg.bionet)/length(deg.bionet))
rel.bionet$module<-"Bionet"
colnames(rel.bionet)<-c("degree", "density", "module")
rel.probal<-data.table(table(deg.probal)/length(deg.probal))
rel.probal$module<-"pKPM"
colnames(rel.probal)<-c("degree", "density", "module")
rel.web<-data.table(table(deg.web)/length(deg.web))
rel.web$module<-"KPM K=10"
colnames(rel.web)<-c("degree", "density", "module")
rel.diff<-data.table(table(deg.diff)/length(deg.diff))
rel.diff$module<-"DE Genes"
colnames(rel.diff)<-c("degree", "density", "module")
rel.all<-rbind(rel.union, rel.bionet, rel.probal, rel.web, rel.diff)
rel.all$degree<-as.numeric(rel.all$degree)
colnames(rel.all)<-c("degree", "density", "Module")

g.comp<-ggplot(rel.all, aes(x=degree, y=density, col=Module))+geom_line()+
  scale_color_discrete(aesthetics = "module")+scale_x_log10()+
  xlab("Node degree")+ylab("#Nodes/Module Size")+
  scale_fill_viridis_d(aesthetics ="color")+
  theme(legend.position = "bottom")+
  ggtitle("Comparison of Node Degree Densities")

g.web<-ggplot(data.frame(n=deg.web), aes(x=n))+geom_histogram(bins=50)+ggtitle("Classical KPM")+
  ylab("#Nodes")+xlab("Node Degree in String900")
g.bionet<-ggplot(data.frame(n=deg.bionet), aes(x=n))+geom_histogram(bins=100)+ggtitle("Bionet")+
  ylab("#Nodes")+xlab("Node Degree in String900")
g.diff<-ggplot(data.frame(n=deg.diff), aes(x=n))+geom_histogram(bins=50)+ggtitle("DE Genes")+
  ylab("#Nodes")+xlab("Node Degree in String900")
g.probal<-ggplot(data.frame(n=deg.probal), aes(x=n))+geom_histogram(bins=50)+ggtitle("probabilistic KPM")+
  ylab("#Nodes")+xlab("Node Degree in String900")
g.union<-ggplot(data.frame(n=deg.union), aes(x=n))+geom_histogram(bins=50)+ggtitle("pKPM Union Node Set")+
  ylab("#Nodes")+xlab("Node Degree in String900")
plot(g.web)
plot(g.bionet)
plot(g.diff)
plot(g.probal)
plot(g.union)

pp<-plot_grid(g.probal,  g.bionet, g.diff, g.union, g.web, g.comp, ncol=3)
ggsave(pp,file=file.path("~/Masterarbeit/evalutation_all_methods/", "degree_dist_modules.pdf"), width=20, height = 20)


"ENSG00000197386" %in% bionet.nodes$nodes
"ENSG00000197386" %in% diff[padj<0.01]$geneName
"ENSG00000197386" %in% web.nodes$nodes
"ENSG00000197386" %in% probal.1$nodes
"ENSG00000197386" %in% union.nodes

adj.web<-unique(unlist(sapply(web.nodes$nodes, function(x) adj(g, x))))
adj.bionet<-unique(unlist(sapply(bionet.nodes$nodes, function(x) adj(g, x))))
adj.union<-unique(unlist(sapply(union.nodes, function(x) adj(g, x))))
adj.probal<-unique(unlist(sapply(probal.1$nodes, function(x) adj(g, x))))
#adj.web<-unique(unlist(sapply(web.nodes$nodes, function(x) adj(g, x))))
"ENSG00000197386" %in% adj.bionet
"ENSG00000197386" %in% adj.probal
"ENSG00000197386" %in% adj.union 
"ENSG00000197386" %in% adj.web


summary<-data.table(method = c("pKPM", "pKPM (union)", "Bionet", "classical KPM"),
           size = c(nrow(probal.1), length(union.nodes), nrow(bionet.nodes), nrow(web.nodes)),
           neighbourhood = c(length(adj.probal), length(adj.union), length(adj.bionet), length(adj.web)),
           containsHTT = c("ENSG00000197386" %in% probal.1$nodes,
                           "ENSG00000197386" %in% union.nodes,
                           "ENSG00000197386" %in% bionet.nodes$nodes,
                           "ENSG00000197386" %in% web.nodes$nodes),
           neighbourHTT = c("ENSG00000197386" %in% adj.probal,
                            "ENSG00000197386" %in% adj.union ,
                            "ENSG00000197386" %in% adj.bionet,
                            "ENSG00000197386" %in% adj.web)
                           )
fwrite(summary, "~/Masterarbeit/evalutation_all_methods/module_sizes.tsv", sep = "\t")

#GO enrichment
library(topGO)
select<-function (allScore) {
  return(allScore < 0.1)
  #return(allScore %in% bionet.nodes)
}

golist<-go_list()

#bionet
geneL<-diff$padj
names(geneL)<-diff$geneName
geneL[which(!( diff$geneName %in% bionet.nodes$nodes))]<-1
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneL,
                    nodeSize = 10,
                    geneSelectionFun= select,
                    annotationFun=annFUN.gene2GO,
                    gene2GO = golist)
result_bionet_only<-runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
terms_bionet_only<-annot(result_bionet_only)
fwrite(terms_bionet_only, file=file.path(go_dir, "bionet.tsv"), col.names = F, sep="\t")

#web-kpm
geneL<-diff$padj
names(geneL)<-diff$geneName
geneL[which(!( diff$geneName %in% web.nodes$nodes))]<-1
sampleGOdata.web <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneL,
                    nodeSize = 10,
                    geneSelectionFun= select,
                    annotationFun=annFUN.gene2GO,
                    gene2GO = golist)
result.web<-runTest(sampleGOdata.web, algorithm = "parentchild", statistic = "fisher")
terms.web<-annot(result.web)
fwrite(terms.web, file=file.path(go_dir, "KPM_k5.tsv"), col.names = F, sep="\t")

# probal only
geneL<-diff$padj
names(geneL)<-diff$geneName
geneL[which(!( diff$geneName %in% probal.1$nodes))]<-1
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneL,
                    nodeSize = 10,
                    geneSelectionFun= select,
                    annotationFun=annFUN.gene2GO,
                    gene2GO = golist)
result_pkpm<-runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
terms_pkpm<-annot(result_pkpm)
