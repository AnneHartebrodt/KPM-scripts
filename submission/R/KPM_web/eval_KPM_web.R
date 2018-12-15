require(data.table)
require(ggplot2)
require(viridis)
require(VennDiagram)
require(cowplot)

kpm.all<-NULL

dir_list<-c("~/Masterarbeit/Greedy_KPM_test/0.001/")
file_list<-c("/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/20/Pathway-01-INTERACTIONS-.txt",
             "/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/40/Pathway-01-INTERACTIONS-(2).txt",
             "/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/10/Pathway-01-INTERACTIONS-(1).txt",
             "/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/5/Pathway-01-INTERACTIONS-(1).txt",
             "/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/2/Pathway-01-INTERACTIONS-.txt",
             "/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/1/Pathway-01-INTERACTIONS-.txt",
             "/home/anne/Masterarbeit/Greedy_KPM_test/KPM-4.0/results/0/Pathway-01-INTERACTIONS-.txt")

dir<-c(20,40,10,5,2,1,0)
for(file in 1:length(file_list)){
  kpm.run<-fread(file_list[file], header = F)
  kpm.run$file<-file_list[file]
  kpm.run$dir<-dir[file]
  if(is.null(kpm.all)){
    kpm.all<-kpm.run
  }
  else{
    kpm.all<-rbind(kpm.run, kpm.all)
  }
}
#}
#colnames(K40)<-colnames(kpm.all)
#kpm.all<-rbind(kpm.all, K40)
sizes<-sapply(unique(kpm.all$file), function(x) length(unique(c(kpm.all[file==x]$V1, kpm.all[file==x]$V2))))
#sizes<-apply(unique(kpm.all[,.(file,dir)]),1,function(x) length(unique(c(kpm.all[file==x[1] & dir==x[2]]$V1, kpm.all[file==x[1] & dir==x[2]]$V3))))

dt<-unique(kpm.all[,.(file, dir)])
dt$sizes<-sizes

pathway_id<-"hsa05016"
kegg<-KEGGREST::keggGet(pathway_id)
kegglist<-unlist(kegg[[1]]$GENE)
ind1<-seq(1,386,2)
ind2<-seq(2,386,2)
kegglist<-data.table(gene_name=kegglist[ind1], desc=kegglist[ind2])


map.file<-"Masterarbeit/data/identifier_mappings/mart_export.txt"
map<-fread(map.file)
map<-unique(map[,.(`Gene stable ID`, `NCBI gene ID`)])
map<-map[!is.na(`NCBI gene ID`)]

keggi<-function(x){
  gene_list<-unique(c(kpm.all[file==x[1] & dir==x[2]]$V1, kpm.all[file==x[1] & dir==x[2]]$V3))
  gene_NCBI<-map[`Gene stable ID` %in% gene_list]$`NCBI gene ID`
  genes<-gene_NCBI[which( gene_NCBI %in% kegglist$gene_name)]
  return(length(gene_NCBI[which( gene_NCBI %in% kegglist$gene_name)]))
}
keggi2<-function(x){
  gene_list<-unique(c(kpm.all[file==x[1] & dir==x[2]]$V1, kpm.all[file==x[1] & dir==x[2]]$V3))
  gene_NCBI<-map[`Gene stable ID` %in% gene_list]$`NCBI gene ID`
  genes<-gene_NCBI[which( gene_NCBI %in% kegglist$gene_name)]
  return(genes)
}

sizes<-apply(unique(kpm.all[,.(file)]),1,function(x) keggi(x))
genes<-unique(unlist(apply(unique(kpm.all[,.(file, dir)]),1,function(x) keggi2(x))))

ol<-kegglist[gene_name %in% genes]
fwrite(ol,"~/Masterarbeit/Greedy_KPM_test/ol_KEGG.tsv")

dt$overlapKEGG<-sizes

dt<-dt[file!="KPM0_2.sif"]
#dt[,k:=as.factor(as.numeric(gsub(".sif", "",gsub("KPM", "", file))))]
dt[, dir:=as.factor(gsub("/", "",gsub("~/Masterarbeit/Greedy_KPM_test/", "", dir)))]



dt$dir<-as.integer(as.character(dt$dir))
res_sizes<-ggplot(dt, aes(x=dir, y=sizes))+
  geom_point(aes())+geom_line()+
  scale_fill_manual(aesthetics = "color", values = c("#6baed6", "#08519c"))+
  xlab("k")+ylab("Size of Best Module")
plot(res_sizes)
ggsave(res_sizes, file=file.path("~/Masterarbeit/Greedy_KPM_test/sizes.pdf"))

fwrite(dt,"~/Masterarbeit/Greedy_KPM_test/summary.tsv")

m<-melt(kpm.all[, c(1,3:5)], id.vars = c("file", "dir"))
fwrite(m, "~/Masterarbeit/Greedy_KPM_test/nodes.tsv" )

area1<-length(unique(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2)))
area2<-length(unique(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2)))
area3<-length(unique(c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)))
area4<-length(unique(c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#area5<-nrow(bionet.700.nodes)

n12<-length(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2)))
n13<-length(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)))
n14<-length(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#n15<-length(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), bionet.700.nodes$nodes))

n23<-length(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)))
n24<-length(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#n25<-length(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), bionet.700.nodes$nodes))

n34<-length(intersect(c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#n35<-length(intersect(c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2), bionet.700.nodes$nodes))

#n45<-length(intersect(c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2), bionet.700.nodes$nodes))

n123<-length(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2)), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)))
n124<-length(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2)), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#n125<-length(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2)), bionet.700.nodes$nodes))
n134<-length(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#n135<-length(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), bionet.700.nodes$nodes))
#n145<-length(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)), bionet.700.nodes$nodes))

#n235<-length(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), bionet.700.nodes$nodes))
n234<-length(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)))
#n245<-length(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2) ), bionet.700.nodes$nodes))

#n345<-length(intersect(intersect(c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)),bionet.700.nodes$nodes))

n1234<-length(intersect(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2)))
#n1235<-length(intersect(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), bionet.700.nodes$nodes), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2)))
#n1245<-length(intersect(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)), bionet.700.nodes$nodes), c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2)))
#n1345<-length(intersect(intersect(intersect(c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)), bionet.700.nodes$nodes), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)))

#n2345<-length(intersect(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)), bionet.700.nodes$nodes), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)))

#n12345<-length(intersect(intersect(intersect(intersect(c(kpm.all[dir==10]$V1,kpm.all[dir==10]$V2), c(kpm.all[dir==2]$V1,kpm.all[dir==2]$V2)), bionet.700.nodes$nodes), c(kpm.all[dir==0]$V1,kpm.all[dir==0]$V2)), c(kpm.all[dir==5]$V1,kpm.all[dir==5]$V2)))

pdf(file.path("~/Masterarbeit/Greedy_KPM_test/", "KPM_k_overlap.pdf"))
ggplot()+geom_blank()+ggtitle("Overlap of KPM results with different k")
draw.quad.venn(area1, area2, area3, area4, n12, n13, n14,
               n23, n24,  n34, n123, n124, n134,  n234, n1234,
               category = c("k=5", "k=10", "k=0", "k=2"),
               cat.pos = c(350, 10,0,0), cat.dist =c(-0.25, -0.24,-0.15,-0.15) , 
               fill = c("#eff3ff","#bdd7e7","#6baed6","#3182bd" ), rotation.degree = 180
)
dev.off()
