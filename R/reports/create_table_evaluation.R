require(data.table)
require(ggplot2)
require(gridExtra)
require(viridis)
require(cowplot)
require(stringi)
require(RColorBrewer)
require(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="batch evaluation folder", metavar="character"),
  make_option(c("-n", "--networks"), type="character", default=NULL, 
              help="goldstandard network folder", metavar="character"),
  make_option(c("-g", "--gold"), type="character", default=NULL, 
              help="goldstandard expression data folder", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output folder name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

filelocation<-opt$file
goldstandard<-opt$gold
outdir<-opt$out
direct<-opt$networks


#filelocation<-"~/Masterarbeit/final2/biogrid/random/all/sum/max_lag"
#goldstandard<-"/home/anne/Masterarbeit/test_runs/sample_data/String900/"
#direct<-"/home/anne/Masterarbeit/test_runs/sample_networks/String900/"
#outdir<-"/home/anne/Masterarbeit/test_runs/results/String900/"

#################
# FUNCTIONS
#################

# jaccard distance
jaccard <- function(g, f1) {
  sub1 <- f1[V1 == g, 2:4]
  sub1 <- sub1[order(sub1$V2)]
  su <- unique(c(sub1$V2, sub1$V4))
  intersect_1 <- length(which(su %in% gold$gene))
  union_1 <- length(unique(c(su, gold$gene)))
  jaccard <- intersect_1 / union_1
  return(list(jaccard, union_1, intersect_1, g))
}
extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}

fdr <- data.frame(V1 = c(), V2 = c(), V3 = c(), V4=c())
dd<-dir(filelocation)
#dd<-grep("2453", dd, value = T)
for(pdir in dd){
  for (d in dir(file.path(filelocation,pdir))) {
    if(file.exists(file.path(filelocation, pdir, d, "teststats.txt"))){
      fdr_temp <- fread(file.path(filelocation, pdir, d, "teststats.txt"))
    }
    fdr_temp$V3 <- d
    fdr_temp$V4<-pdir
    fdr <- rbind(fdr, fdr_temp)
    
    
  }
}

annotation_table<-data.table(unique(fdr$V4))
annotation_table<-annotation_table[,tstrsplit(V1, "_")]
annotation_table$orig<-unique(fdr$V4)
if(ncol(annotation_table)==14){
colnames(annotation_table)<-c("Network", "Aggr_meth", "FDR", "Perc_perm", "Seed", "Perm_meth", "Sele_meth", "Perm_hi_de", "Deg_hi", "background", "break1", "break2", "window","orig")
annotation_table$Perm_meth<-gsub("degreeawarenodeswap", "degreeaware", annotation_table$Perm_meth)
#fdr<-merge(fdr, annotation_table, by.x = "V4", by.y= "orig")
annotation_table<-annotation_table[,-c("break2")]
}
if(ncol(annotation_table)==13){
  colnames(annotation_table)<-c("Network", "Aggr_meth", "FDR", "Perc_perm", "Seed", "Perm_meth", "Sele_meth", "Perm_hi_de", "Deg_hi", "background", "break1", "window","orig")
  annotation_table$Perm_meth<-gsub("degreeawarenodeswap", "degreeaware", annotation_table$Perm_meth)
  #fdr<-merge(fdr, annotation_table, by.x = "V4", by.y= "orig")
}
fwrite(fdr, file=file.path(outdir, "fdr.tsv"), sep="\t")

network_jaccard <-
  data.frame(
    jaccard_dists = c(),
    union = c(),
    intersection = c(),
    graph = c(),
    filename = c(),
    run = c()
  )
complete_graphs<-data.frame()
counter=1
for(pdir in dir(filelocation)){
  for (gf in dir(file.path(filelocation,pdir))) {
    if(file.exists(paste0(goldstandard, "/", gsub("result", "expression_type", gf), ".tsv")))
      gold <-fread(paste0(goldstandard, "/", gsub("result", "expression_type", gf), ".tsv"))
    gold <- gold[V2 %in% c("differential")]
    gold <- gold[, 1]
    colnames(gold) <- c("gene")
    
    d.file <- gsub("expression_type", "result", gf)
    if(!file.exists(file.path(filelocation, pdir, d.file, "file.graph"))){
      next
    }
    if(!file.exists(file.path(filelocation, pdir, gf, "file.stat"))){
      next
    }
    
    f1 <- fread(file.path(filelocation, pdir, d.file, "file.graph"))
    #f1[, .N, by = V1]
    
    f1$file<-gf
    f1$run<-pdir
    complete_graphs<-rbind(complete_graphs,f1)
    
    #su <- unique(c(f1$V2, f1$V4))
    #intersect_1 <- length(which(su %in% gold$gene))
    #union_1 <- length(unique(c(su, gold$gene)))
    #index_j <- intersect_1 / union_1
    #overall_jaccard[[gf]] <- index_j
    
    jaccard_dists <- sapply(unique(f1$V1), function(g) jaccard(g, f1))
    jaccard_dists<-data.table(jaccard_dists =jaccard_dists[1,], union = jaccard_dists[2,], intersection=jaccard_dists[3,], graph=jaccard_dists[4,])
    jaccard_dists$filename<-gf
    jaccard_dists$run<-pdir
    
    print(counter)
    counter<-counter+1
    print(file.path(filelocation, pdir, gf, "file.stat"))
   
    pval <- fread(file.path(filelocation, pdir, gf, "file.stat"))
    nodes <- fread(file.path(filelocation,pdir,  gf, "file.nodes"))
    gg<-nodes[, .N, by=V2]
    colnames(gg)<-c("V2", "nrnodes")
    myframe<-merge(pval, gg, by.x = "V1", by.y="V2")
    myframe[,norm:=V3/nrnodes]
    #myframe<-myframe[order(norm)]
    myframe$filename<-gf
    myframe$run<-pdir
    jaccard_dists$graph<-unlist(jaccard_dists$graph)
    jaccard_dists<-merge(jaccard_dists, myframe, by.x =c("graph", "filename", "run"), by.y = c("V1", "filename", "run"))
    network_jaccard<-rbind(network_jaccard, jaccard_dists)
  }
}

network_jaccard$jaccard_dists <-as.numeric(network_jaccard$jaccard_dists)
network_jaccard$union <- as.numeric(network_jaccard$union)
network_jaccard$intersection <-as.numeric(network_jaccard$intersection)

#network_jaccard<-merge(network_jaccard, annotation_table, by.x="run", by.y="orig")


sample_graphs<-dir(direct)
sample_graphs<-grep(".nodes", sample_graphs, value = T)
complete_networks<-data.frame(node=c(), network_size=c(), node_value=c(), degree=c())
for(sample in sample_graphs){
  sub<-fread(file.path(direct, sample))
  complete_networks<-rbind(complete_networks, sub)
}


complete_graphs<-as.data.table(complete_graphs)
x<-unlist(stri_split(complete_graphs$file,fixed="_", simplify = T))
x<-as.data.frame(x)
x<-x[,c(1,2,4,6)]
colnames(x)<-c("res", "g1", "g2", "g3")
complete_graphs<-cbind(complete_graphs, x[, 2:4])
complete_graphs$g1<-as.numeric(as.character(complete_graphs$g1))
complete_graphs$g2<-as.numeric(as.character(complete_graphs$g2))
complete_graphs$g3<-as.numeric(as.character(complete_graphs$g3))

d1<-c()
d2<-c()
d3<-c()
count<-complete_graphs[, .N, by=c("V1","file", "run", "g1", "g2", "g3")]
for(t in 1:(nrow(count))){
  triple<-count[t,]
  print(t)
  s<-complete_graphs[V1==triple$V1 & file==triple$file& run==triple$run]
  nod<-data.table(nodes=unique(c(s$V2, s$V4)))
  d1<-c(d1,length(which(nod$nodes %in% unique(complete_networks[V2==triple$g1]$V1))))
  d2<-c(d2,length(which(nod$nodes %in% unique(complete_networks[V2==triple$g2]$V1))))
  d3<-c(d3,length(which(nod$nodes %in% unique(complete_networks[V2==triple$g3]$V1))))
  
}
d<-data.table(d1,d2,d3)
count<-cbind(count, d)
colnames(count)<-c(colnames(count)[1:6], "nrinteractions", colnames(count)[8:10])

count$perc1<-count$d1/count$g1
count$perc2<-count$d2/count$g2
count$perc3<-count$d3/count$g3

network_jaccard<-merge(network_jaccard, count, by.=c("run", "filename", "graph"), by.y=c("run", "file", "V1"))

fwrite(annotation_table, file=file.path(outdir, "annotation.tsv"), sep="\t")
fwrite(network_jaccard, file= file.path(outdir, "evaluation.tsv"), sep="\t")
