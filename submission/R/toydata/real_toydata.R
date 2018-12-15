require(data.table)
require(data.table)
require(VennDiagram)
require(ggplot2)
require(bvenn)
require(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="DESeq outfile", metavar="character"),
  make_option(c("-n", "--network"), type="character", default=NULL, 
              help="network name", metavar="character"),
  make_option(c("-m", "--netfile"), type="character", default=NULL, 
              help="network filename", metavar="character"),
  make_option(c("-b", "--basedir"), type="character", default=NULL, 
              help="base dir", metavar="character"),
  make_option(c("-k", "--netdir"), type="character", default=NULL, 
              help="network directory", metavar="character"),
  make_option(c("-d", "--datasets"), type="character", default=NULL, 
              help="network directory", metavar="character"),
  make_option(c("-r", "--nrdatasets"), type="character", default=NULL, 
              help="network directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pvals<-opt$file
network<-opt$network
network.filename<-opt$netfile
base.dir<-opt$basedir



## Load data: networks and mapping files

netdir<-opt$netdir
net <- fread(file.path(netdir, network.filename))

out.dir<-file.path(base.dir,"sample_data", network)


dd<-as.integer(opt$nrdatasets)
nes<-as.integer(opt$datasets)
#dd<-5
#nes<-2
#network<-"String900"
#network<-"biogrid"
#network.filename<-"Homo_sapiens_biogrid.tsv"
#netdir<-"Masterarbeit/data/networks/StringDB/"
#pvals<-"Masterarbeit/differential_out/dataHD/DESEQ/results.tsv"
#network.filename<-"Homo_Sapiens_String_NCBI900.tsv"
#base.dir<-file.path("/home/anne/Masterarbeit/test_runs")

dir.create(out.dir, recursive = T)
dir.create(file.path(out.dir, "plots"))
file.location <- file.path(base.dir,"sample_networks", network)
files <- dir(file.path(file.location))

pvals<-fread(pvals)
pvals<-pvals[order(pvalue)]
pvals<-pvals[!(is.na(pvalue))]
pvals<-pvals[!(is.na(padj))]
pvals$zscore<-qnorm(pvals$pvalue)
#remove graph of size 0
#files<-grep("0.", files, invert = T, value = T)
graphs <- grep(".graph", files, perl = T)
#graphs <-graphs[as.numeric(gsub(".graph", "", grep(".graph", files, perl = T, value = T)))>15]


for (ind in 1:dd) {
  # randomly sample 3 out of the previously made, non overlapping subnetworks
  sam <- sample(graphs, nes, replace = F)
  sam_size<-sapply(sam, function(x) gsub(".graph", "", files[x]))
  
  grap_list <- list()
  union_list <- list()
  counter = 1
  
  for (g in sam) {
    gr <- fread(file.path(file.location, files[g]))
    grap_list[[counter]] <- gr
    union <- unique(c(gr$V2, gr$V4))
    union_list[[counter]] <- union
    counter = counter + 1
  graph_table <- rbindlist(grap_list)
  }
  # for (i in 1:(length(union_list)-1)) {
  #   for (j in (i + 1):length(union_list)) {
  #     print(paste0(i, " ", j, " ", length(
  #       which(union_list[[j]] %in% union_list[[i]])
  #     )))
  #   }
  # }
  
  differential <- unique(c(graph_table$V2, graph_table$V4))
  noise<-sample(unique(c(net$V1, net$V3))[!unique(c(net$V1, net$V3)) %in% differential], 1500, replace = F)
  nondifferential <-
    unique(c(net$V1, net$V3))[!unique(c(net$V1, net$V3)) %in% c(differential, noise)]
  
  n<-nrow(pvals[pvalue<0.001])
  diff <-
    lapply(differential, function(x) sample(abs(pvals$pvalue[1:n]), 8, replace = F))
  diff <- as.data.table(diff)
  diff <- as.data.table(t(diff))
  
  n<-nrow(pvals[pvalue<0.05])
  noi <-
    lapply(noise, function(x) sample(abs(pvals$pvalue[1:n]), 8, replace = F))
  noi <- as.data.table(noi)
  noi <- as.data.table(t(noi))
  
  nodiff <- lapply(nondifferential, function(x) runif(min=0, max=1,n = 8))
  nodiff <- as.data.table(nodiff)
  nodiff <- as.data.table(t(nodiff))
  
  #write files
  toydata <- rbind(diff, noi, nodiff)
  toydata$gene_name <- c(differential, noise, nondifferential)
  toydata <- data.table(toydata$gene_name, toydata[, 1:8])
  
  type <- data.table(gene = toydata$V1,
                     type = c(
                       rep("differential", length(differential)),
                       rep("noise", length(noise)),
                       rep("nondifferential", length(nondifferential))
                     ))
  filename_samplewise <-
    file.path(
      out.dir,
      paste0("samplewise_",
             paste0(sam_size, collapse = "_"),
             ".tsv")
    )
  write.table(
    toydata,
    filename_samplewise,
    sep = "\t",
    col.names = F,
    row.names = F
  )
  filename_general <-
    file.path(
      out.dir,
      paste0("general_",
             paste0(sam_size, collapse = "_"),
             ".tsv")
    )
  write.table(
    toydata[, 1:2],
    filename_general,
    sep = "\t",
    col.names = F,
    row.names = F
  )
  write.table(
    type,
    file.path(
      out.dir,
      
      paste0("expression_type_",
             paste0(sam_size, collapse = "_"),
             ".tsv"
      )),
    sep = "\t",
    col.names = F ,
    row.names = F
  )
  
  colnames(toydata)<-c("gene", "p", seq(1,7,1))
  g1<-ggplot(toydata, aes(x = p, fill=I("#00204DFF"))) +geom_histogram(bins = 100)+
    ggtitle("P-value distribution of Toydata genes")+
    xlab("P-Value")+ylab("Count")+scale_color_viridis_d(option="E")
  ggsave(filename = file.path(out.dir, "plots", paste0("histogram",
                                              paste0(sam_size, collapse = "_"),
                                              ".pdf"
  )))
  
}
