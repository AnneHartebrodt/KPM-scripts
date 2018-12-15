require(data.table)

pvals<-fread("Masterarbeit/differential_out/dataHD/DESEQ/results.tsv")
#pvals<-fread("/home/anne/Masterarbeit/data/HDarray/processed/pvalue_diff_exp_norm.txt")
#colnames(pvals)<-c("V1", "pvalue")

pvals<-pvals[order(pvalue)]
pvals<-pvals[!(is.na(pvalue))]
pvals<-pvals[!(is.na(padj))]
pvals$zscore<-qnorm(pvals$pvalue)

require(data.table)
require(VennDiagram)
require(ggplot2)
require(bvenn)


#network<-"StringDB"
network<-"biogrid"
network.filename<-"Homo_sapiens_biogrid.tsv"
#network.filename<-"Homo_Sapiens_String.tsv"
base.dir<-file.path("/home/anne/Masterarbeit/simplistic/")
file.location <- file.path(base.dir,"sample_networks", network)
## Load data: networks and mapping files
net <- fread(file.path("~/Masterarbeit/data/networks/", network, network.filename))
out.dir<-file.path(base.dir,"sample_data", network)

dir.create(out.dir, recursive = T)
files <- dir(file.location)
#remove graph of size 0
files<-grep("0.", files, invert = T, value = T)
graphs <- grep(".graph", files, perl = T)
#graphs <-graphs[as.numeric(gsub(".graph", "", grep(".graph", files, perl = T, value = T)))>15]
m1<-0.0001
for(m2 in c(0.5, 0.4, 0.3, 0.2)){
  print(m2)
for (ind in 1:5) {
  # randomly sample 3 out of the previously made, non overlapping subnetworks
  for(s in c(1,3,5)){
  sam <- sample(graphs, s, replace = F)
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
  }
  graph_table <- rbindlist(grap_list)
  
  # for (i in 1:(length(union_list)-1)) {
  #   for (j in (i + 1):length(union_list)) {
  #     print(paste0(i, " ", j, " ", length(
  #       which(union_list[[j]] %in% union_list[[i]])
  #     )))
  #   }
  # }
  
  differential <- unique(c(graph_table$V2, graph_table$V4))
  nondifferential <-
    unique(c(net$V1, net$V3))[!unique(c(net$V1, net$V3)) %in% differential]
  
  n<-nrow(pvals[pvalue<0.01])
  diff <-
    lapply(differential, function(x) abs(rnorm(8, mean = m1, sd = 0.000001)))
  diff <- as.data.table(diff)
  diff <- as.data.table(t(diff))
  
  nodiff <- lapply(nondifferential, function(x) abs(rnorm(8, mean = m2, sd = 0.1)))
  nodiff <- as.data.table(nodiff)
  nodiff <- as.data.table(t(nodiff))
  
  #write files
  toydata <- rbind(diff, nodiff)
  toydata$gene_name <- c(differential, nondifferential)
  toydata <- data.table(toydata$gene_name, toydata[, 1:8])
  
  type <- data.table(gene = toydata$V1,
                     type = c(
                       rep("differential", length(differential)),
                       rep("nondifferential", length(nondifferential))
                     ))
  dir.create(file.path(out.dir,
              paste0(m1, "_", m2)))
  
  filename_samplewise <-
    file.path(
      out.dir, paste0(m1, "_", m2),
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
      paste0(m1, "_", m2),
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
      paste0(m1, "_", m2),
      paste0("expression_type_",
             paste0(sam_size, collapse = "_"),
             ".tsv"
      )),
    sep = "\t",
    col.names = F ,
    row.names = F
  )
}
#ggplot(toydata, aes(x = V2)) + geom_histogram(bins = 100)
}
}
