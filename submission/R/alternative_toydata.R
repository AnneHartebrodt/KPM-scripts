require(data.table)
require(VennDiagram)
require(ggplot2)
require(bvenn)

network<-"StringDB"
#network.filename<-"Homo_sapiens_biogrid.tsv"
network.filename<-"Homo_Sapiens_String.tsv"
base.dir<-file.path("/home/anne/Masterarbeit/pipeline_new")
file.location <- file.path(base.dir,"sample_networks", network)
## Load data: networks and mapping files
net <- fread(file.path("~/Masterarbeit/data/networks/", network, network.filename))
out.dir<-file.path(base.dir,"sample_data", network)

dir.create(out.dir, recursive = T)
files <- dir(file.location)
#remove graph of size 0
files<-grep("0.", files, invert = T, value = T)
graphs <- grep(".graph", files, perl = T)


for (ind in 1:10) {
  # randomly sample 3 out of the previously made, non overlapping subnetworks
  sam <- sample(graphs, 3, replace = F)
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
  
  for (i in 1:(length(union_list)-1)) {
    for (j in (i + 1):length(union_list)) {
      print(paste0(i, " ", j, " ", length(
        which(union_list[[j]] %in% union_list[[i]])
      )))
    }
  }
  
  differential <- unique(c(graph_table$V2, graph_table$V4))
  nondifferential <-
    unique(c(net$V1, net$V3))[!unique(c(net$V1, net$V3)) %in% differential]
  
  diff <-
    lapply(differential, function(x)
      abs(rnorm(8, mean = 0.0, sd = 0.005)))
  diff <- as.data.table(diff)
  diff <- as.data.table(t(diff))
  
  nodiff <- lapply(nondifferential, function(x)
    abs(runif(8, 0, 1)))
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
}
#ggplot(toydata, aes(x = V2)) + geom_histogram(bins = 100)
