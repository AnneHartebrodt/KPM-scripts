require(data.table)
require(graph)
require(ggplot2)
require(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="sif submoule", metavar="character"),
  make_option(c("-n", "--network"), type="character", default=NULL, 
              help="Network", metavar="character"),
  make_option(c("-p", "--pathway"), type="character", default=NULL, 
              help="pathway id", metavar="character"),
  make_option(c("-m", "--mapping"), type="character", default=NULL, 
              help="id mapping file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output folder name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#network<-"Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_NCBI700.tsv"
network<-opt$network
map<-opt$mapping
#map<-" "
#pathway_id<-"hsa05016"
pathway_id<-opt$pathway
module<-opt$file
#module<-"Masterarbeit/bionet_out/dataHD/StringDB700/module.sif"
out<-opt$out
#out<-"Masterarbeit/bionet_out/dataHD/StringDB700/overlap.tsv"

#Load String Database
String700<-fread(network, header = F)
#load mapping
mapping<-fread(map)

#Load Kegg pathway for Huntingtons' Disease
kegg<-KEGGREST::keggGet(pathway_id)
kegglist<-unlist(kegg[[1]]$GENE)
ind1<-seq(1,386,2)
ind2<-seq(2,386,2)
kegglist<-data.table(gene_name=kegglist[ind1], desc=kegglist[ind2])

#Load module
moduleString700<-fread(module, header = F)
colnames(moduleString700)<-c("gene1", "pp", "gene2")
unic_nodes<-unique(c(moduleString700$gene1,moduleString700$gene2))
uniq_interactions<-nrow(unique(moduleString700))

#map unique genes to NCBI
mapping<-unique(mapping[,.(`Gene stable ID`, `NCBI gene ID`)])
unic_nodes_NCBI<-mapping[`Gene stable ID` %in% unic_nodes]$`NCBI gene ID`

#overlap of genes in KEGG pathway and Bionet module
overlap<-kegglist[as.numeric(kegglist$gene_name) %in% unic_nodes_NCBI]

#create graph object from submodule containing the overlap genes
sub<-String700[V1 %in% overlap$gene_name | V3 %in% overlap$gene_name]
sub$V1<-as.character(sub$V1)
sub$V3<-as.character(sub$V3)
g<-graphNEL()
g<-addNode(node = unique(as.character(c(sub$V1, sub$V3))), g)
g<-addEdge(as.vector(sub$V1), as.vector(sub$V3), g, 1)

#calculate the degree of the overlap nodes
overlap$degree<-degree(g)[as.character(overlap$gene_name)]
overlap[,c("Gene Name", "Description"):=tstrsplit(desc, ";")]

colnames(overlap)<-c("NCBI gene ID", "Gene Name and Description", "Degree", "Name",
                     "Description")
overlap<-overlap[, c("NCBI gene ID", "Degree", "Name",
            "Description")]

fwrite(file=file.path(out, "overlap.tsv"), overlap)
