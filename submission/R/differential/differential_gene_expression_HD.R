require(data.table)
require(ggplot2)
require(DESeq2)
require(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
library(stringr)
require(optparse)
require(cowplot)
require(gProfileR)
require(graph)
require(ggrepel)

## Remove version numbers of Transcripts form files
## cd "/home/anne/Documents/Master/MA/data/dataHD/quants/"
## bash ../../../scripts/salmon/cut_name.txt 

option_list = list(
  make_option(c("-q", "--quant"), type="character", default=NULL, 
              help="quant dir", metavar="character"),
  make_option(c("-m", "--map"), type="character", default=NULL, 
              help="identifier mapping", metavar="character"),
  make_option(c("-r", "--run"), type="character", default=NULL, 
              help="runtable", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-t", "--kpmin"), type="character", default=NULL, 
              help="kpm input directory", metavar="character"),
  make_option(c("-n", "--network"), type="character", default=NULL, 
              help="input network", metavar="character"),
  make_option(c("-p", "--pathway"), type="character", default=NULL, 
              help="input network", metavar="character"),
  make_option(c("-a", "--map2"), type="character", default=NULL, 
              help="input network", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

runtable<-opt$run
quants.dir<-opt$quant
out.dir<-opt$out
map.file<-opt$map
network<-opt$network
pathway_id<-opt$pathway
kpm.in<-opt$kpmin
maap<-opt$map2

#runtable<-"Masterarbeit/data/salmon/SraRunTable(2).txt"
#quants.dir<-"/home/anne/Documents/Master/MA/data/dataHD/quants/"
#out.dir<-"/home/anne/Documents/Master/MA/differential_out/dataHD/DESEQ_new/"
#network<-"Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_NCBI900.tsv"
#pathway_id<-"hsa05016"
#map.file<-"Masterarbeit/data/identifier_mappings/mart_export.txt"
#maap<-"Masterarbeit/data/identifier_mappings/mart_export_names.txt"
#kpm.in<-"Masterarbeit/real_data/data/"
filename<-file.path(out.dir,"results.tsv")

dir.create(out.dir, recursive = T)

##Transcript mapping
txdb<-EnsDb.Hsapiens.v86
tx2gene<-transcripts(txdb, return.type = "DataFrame")
tx2gene<-data.table(TXNAME = tx2gene$tx_id, GENEID = tx2gene$gene_id)


#tx2gene<-merge(tx2gene, mapping, by.x="GENEID", by.y="Gene stable ID")
#tx2gene<-tx2gene[,c("TXNAME","NCBI gene ID")]
#colnames(tx2gene)<-c("TXNAME","GENEID" )

# get all the filenames from the directory
dirlist<- list.dirs(quants.dir, recursive = F)
names<-sapply(basename(dirlist), function(x) gsub("_quant", "", x))
runs<-sapply(1:length(names), function(x) file.path(dirlist[x],"quant.sf"))



##Read and create sample information
age<-function(a){
  if(a<=45){
    return (1)
  }
  else if(a>45 & a<=60){
    return (2)
  }
  else if(a>60 & a<=75){
    return (3)
  }
  else{
    return(4)
  }
}
diag<-function(a){
  if(a=="Huntington's Disease"){
    return ("HD")
  }
  else{
    return("NN")
  }
}

##create annotation, read file, recode age and diagnosis to usable values
c<-fread(runtable)
c<-c[, c("Run", "age_of_death", "age_of_onset", "diagnosis", "rin")]
c$age_class<-sapply(c$age_of_death, function(x) age(x))
samples<-data.table(sample=names, run = runs)
ex<-sapply(samples$run, function(x) file.exists(x))
samples<-samples[ex]
samples<-merge(samples, c, by.x="sample", by.y="Run")
samples<-as.data.frame(samples)
rownames(samples)<-samples$sample
samples$diagnosis<-sapply(samples$diagnosis, function(x) diag(x))
samples$diagnosis<-as.factor(samples$diagnosis)
samples$age_class<-as.factor(samples$age_class)


##import count data
txi.salmon <- tximport(samples$run, type = "salmon", tx2gene = tx2gene)
samples<-samples[,c(5,6,7)]


#run DE analysis
ddsTxi<-DESeqDataSetFromTximport(txi.salmon, colData = samples, design = as.formula("~age_class+rin+diagnosis"))
ddsTxi$diagnosis<-relevel(ddsTxi$diagnosis, ref="NN")
dds<-DESeq(ddsTxi)
resLFC <- lfcShrink(dds, coef=2)
res<-results(dds)
resOrdered <- res[order(res$padj),]

#write results to file

write.table(resOrdered, filename , sep = "\t" )


#map ENSG to Entrez gene ids and create additional result file

#filename<-"/home/anne/Masterarbeit/differential_out/dataHD/DESEQ_new//results.tsv"
#out.dir<-"/home/anne/Masterarbeit/differential_out/dataHD/DESEQ/"
#maap<-fread()
#maap<-möp[-which(is.na(maap$`Gene stable ID`))]
maap<-fread(maap)
re<-fread(filename)
colnames(re)<-c("geneName","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
re<-merge(re, maap, by.x="geneName" ,by.y="Gene stable ID", all.x=T)
na<-which(is.na(re$`Gene name`))
re$`Gene name`[na]<-re$geneName[na]

re<-re[-which(is.na(padj))]
re$'P-Value'<-as.factor(as.integer(re$padj<0.001&abs(re$log2FoldChange)>2))
levels(re$'P-Value')<-c("<0.001\nabs(FC)<2",">0.001\n(abs(FC)>2")
g1<-ggplot(re, aes(x=log2FoldChange, y=-log10(padj), label=geneName))+
  geom_point(aes(color=`P-Value`))+
  xlab("Log2 Fold Change")+ylab("-Log10(Adjusted P-Value)")+
  ggtitle("Differential Gene Expression Analysis: Significant Genes")+
  scale_color_viridis_d(option="E")+theme(legend.title = element_blank())+
 geom_text_repel(aes(label=ifelse((re$padj<0.001&abs(re$log2FoldChange)>10)|re$padj<10^-10,`Gene name`,'')), 
           size=2)+
  ylim(c(0,20))+theme(legend.position = "bottom")
#g1
ggsave(filename = file.path(out.dir, "volcano.pdf"), plot = g1,
       width = 20, height = 20, units = "cm")

g2<-ggplot(re, aes(x=pvalue, fill=I("#00204DFF")))+geom_histogram(bins=100)+
  ggtitle("P-value distribution of Huntington's Disease DE analysis")+
  xlab("P-Value")+ylab("Count")+scale_color_viridis_d(option="E")
ggsave(filename = file.path(out.dir, "histp.pdf"), plot = g2,
       width = 20, height = 20, units = "cm")

g3<-ggplot(re, aes(x=padj, fill=I("#00204DFF")))+geom_histogram(bins = 100)+
  ggtitle("P-value (adjusted) distribution of Huntington's Disease DE analysis")+
  xlab("P-Value")+ylab("Count")+scale_color_viridis_d(option="E")
ggsave(filename = file.path(out.dir, "histpadj.pdf"), plot = g3, 
       width = 20, height = 20, units = "cm")




re<-fread(filename)
colnames(re)<-c("geneName","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
#identifier mapping
mapping<-fread(map.file)

mapping<-mapping[!is.na(`Protein stable ID`)]
mapping<-mapping[!is.na(`NCBI gene ID`)]
mapping<-mapping[!is.na(`Gene stable ID`)]

mapping<-unique(mapping[,.(`Gene stable ID`, `NCBI gene ID`)])
re<-merge(re, mapping, by.x="geneName", by.y="Gene stable ID")
re<-re[,c(8,2:7)]

# deduplicate using max padj = conservative
max<-re[,.SD[which.min(padj)], by="NCBI gene ID"]
file.NCBI<- file.path(out.dir, "results_NCBI.tsv")
write.table(max, file.path(out.dir, "results_NCBI.tsv"), sep = "\t" , row.names = F)


# Write KPM compatible file

re<-fread(filename)
colnames(re)<-c("geneName","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
dat<-re[,.(geneName, padj)]
na<-which(is.na(dat$padj))
dat<-dat[-na]
fwrite(file= file.path(kpm.in, "expression_values.txt"), x = dat, sep="\t", col.names = F)
dat$binary<-as.integer(dat$padj<0.001)
bin<-dat[,.(geneName, binary)]
fwrite(file= file.path(kpm.in, "expression_values_bin.txt"), x = bin, sep="\t", col.names = F)



#file.NCBI<-"/home/anne/Masterarbeit/differential_out/dataHD/DESEQ_final/results_NCBI.tsv"
re<-fread(file.NCBI)
colnames(re)<-c("geneName","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
#kpm.in<-"Masterarbeit/kpm_hd/data/"
#pathway_id<-"hsa05016"
#out.dir<-"/home/anne/Masterarbeit/differential_out/dataHD/DESEQ/"
#network<-"Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_NCBI900.tsv"

#Load Kegg pathway for Huntingtons' Disease
kegg<-KEGGREST::keggGet(pathway_id)
kegglist<-unlist(kegg[[1]]$GENE)
ind1<-seq(1,386,2)
ind2<-seq(2,386,2)
kegglist<-data.table(gene_name=kegglist[ind1], desc=kegglist[ind2])

top<-re[padj<0.001]

overlap<-kegglist[kegglist$gene_name %in% top$geneName]
#write.table(overlap, file.path(out.dir, "KEGG_overlap.tsv"), quote = F, row.names = F)
#Load String Database
String700<-fread(network)

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
#out.dir<-"/home/anne/Masterarbeit/differential_out/dataHD/DESEQ/"
fwrite(file = file.path(out.dir, "KEGG_overlap.tsv"), overlap)

# re<-fread(filename)
# colnames(re)<-c("geneName","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
# top<-re[padj<0.001]
# gsig<-gprofiler(top$geneName, organism = "hsapiens", ordered_query = T)
# gsig$size<-sapply(gsig$intersection, function(x) length(unlist(strsplit(x, ","))))
# 

re<-fread(file.NCBI)
colnames(re)<-c("geneName","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
r1<-re[, .(geneName, padj)]
r1$rank<-order(r1$padj)
r2<-re[,.(geneName, log2FoldChange)]
r2$rank<-order(r2$log2FoldChange)
rall<-merge(r1,r2, by="geneName")
rall$rankproduct<-rall$rank.x*rall$rank.y
fwrite(file = file.path(out.dir, "rank_product.tsv"), rall)
