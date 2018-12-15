library(BioNet)
library(DLBCL)
require(data.table)
require(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="differential expression dataset file name", metavar="character"),
  make_option(c("-n", "--network"), type="character", default=NULL, 
              help="network file in sif format", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

outdir<-opt$out
infile<-opt$file
innet<-opt$network
#outdir<-"Masterarbeit/bionet_out/HD/StringDB/"
#infile<-"Masterarbeit/data/differential_out/HD/DESEQ/results.tsv"
#innet<-"Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS700.tsv"

dir.create(outdir, recursive = T)

start<-Sys.time()
#create graphNEL of Network
dat<-fread(innet, header = F)
dat$V1<-as.character(dat$V1)
dat$V3<-as.character(dat$V3)
dat<-dat[!(V1==V2)]
g<-graphNEL()

g<-addNode(node = unique(c(dat$V1, dat$V3)), g)
g<-addEdge(as.vector(dat$V1), as.vector(dat$V3), g, 1)

#pvana<-fread("Masterarbeit/data/DESeq_out/DESeq_Pvals_Entrez.tsv")
pvalues<-fread(infile)
pvalue<-as.numeric(pvalues$pvalue)
names(pvalue)<-pvalues$V1
pvalue<-pvalue[!is.na(pvalue)]

# Program throws an error if you try to select non exitent nodes
# from the graph, so make sure, that they are actually in the graph
pvalue<-pvalue[names(pvalue) %in% names(nodeData(g))]

g<-subGraph(names(pvalue), g)

pdf(file.path(outdir, "BUM_fit.pdf"))
fb <- fitBumModel(pvalue, plot = TRUE)
dev.off()

pdf(file.path(outdir,"log_likelihood_surface.pdf"))
plotLLSurface(pvalue, fb)
dev.off()

#score Nodes
scores <- scoreNodes(g, fb, fdr = 0.001)

#Run Heinz, only heuristic solution works due to installation problems.
starttime<-Sys.time()
print(starttime)
#module<-runHeinz(g, scores, heinz.folder = "Masterarbeit/data/h")
module <- runFastHeinz(g, scores)
endtime<-Sys.time()
print(endtime)

# log runtime 
duration<-endtime - starttime
cat(paste0("Running time of heinz:", "\t" ,duration, "\n"), file = file.path(outdir, "logger.txt"))


#plotModule(module, scores = scores)
saveNetwork(module, file=file.path(outdir,"module.sif"), type="sif")
end<-Sys.time()
d<-end-start
cat(paste0("Total run time:", "\t" ,d, "\n"), file = file.path(outdir, "logger.txt"), append = T)
