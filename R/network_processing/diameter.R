require(igraph)
require(data.table)

string<-fread("Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS900.tsv")
s<-as.matrix(string[,c(1,3)])
g<-graph_from_edgelist(s, directed = F)

diam_900<-diameter(g)

string<-fread("Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS700.tsv")
s<-as.matrix(string[,c(1,3)])
g<-graph_from_edgelist(s, directed = F)
diam_700<-diameter(g)
