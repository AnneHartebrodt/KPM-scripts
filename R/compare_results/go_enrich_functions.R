go_list<-function(){
g<-fread("~/Masterarbeit/data/go/goa_human.gaf", sep="\t", header = F, skip = "!" )
map.file<-"~/Masterarbeit/data/identifier_mappings/mart_export_names.txt"
maap<-fread(map.file)
#only curated entries
go_evd<-c("IDA", "EXP", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP","ISS",
          "ISO", "ISA", "ISAM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "TAS", "NAS")
go_man<-g[V7 %in% go_evd]
#remove non assosciations
go_not<-c("NOT")
go_man<-go_man[!(V4 %in% go_not)]
go_man<-go_man[,.(V2,V3,V5, V11)]
go_man<-merge(go_man, maap, by.x= "V3", by.y="Gene name")

list<-sapply(unique(go_man$`Gene stable ID`), function(x) go_man[`Gene stable ID`==x]$V5)
names(list)<-unique(go_man$`Gene stable ID`)
return(list)
}


annot<-function(result){
go_annot<-fread("~/Masterarbeit/data/go/go_terms.tsv", header=F)
terms<-data.table(name=names(result@score),score=result@score)
terms<-merge(terms, go_annot, by.x="name", by.y = "V1")
return(terms)
}
#GeneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))