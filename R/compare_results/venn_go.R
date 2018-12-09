go.dir<-"/home/anne/Masterarbeit/evalutation_all_methods/go/"
probal.1.go.file<-fread(file.path(go.dir, "pkpm_best.tsv"), sep="\t")
bionet.go.file<-fread(file.path(go.dir, "bionet.tsv"), sep="\t")
diff.go.file<-fread(file.path(go.dir, "deseq_top.tsv"), sep="\t")
web.go.file<-fread(file.path(go.dir, "KPM_k5.tsv"), sep="\t")
#bionet.700.go
probal.1.go<-probal.1.go.file[V2<0.01]$V1
bionet.go<-bionet.go.file[V2<0.01]$V1
diff.go<-diff.go.file[V2<0.01]$V1
web.go<-web.go.file[V2<0.01]$V1




bionet.kpm<-setdiff(intersect(bionet.go, web.go), diff.go)
bionet.kpm<-bionet.go.file[V1 %in% bionet.kpm][order(V2)][1:20]
onto2<-web.go.file[V1 %in% bionet.kpm][order(V2)][1:20]
#onto1<-merge(onto1, onto2, by="V1", all=T)
fwrite(bionet.kpm,file=file.path(go_dir, "bionet.web.tsv"), sep="\t")

diff.pkpm<-setdiff(setdiff(intersect(diff.go, probal.1.go), bionet.go),web.go)
diff.pkpm<-probal.1.go.file[V1 %in% diff.pkpm][order(V2)][1:20]
onto22<-diff.go.file[V1 %in% diff.pkpm][order(V2)][1:20]
#merge(onto11, onto22, by="V1", all=T)
fwrite(diff.pkpm,file=file.path(go_dir, "diff.pkpm.tsv"), sep="\t")


bionet.only<-setdiff(setdiff(setdiff(bionet.go, web.go),diff.go),probal.1.go)
bio.onl<-bionet.go.file[V1 %in% bionet.only][order(V2)]
fwrite(bio.onl,file=file.path(go_dir, "bionet.only.tsv"), sep="\t")


web.only<-setdiff(setdiff(setdiff(web.go, bionet.go),diff.go),probal.1.go)
web.onl<-web.go.file[V1 %in% web.only][order(V2)]
fwrite(web.onl,file=file.path(go_dir, "web.only.tsv"), sep="\t")


bio.pkpm.diff<-setdiff(intersect(intersect(bionet.go, diff.go), probal.1.go), web.go)
bio.pkpm.diff<-probal.1.go.file[V1 %in% bio.pkpm.diff][order(V2)]
fwrite(bio.pkpm.diff,file=file.path(go_dir, "bio.pkpm.diff.tsv"), sep="\t")


diff.pkpm.web<- setdiff(intersect(intersect(diff.go, probal.1.go), web.go), bionet.go)
diff.pkpm.web<-probal.1.go.file[V1 %in% diff.pkpm.web][order(V2)]
fwrite(diff.pkpm.web,file=file.path(go_dir, "diff.pkpm.web.tsv"), sep="\t")


bionet.pkpm.diff.web<-intersect(intersect(intersect(diff.go, web.go),probal.1.go), bionet.go)
onto111<-bionet.go.file[V1 %in% bionet.pkpm.diff.web][order(V2)][1:20]
onto222<-probal.1.go.file[V1 %in% bionet.pkpm.diff.web][order(V2)][1:20]
onto333<-diff.go.file[V1 %in% bionet.pkpm.diff.web][order(V2)][1:20]
onto444<-web.go.file[V1 %in% bionet.pkpm.diff.web][order(V2)][1:20]

onto_m<-merge(onto111, onto222, by="V1", all=T)
onto_m$bionet<-as.integer(!is.na(onto_m$V2.x))
onto_m$pkpm<-as.integer(!is.na(onto_m$V2.y))
onto_m$V3.x[is.na(onto_m$V3.x)]<-onto_m$V3.y[is.na(onto_m$V3.x)]
onto_m$pval<-sapply(1:nrow(onto_m), function(x) max(onto_m$V2.x[x], onto_m$V2.y[x], na.rm = T))
onto_m<-onto_m[, c(1,3,6:8)]

onto_m2<-merge(onto333, onto444, by="V1", all=T)
onto_m2$diff<-as.integer(!is.na(onto_m2$V2.x))
onto_m2$web<-as.integer(!is.na(onto_m2$V2.y))
onto_m2$V3.x[is.na(onto_m2$V3.x)]<-onto_m2$V3.y[is.na(onto_m2$V3.x)]
onto_m2$pval<-sapply(1:nrow(onto_m2), function(x) max(onto_m2$V2.x[x], onto_m2$V2.y[x], na.rm = T))
onto_m2<-onto_m2[, c(1,3,6:8)]
onto_m<-merge(onto_m, onto_m2, by="V1")
onto_m$pval<-sapply(1:nrow(onto_m), function(x) max(onto_m$pval.x[x], onto_m$pval.y[x], na.rm = T))
onto_m<-onto_m[, c(1:4, 7:8,10)]

colnames(onto_m)<-c("GO:term", "Description", "Bionet", "pKPM", "DE Genes", "KPM k=5", "max(p-value)")
fwrite(onto_m,file=file.path(go_dir, "all.go.tsv"), sep="\t")


area1<-length(probal.1.go)
area2<-length(bionet.go)
area3<-length(diff.go)
area4<-length(web.go)
#area5<-length(bionet.700.go)

n12<-length(intersect(probal.1.go, bionet.go))
n13<-length(intersect(probal.1.go, diff.go))
n14<-length(intersect(probal.1.go, web.go))
#n15<-length(intersect(probal.1.go, bionet.700.nodes$nodes))

n23<-length(intersect(bionet.go, diff.go))
n24<-length(intersect(bionet.go, web.go))
#n25<-length(intersect(bionet.go, bionet.700.nodes$nodes))

n34<-length(intersect(diff.go, web.go))
#n35<-length(intersect(diff.go, bionet.700.nodes$nodes))

#n45<-length(intersect(web.go, bionet.700.nodes$nodes))

n123<-length(intersect(intersect(probal.1.go, bionet.go), diff.go))
n124<-length(intersect(intersect(probal.1.go, bionet.go), web.go))
#n125<-length(intersect(intersect(probal.1.go, bionet.go), bionet.700.nodes$nodes))
n134<-length(intersect(intersect(probal.1.go, diff.go), web.go))
#n135<-length(intersect(intersect(probal.1.go, diff.go), bionet.700.nodes$nodes))
#n145<-length(intersect(intersect(probal.1.go, web.go), bionet.700.nodes$nodes))

#n235<-length(intersect(intersect(bionet.go, diff.go), bionet.700.nodes$nodes))
n234<-length(intersect(intersect(bionet.go, diff.go), web.go))
#n245<-length(intersect(intersect(bionet.go, web.go ), bionet.700.nodes$nodes))

#n345<-length(intersect(intersect(diff.go, web.go),bionet.700.nodes$nodes))

n1234<-length(intersect(intersect(intersect(bionet.go, web.go), diff.go), probal.1.go))
#n1235<-length(intersect(intersect(intersect(bionet.go, bionet.700.nodes$nodes), diff.go), probal.1.go))
#n1245<-length(intersect(intersect(intersect(probal.1.go, web.go), bionet.700.nodes$nodes), bionet.go))
#n1345<-length(intersect(intersect(intersect(probal.1.go, web.go), bionet.700.nodes$nodes), diff.go))
#n2345<-length(intersect(intersect(intersect(bionet.go, web.go), bionet.700.nodes$nodes), diff.go))
#n12345<-length(intersect(intersect(intersect(intersect(bionet.go, web.go), bionet.700.nodes$nodes), diff.go), probal.1.go))

# 
# tt3 <- ttheme_minimal(base_size = 8, core=list(fg_params=list(hjust=0, x=0.1)),
#                       rowhead=list(fg_params=list(hjust=0, x=0)),
#                       padding = unit(c(2,2), "mm"))
# tt3<-ttheme_minimal( base_size = 8, padding = unit(c(2,2), "mm"))
# 
# g1<-ggplot()+annotation_custom(tableGrob(bio.onl[1:5,c(1,3)], theme=tt3, cols = NULL, rows = NULL),ymin = 0.5)
# g2<-ggplot()+annotation_custom(tableGrob(bio.pkpm.diff[1:5,c(1,3)], theme=tt3, cols = NULL, rows = NULL),ymin = 0.5)
# g3<-ggplot()+annotation_custom(tableGrob(onto_m[1:5,c(1,2)], theme=tt3, cols = NULL, rows = NULL),
#                                xmin = 0.3, ymin = 0.5)
# g4<-ggplot()+annotation_custom(tableGrob(bionet.kpm[1:5,c(1,3)], theme=tt3, cols = NULL, rows = NULL))
# g6<-ggplot()+annotation_custom(tableGrob(diff.pkpm.web[1:5,c(1,3)], theme=tt3, cols = NULL, rows = NULL))
# g7<-ggplot()+annotation_custom(tableGrob(web.onl[1:5,c(1,3)], theme=tt3, cols = NULL, rows = NULL),ymin = -0.5)
# g8<-ggplot()+geom_blank()
# g9<-ggplot()+annotation_custom(tableGrob(diff.pkpm[1:5,c(1,3)], theme=tt3, cols = NULL, rows = NULL),ymin = -0.5)

pdf(file.path("~/Masterarbeit/evalutation_all_methods/", paste0("vennGO.pdf")))
ggplot()+geom_blank()+ggtitle("Overlap between the GO term sets in the best solutions")
grob.venn<-draw.quad.venn(area1, area2, area3, area4, n12, n13, n14,
                          n23, n24,  n34, n123, n124, n134,  n234, n1234,
                          category = c("pKPM", "Bionet", "DE Genes", "KPM\nk=5"),
                          cat.pos = c(350, 10,0,0), cat.dist =c(-0.25, -0.24,-0.15,-0.15) , 
                          fill = c("#eff3ff","#bdd7e7","#6baed6","#3182bd" ), rotation.degree = 180, )
dev.off()

 




