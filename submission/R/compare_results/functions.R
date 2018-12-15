#require(data)
extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  #grobs <- ggplot_gtable(gg)
  
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}
make_anno<-function(annotation_table, network_jaccard){
  print(colnames(annotation_table))
  annotation_table$Aggr_meth<-gsub("meanLog", "mean(-{log(p)})", annotation_table$Aggr_meth)
  annotation_table$Aggr_meth<-gsub("sum", "sum({p})", annotation_table$Aggr_meth)
  annotation_table$Aggr_meth<-gsub("normDegSum", "sum({p*deg(p)})", annotation_table$Aggr_meth)
  annotation_table$Perm_meth<-gsub("edgerewire", "edge-rewire", annotation_table$Perm_meth)
  annotation_table$Perm_meth<-gsub("degreeaware", "degree-aware\nnodeswap", annotation_table$Perm_meth)
  annotation_table$break1<-gsub("largest", "largest\nsubnetwork", annotation_table$break1)
  annotation_table$break1<-gsub("sliding", "sliding window", annotation_table$break1)
  annotation_table$break1<-gsub("maximum", "max. distance\nfrom theshold", annotation_table$break1)
  colnames(annotation_table)<-c("Network","Aggregation\nMethod", "FDR", "%Network\nPermutation", 
                                 "Seed", "Permutation\nMethod", "Selection\nmethod",
                                "#High Degree\nNodes Permuted",
                                 "Degree of\nHi.deg Node", "Background",
                                "Termination\ncriterion", "Sliding\nwindow size", "orig")
  aa<-melt(annotation_table, id.vars = "orig")
  aa$variable<-as.factor(aa$variable)
  aa$group <- paste0(aa$variable, aa$value)
  a <-sapply(colnames(annotation_table), function(x) length(unlist(unique(annotation_table[,x, with=F])))==1)
  a<-names(which(!a))
  a<-grep("orig", a, value = T, invert = T)
  print(a)

  if(length(a)!=0){
    annotation_table<-annotation_table[orig %in% network_jaccard[, .N, by="run"]$run]
    #palettes<-c("BuGn", "Reds", "OrRd", "PuBu", "PuRd", "PuBu","BuGn", "Reds", "OrRd", "PuBu")
    n<-10
    df<-data.table(col1=c(viridis(n)), col2=c(viridis(n,0.2)),col3=c( viridis(n,0.4)),col4=c( viridis(n,0.6)), col5=c(viridis(n,0.8)))
    #df<-data.table(col1=c(plasma(n/2), viridis(n/2)), col2=c(plasma(n/2,0.2), viridis(n/2,0.2)),col3=c(plasma(n/2,0.4), viridis(n/2,0.4)),col4=c(plasma(n/2,0.6), viridis(n/2,0.6)), col5=c(plasma(n/2,0.8), viridis(n/2,0.8)))
    df<-rbind(df[seq(1,n,2)], df[seq(2,n,2)])
    blank_plot<-list()
    tile<-list()
    for(ai in 1:length(a)){
      #pal<-brewer.pal(7, palettes[ai])
      #pal<-pal[2:7]
      annotation_table[, a[ai]:=as.factor(get(a[ai]))]
      print(a[ai])
      print(annotation_table[,c(get(a[ai]))])
      ggg<-ggplot(annotation_table)+geom_tile(aes(x=orig, y="", fill=as.factor(get(a[ai]))),position = "identity")+
        scale_color_manual(aesthetics = "fill",values = as.character(unlist(df[ai])))+
      guides(fill=guide_legend(title=a[ai]))
      #plot(ggg)
      leg<-extractLegend(ggg)
      x<-which(colnames(annotation_table)==a[ai])
      m<--1+(length(unique(unlist(annotation_table[,..x])))*0.2)
      blank_plot[[ai]]<-plot_to_gtable(ggplot()+geom_blank()+annotation_custom(leg, ymin=m))
      tile[[ai]]<-plot_to_gtable(ggg+theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),plot.margin=unit(c(0,0,0,0), unit="cm"),  axis.line = element_blank())+guides(fill=F))
    #plot(tile[[ai]])
      }
    print("creating plot")
    #plot.margin=unit(c(0,0,0,0), unit="cm"),
    legend<-do.call(plot_grid, c(blank_plot, nrow=1, align="none"))
    plo<-do.call(plot_grid, c(tile, align="v", ncol=1))
    tt<<-tile
  }
  else{
    plo<-ggplot()+geom_blank()
    legend<-ggplot()+geom_blank()
  }
  return(list(plo, legend))
}
makePlot<-function(annotation_table, network_jaccard){
  g1 <-ggplot(network_jaccard, aes(y = jaccard_dists, x = run)) + geom_boxplot() +
    ggtitle("Jaccard Indices of Foreground and Solution Subnetworks") + ylab("Jaccard distance") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
          plot.margin=unit(c(0,0,0,0), unit="cm"), plot.title = element_text(size=30), 
          axis.text=element_text(size=30), axis.title =element_text(size=20))+
    guides(fill=F)
  
  a<-make_anno(annotation_table, network_jaccard)
  plo<-a[[1]]
  legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.6,0.2,0.3)))
}
makePlot2<-function(annotation_table, network_jaccard){
  g1 <-ggplot(network_jaccard, aes(y = percMean, x = run)) + geom_boxplot() +
    ggtitle("Fraction of Foreground Genes Found in Solution") + ylab("Fraction") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
          plot.margin=unit(c(0,0,0,0), unit="cm"), plot.title = element_text(size=30),
          axis.text=element_text(size=30),  axis.title =element_text(size=20))+
    guides(fill=F)
  a<-make_anno(annotation_table, network_jaccard)
  plo<-a[[1]]
  legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.5,0.25,0.25)))
}
