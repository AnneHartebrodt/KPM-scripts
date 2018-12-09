#require(data)
extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  #grobs <- ggplot_gtable(gg)
  
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}
make_anno<-function(annotation_table, network_jaccard){
  aa<-melt(annotation_table, id.vars = "orig")
  aa$variable<-as.factor(aa$variable)
  aa$group <- paste0(aa$variable, aa$value)
  a <-sapply(colnames(annotation_table), function(x) length(unlist(unique(annotation_table[,x, with=F])))==1)
  a<-names(which(!a))
  a<-grep("orig", a, value = T, invert = T)
  
  if(length(a)!=0){
    annotation_table<-annotation_table[orig %in% network_jaccard[, .N, by="run"]$run]
    #palettes<-c("BuGn", "Reds", "OrRd", "PuBu", "PuRd", "PuBu","BuGn", "Reds", "OrRd", "PuBu")
    n<-10
    df<-data.table(col1=c(plasma(n/2), viridis(n/2)), col2=c(plasma(n/2,0.2), viridis(n/2,0.2)),col3=c(plasma(n/2,0.4), viridis(n/2,0.4)),col4=c(plasma(n/2,0.6), viridis(n/2,0.6)), col4=c(plasma(n/2,0.8), viridis(n/2,0.8)))
    df<-rbind(df[seq(1,n,2)], df[seq(2,n,2)])
    blank_plot<-list()
    tile<-list()
    for(ai in 1:length(a)){
      #pal<-brewer.pal(7, palettes[ai])
      #pal<-pal[2:7]
      annotation_table[, a[ai]:=as.factor(get(a[ai]))]
      ggg<-ggplot(annotation_table)+geom_tile(aes_string(x="orig", y="''", fill=a[ai]),position = "identity")+
        scale_color_manual(aesthetics = "fill",values = as.character(unlist(df[ai])))+
        guides(fill=guide_legend(title=a[ai]))
      leg<-extractLegend(ggg)
      x<-which(colnames(annotation_table)==a[ai])
      m<--1+(length(unique(unlist(annotation_table[,..x])))*0.2)
      blank_plot[[ai]]<-plot_to_gtable(ggplot()+geom_blank()+annotation_custom(leg, ymin=m))
      tile[[ai]]<-ggg+theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),plot.margin=unit(c(0,0,0,0), unit="cm"),  axis.line = element_blank())+guides(fill=F)
    }
    print("creating plot")
    #plot.margin=unit(c(0,0,0,0), unit="cm"),
    legend<-do.call(plot_grid, c(blank_plot, nrow=1, align="none"))
    plo<-do.call(plot_grid, c(tile, align="v", ncol=1))
  }
  else{
    plo<-ggplot()+geom_blank()
    legend<-ggplot()+geom_blank()
  }
  return(list(plo, legend))
}
makePlot<-function(annotation_table, network_jaccard){
  g1 <-ggplot(network_jaccard, aes(y = jaccard_dists, x = run)) + geom_boxplot() +
    ggtitle("Jaccard distances to 'gold standard' toydata genes") + ylab("Jaccard distance") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
          plot.margin=unit(c(0,0,0,0), unit="cm"), plot.title = element_text(size=25))+
    guides(fill=F)
  
  a<-make_anno(annotation_table, network_jaccard)
  plo<-a[[1]]
  legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.6,0.2,0.3)))
}
makePlot2<-function(annotation_table, network_jaccard){
  g1 <-ggplot(network_jaccard, aes(y = percMean, x = run)) + geom_boxplot() +
    ggtitle("Fraction of 'Gold Standard' Toydata Genes Found in Solution") + ylab("Fraction") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), unit="cm"))+guides(fill=F)
  a<-make_anno(annotation_table, network_jaccard)
  plo<-a[[1]]
  legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.5,0.25,0.25)))
}
