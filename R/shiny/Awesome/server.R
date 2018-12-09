#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
require(data.table)
require(ggplot2)
require(gridExtra)
require(viridis)
require(cowplot)
require(stringi)
require(RColorBrewer)
require(plotly)
require(colorspace)

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
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,0), unit="cm"))+guides(fill=F)
  a<-make_anno(annotation_table, network_jaccard)
 plo<-a[[1]]
 legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.5,0.25,0.25)))
}
makeTable<-function(annotation_table, network_jaccard){
  table_new<-network_jaccard[,.SD[which.min(V3)], by=c("run", "filename")]
  table_new$id<-"min_V3"
  table<-rbind(table_new)
  table_new<-network_jaccard[,.SD[which.max(nrnodes)], by=c("run", "filename")]
  table_new$id<-"max_nodes"
  table<-rbind(table, table_new)
  table_new<-network_jaccard[,.SD[which.max(nrinteractions)], by=c("run", "filename")]
  table_new$id<-"max_inter"
  table<-rbind(table, table_new)
  table_new<-network_jaccard[,.SD[which.min(norm)], by=c("run", "filename")]
  table_new$id<-"min_norm"
  table<-rbind(table, table_new)
  return(table)
}

makeBest<-function(annotation_table, network_jaccard){
  table<-makeTable(annotation_table, network_jaccard)
  g1<-ggplot(table, aes(y=jaccard_dists, x=run, fill=id))+geom_boxplot()+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          legend.position = "bottom",plot.margin=unit(c(0,0,0,0), unit="cm"))
  a<-make_anno(annotation_table, network_jaccard)
  
  plo<-a[[1]]
  legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.5,0.25,0.25))) 
}

weirdDots<-function(annotation_table, network_jaccard){
  tab<-makeTable(annotation_table, network_jaccard)
  min<-tab[,.SD[which.min(jaccard_dists)], by=c("run", "id")]
  g1<-ggplot(min, aes(x=run, y=jaccard_dists))+geom_point()+facet_wrap(~id)
  a<-make_anno(annotation_table, network_jaccard)
  plo<-a[[1]]
  legend<-a[[2]]
  return(plot_grid(g1,plo, legend, ncol = 1, align = "v", rel_heights = c(0.5,0.25,0.25))) 
}

select<-function(fdr, list_params){
  fdr<-fdr[Network %in% list_params$Network & 
             Aggr_meth %in% list_params$Aggr_meth &
             FDR %in% list_params$FDR &
             Perc_perm %in% list_params$Perc_perm &
             Seed %in% list_params$Seed &
             Perm_meth %in% list_params$Perm_meth &
             Sele_meth %in% list_params$Sele_meth &
             Perm_hi_de %in% list_params$Perm_hi_de &
             Deg_hi %in% list_params$Deg_hi &
             background %in% list_params$background &
             break1 %in% list_params$break1 &
             window %in% list_params$window
           ]
  return(fdr)
}
plotThresholds<-function(fdr){
  fdr$V1<-as.factor(fdr$V1)
  g1<-ggplot(fdr, aes(x=V1,y=V2))+geom_boxplot()+facet_wrap(~V4)
  return(g1)
}
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=1000*1024^2)
  network_jaccard<-NULL
  annotation_table<-NULL
  list_params<-NULL
  
  observeEvent(input$jaccard_table, {
    network_jaccard<<-fread(input$jaccard_table$datapath, sep="\t")
    #network_jaccard<<-network_jaccard[run!="run"]
    })

  
  observeEvent(input$annotation_table,{
    annotation_table<<-fread(input$annotation_table$datapath, sep="\t")
    annotation_table<<-annotation_table[Network != "Network"]
    list_params<-sapply(colnames(annotation_table), function(x) unique(annotation_table[,x, with=FALSE]))
    names(list_params)<-colnames(annotation_table)
    list_params<-list_params[-which(names(list_params)=="orig")]
    list_params<<-list_params
    checks<-list()
    for( i in 1:length(list_params)) {
      var<-names(i)
      checks[[i]]<-checkboxGroupInput(names(list_params[i]), names(list_params[i]), choiceNames=as.character(list_params[[i]]),
                                      choiceValues=as.character(list_params[[i]]), selected = as.character(list_params[[i]]))
    }
    output$check<-renderUI(checks)
  }
  )
  observeEvent(input$fdr, {
    fdr<<-fread(input$fdr$datapath) 
  })
  
  output$value <- renderPrint({ input$directory
    })
  
  output$select <- renderPrint({ input$checkGroup })
  
  observeEvent(input$update, {
    listig<-list(input$Network, input$Aggr_meth,
                 input$FDR, input$Perc_perm, input$Seed, 
                 input$Perm_meth, input$Sele_meth, 
                 input$Perm_hi_de, input$Deg_hi,
                 input$background, input$break1,input$window)
    names(listig)<-names(list_params)
    print(listig)
    net<-network_jaccard
    net<-merge(net, annotation_table, by.x="run", by.y="orig", allow.cartesian=T)
    net<-select(net, listig)
    print(annotation_table)
    fdr_sub<-merge(fdr, annotation_table, by.x="V4", by.y="orig", allow.cartesian=T)
    print(fdr_sub)
    fdr_sub<-select(fdr_sub, listig)
    
    output$eval<-renderPlot({
      makePlot(annotation_table, net)
    })
    
    output$best<-renderPlot({
      makeBest(annotation_table, net)
    })
    
    output$weirdDots<-renderPlot({
      weirdDots(annotation_table, net)
    })
    
    output$thresholds<-renderPlot({
      validate(need(length(unique(fdr_sub$V4))<4, "Please select fewer parameters"))
      print("plooting")
      
      plot(plotThresholds(fdr_sub))
      print("fin")
    })
  })
  

  observeEvent(input$action, {
    print(str(annotation_table))
    print(str(network_jaccard))
  })
  })







