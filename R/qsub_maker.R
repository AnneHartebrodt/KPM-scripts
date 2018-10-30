require(data.table)

base.dir<-"/home/h/hartebrodt/grid_job_KPM/run"
method<-"normDegSum"
network<-"biogrid"
#network<-"StringDB"
network.filename<-"Homo_sapiens_biogrid.tsv"
#network.filename<-"Homo_Sapiens_String.tsv"
fdr_cutoff <- 0.05
min_net<-0
max_net<-150
max_net_size=FALSE
perc_permute<-25
random_meth<-"degreeawarenodeswap"
ranking_method<-"median"
high_degree_node_permutation<-1

grid_counter<-1

for(background in c("greedy")) {
  for(termination in c("sliding_window", "largest", "maximum_distance")){
    for(method in c("sum", "normDegSum","meanLog")){
      for(def_high_deg in c(250)){
        for(high_deg in c(20)){
          for(random_meth in c("edgerewire")){
            for(seed in c(8998)){
              for(perc_permute in c(30)){
                for(fdr_cutoff in c(0.01)){
                  
                  name<-paste0(network,"_",method,"_",fdr_cutoff,"_",perc_permute, "_", seed,
                               "_", random_meth, "_", ranking_method, "_", high_deg, "_", 
                               def_high_deg, "_", background, "_", termination)
                  #m1<-0.0001
                  #m2<-0.2
                  #name<-paste0(m1, "_", m2)
                  file.dir <- as.character(file.path(base.dir,"sample_data", network))
                  result.dir<-as.character(file.path(base.dir,"result_true",network, name))
                  call.dir<-as.character(file.path(base.dir, "calls" ,network, name))
                  grid.dir<-as.character(file.path(base.dir, "calls_console", network))
                  graph.file<-paste0(base.dir, "../data/networks/", network,"/", network.filename)
                  
                  
                  dir.create(result.dir, recursive = T)
                  dir.create(call.dir, recursive = T)
                  dir.create(grid.dir, recursive = T)
                  file.names<-dir(file.dir)
                  #file.names<-file.names[grep("86_30_72",file.names)]
                  general.files<-grep("general", file.names, value = T)
                  data<-""
                  for(filename_general in general.files){
                    
                    filename_general<-gsub(".tsv", "", filename_general)
                    filename_samplewise<-gsub("general_", "samplewise_", filename_general)
                    
                    if(max_net_size){
                      exp<-gsub("general_", "expression_type_", filename_general)
                      data<-fread(file.path(file.dir, paste0(exp, ".tsv")))
                      max_net<-length(which(data$V2=="differential"))
                      slack<-length(gregexpr("_", exp)[[1]])-2
                      max_net<-max_net+slack*floor(max_net/10)
                    }
                    
                    diagnostic_call <- paste(
                      "-numProc=1",
                      paste0("-matrix1=",
                             file.path(file.dir,
                                       paste0(filename_samplewise, ".tsv"))),
                      "-datasetsFile=/home/anne/Documents/Master/MA/Testing/datasets.txt",
                      "-summaryFile=harkan1.txt -combineOp=OR -pathwaysStatsFile=harkan2.txt",
                      paste0(
                        "-resultsDir=", file.path(result.dir,
                                                  gsub("general", "result", filename_general))
                      ),
                      paste0("-graphFile=",graph.file),
                      "-geneStatsFile=harkan3.txt -K=5 -L1=2 -maxsolutions=5",
                      "-algo=FDR -strategy=FDR -Umove_bens -comparator=LET",
                      "-significance_level=0.01",
                      paste0("-L1_pvalues=", file.path(file.dir,
                                                       filename_general), ".tsv"),
                      "-use_double -mfHeader",
                      "-validation_file=/home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/COAD-VAL-ENTREZ.txt",
                      "-L1_pvaluecutoff=0.01",
                      paste0("-aggregation_method=", method),
                      paste0("-fdr_cutoff=", fdr_cutoff),
                      paste0("-min_network_size=", min_net),
                      paste0("-max_network_size=", max_net),
                      paste0("-perc_perturbation=", perc_permute),
                      paste0("-seed=", seed),
                      paste0("-perturbation_technique=", random_meth),
                      paste0("-ranking_method=",ranking_method),
                      paste0("-nr_high_degree_nodes=", high_deg),
                      paste0("-high_degree_nodes=", def_high_deg),
                      paste0("-termination_criterion=", termination),
                      paste0("-background=", background)
                      
                    )
                    
                    ffile<-file.path(call.dir,
                                     gsub("general", "call", filename_general))
                    write(diagnostic_call, file = ffile)
                    
                    diagnostic_call_console<-paste0("java -jar /home/h/hartebrodt/grid_job_KPM/bin/KPM-5-jar-with-dependencies.jar ~/test.txt ", 
                                                    ffile,
                                                    " ~/test.txt")
                    write(diagnostic_call_console, paste0(grid.dir, "/call.", grid_counter))
                    grid_counter<-grid_counter+1
                  }}}}}}}}}}
