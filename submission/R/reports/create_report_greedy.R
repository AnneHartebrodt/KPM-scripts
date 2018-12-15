require(optparse)


option_list = list(
  make_option(c("-f", "--filelocation"), action="store", default=NA, type='character',
              help="location of FDR output"),
  make_option(c("-g", "--goldstandard"), action="store", default=NA, type='character',
              help="goldstandard file")
)
opt = parse_args(OptionParser(option_list=option_list))

opt$filelocation<-normalizePath(opt$filelocation)
opt$goldstandard<-normalizePath(opt$goldstandard)



rmarkdown::render("~/Masterarbeit/R/reports/report.Rmd", 
                  params = list(filelocation = opt$filelocation, goldstandard = opt$goldstandard),
                  output_file   =file.path(opt$filelocation,"report.html"))
