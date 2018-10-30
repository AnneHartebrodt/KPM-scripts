require(DESeq2)
require(illuminaio)
require(data.table)
require(limma)
require(affy)
require(genefilter)
require(ggplot2)

d<-fread("~/Documents/Master/MA/data/HDarray/GSE3790-GPL96_series_matrix.txt", skip = 50)
dat<-fread("~/Documents/Master/MA/data/HDarray/GSE3790-GPL96_series_matrix.txt", skip = "!series_matrix_table_begin")
genes<-fread("~/Documents/Master/MA/data/HDarray/GSE3790_family.soft", skip = "!platform_table_begin")
d_names<-as.character(d[1,2:ncol(d)])
stopifnot(all(d_names==colnames(dat)[2:ncol(dat)]))
colnames(dat)[2:ncol(dat)]<-as.character(d[1,2:ncol(d)])

# only CN samples
ind<-which(grepl("CN", colnames(d)))
values <-dat[,..ind]
  
annot<-d[,..ind]


control<-grepl("control", annot)
control_tab<-values[,..control]
hdi<-!control
hd<-values[,..hdi]

rowttests<-function(x, fac){
  #t<-t.test(x[1,1],x[1,2:ncol(x)], var.equal = T)$p.value
  l<- sapply(1:nrow(x), function(a) t.test(x[a,1],x[a,2:ncol(x)], var.equal = T)$p.value)
  return(l)
}


ttests<-list()
for (i in colnames(hd)){
  m<-as.matrix(cbind(hd[,..i], control_tab[, 2:ncol(control_tab)]))
  ttests[[i]]<-rowttests(m, fac=as.factor(c(1,rep(2,31))))
}

ttests<-as.data.frame(ttests)
write.table(ttests, file= "~/Documents/Master/MA/data/HDarray/processed/pvalue_diff_exp_raw.txt")

ttestsN<-p.adjust(unlist(ttests), method = "fdr")
testN<-data.table(nam = names(ttestsN),pval=ttestsN)
#testN<-dcast(testN,... ~nam)

write.table(ttestsN, file= "~/Documents/Master/MA/data/HDarray/processed/pvalue_diff_exp_norm.txt")

ttest_list<-data.table(pval=unlist(ttests))
gp<-ggplot(ttest_list, aes(x=pval) )+geom_histogram()+ggtitle("Histogram of p-values for 1 patient vs. all controls using a two-sided t-test")
plot(gp)

gpN<-ggplot(data.frame(pval=ttestsN), aes(x=pval))+geom_histogram()+ggtitle("Histogram of fdr corrected p-values for 1 patient vs. all controls using a two-sided t-test")
plot(gpN)

topL<-list()
topLUn<-list()
for (i in colnames(hd)){
  eset<-NULL
  m<-as.matrix(cbind(hd[,..i], control_tab))
  eset<-ExpressionSet(assayData=as.matrix(m))  
  design<-c(1, rep(2, ncol(control_tab)))
  design <- model.matrix(~ 0+factor(design))
  colnames(design) <- c("group1", "group2")
  fit <- lmFit(eset, design)
  contrast.matrix <- makeContrasts(group2-group1, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  top<-topTable(fit2, number = 22283,coef=1, adjust="BH")
  topLUn[[i]]<-top$P.Value
  topL[[i]]<-top$adj.P.Val
}

topL<-as.data.table(topL)

topL_pList<-data.frame(pval=unlist(topL))
pp<-ggplot(topL_pList, aes(x=pval))+geom_histogram()
plot(pp)

topLUn<-data.table(pval=unlist(topLUn))
pun<-ggplot(topLUn, aes(x=pval))+geom_histogram(aes(y = ..density..))+geom_density(aes(color="red"))
plot(pun)

ggplot()+geom_line(dexp(0:1, rate   =111 )+22600)

pun$data

build<-ggplot_build(pun)

d<-fitdistr(build$data[[1]]$density[1:25]-0.85, "Poisson")
mean(build$data[[1]]$density[15:25])

curve(dexp(x, rate   =1/(d$estimate*0.5) ))

map<-genes[, c("ID", "ENTREZ_GENE_ID")]

