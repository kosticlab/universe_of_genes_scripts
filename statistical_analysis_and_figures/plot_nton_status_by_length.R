###plot nton status length data
##20190221
#Tierney

library(ggplot2)
library(cowplot)

setwd('~/Dropbox (HMS)/orfletons/revisions/nton_status_by_length')

plotdata_withcoverage<-function(d,name){
  #plot gene length distribution
  pdf('./nton_status_contig_coverage_oral_logged.pdf')
  ggplot(d,aes(x=FOLD.COVERAGE.CONTIG,fill=as.factor(GENE.NTON.STATUS)))+geom_histogram(stat='count',bins = 500)+scale_y_log10()+xlab('Contig coverage') + ylab('Frequency')+ggtitle(paste('Frequency of contig coverage by singleton status:',name))
  dev.off()
}

plotdata<-function(d,name){
  #plot gene length distribution
  pdf(paste('./nton_status_gene_length_',name,'.pdf',sep=''))
  ggplot(d,aes(x=as.numeric(as.character(GENE.LENGTH)),fill=as.factor(GENE.NTON.STATUS)))+geom_histogram(bins = 500)+xlab('Gene length')+scale_y_log10()+ylab('Frequency')+ggtitle(paste('Frequency of gene lengths by singleton status:',name))+theme(legend.position = "none") 
  dev.off()
  
  #plot contig length distribution
  pdf(paste('./nton_status_contig_length_',name,'.pdf',sep=''))
  ggplot(d,aes(x=as.numeric(as.character(CONTIG.LENGTH)),fill=as.factor(GENE.NTON.STATUS)))+geom_histogram(bins = 500)+xlab('Contig length')+scale_y_log10() +ggtitle(paste('Frequency of contig lengths by singleton status:',name))+ylab('Frequency')+theme(legend.position = "none") 
  dev.off()
}

d=read.csv('./nton_by_gene_contig_length_with_coverage_data.csv')
plotdata_withcoverage(d,'oral_coverage')
d=read.csv('nton_by_gene_contig_length_gut_oral.csv')
plotdata(d,'oral')
d=read.csv('./nton_by_gene_contig_length_gut.csv')
plotdata(d,'gut')
















