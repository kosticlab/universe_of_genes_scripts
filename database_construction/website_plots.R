
library(ggplot2)
library(cowplot)
library(hist)

setwd('~/Dropbox (HMS)/orfletons/revisions/db_figs/')

lengths=read.csv('gene_lengths_dist_oral_gut.csv')

jpeg('length_plot.jpg',width=5000,height=5000,res = 500)
ggplot(NULL) + geom_bar(alpha=.7,stat='identity',data=lengths,aes(y=as.numeric(lengths$X0),x=as.numeric(rownames(lengths)),fill='Gut'))+ geom_bar(alpha=.3,stat='identity',data=lengths,aes(y=as.numeric(lengths$X0.1),x=as.numeric(rownames(lengths)),fill='Oral')) +xlab('Gene length')+ylab('Frequency')+theme(legend.title=element_blank())
dev.off()

#taxa=read.csv('taxa_id_dist_oral_gut.csv')
#cols=colSums(taxa,na.rm = T)
#taxa2=taxa[2:3]/cols[2:3]
#taxa3=taxa2[order(taxa2[,1],decreasing = T)[1:25],]


#taxa4=taxa2[order(taxa2[,2],decreasing = T)[1:25],]


#jpeg('length_plot.jpg',width=5000,height=5000,res = 500)
#ggplot(NULL) + geom_bar(alpha=.7,stat='identity',data=lengths,aes(y=as.numeric(lengths$X0),x=as.numeric(rownames(lengths)),fill='Gut'))+ geom_bar(alpha=.3,stat='identity',data=lengths,aes(y=as.numeric(lengths$X0.1),x=as.numeric(rownames(lengths)),fill='Oral')) +xlab('Gene length')+ylab('Frequency')+theme(legend.title=element_blank())
#dev.off()