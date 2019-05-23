library(readxl)
library(ggplot2)
library(cowplot)
library(ggpubr)

setwd('~/Dropbox (HMS)/orfletons/revisions/')
d=read_excel('sample_summary_data.xls')

pdf('read_count_v_singleton_count.pdf')
ggplot(d,aes(x=d$`Total Read Count`,y=d$`Singletons (50% ID)`))+geom_point()+xlab('Total Read Count')+ ylab('Singleton Count')+ggtitle('Read Count vs. Singleton Count')+stat_cor(method='spearman',label.x=250000000)
dev.off()

pdf('read_count_v_singleton_fraction.pdf')
ggplot(d,aes(x=d$`Total Read Count`,y=d$`Fraction of singleton-genes (50% ID)`))+geom_point()+xlab('Total Read Count')+ ylab('Fraction of genes that were singletons')+ggtitle('Read Count vs. Singleton Fraction')+stat_cor(method='spearman',label.x=250000000)
dev.off()

pdf('singleton_fraction_hist.pdf')
ggplot(d,aes(x=d$`Fraction of singleton-genes (50% ID)`))+geom_histogram(bins=500)+xlab('Singleton gene fraction per sample')+ ylab('Frequency')+ggtitle('Singleton Fraction Histogram')
dev.off()

pdf('singleton_count_hist.pdf')
ggplot(d,aes(x=d$`Singletons (50% ID)`))+geom_histogram(bins=500)+xlab('Singleton gene count per sample')+ ylab('Frequency')+ggtitle('Singleton Count Histogram')
dev.off()