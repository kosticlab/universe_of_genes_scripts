###run correlations and plotting for gene by gene data
##20190221
#Tierney

library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggridges)

process_data<-function(df,assembler){
  toTest<-c("Gene.Length","Contig.Length","Avg_ffe1old_gene","Avg_fold_contig")
  plotList<-list()
  df$False.Positive..method.1.<-as.factor(ifelse(as.numeric(df$False.Positive..method.1.) > 0,'False Positive','True Positive'))
  df$False.Positive..method.1. <- factor(df$False.Positive..method.1., levels = c("True Positive","False Positive"))
  for(t in toTest)
    local({
      t<-t
      p<-ggplot(df,aes(x=as.numeric(df[,t]),fill=df$False.Positive..method.1.))+geom_histogram(bins=250)+xlab(gsub('[._]+',' ',toupper(t)))+ylab('Frequency')+ theme(legend.title = element_blank())
      plotList[[t]]<<-p
    })
  pdf(paste(gsub('.csv','',assembler),'_histograms.pdf',sep=''),width = 20, height = 15)
  do.call("grid.arrange", c(plotList, ncol=2))
  dev.off()
}

setwd('~/Dropbox (HMS)/orfletons/revisions/synthetic_data/')

a<-read.csv('coverage_0_1_output/gene_by_gene_metaspades.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_metaspades.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_metaspades.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
ms<-rbind(a,b,c)
ms$assembler_type<-rep('metaSPAdes',nrow(ms))
#process_data(ms,'metaSPAdes')

a<-read.csv('coverage_0_1_output/gene_by_gene_megahit_default.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_megahit_default.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_megahit_default.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
mh_df<-rbind(a,b,c)
mh_df$assembler_type<-rep('MEGAHIT default',nrow(mh_df))
#process_data(mh_df,'megahit_default')

a<-read.csv('coverage_0_1_output/gene_by_gene_megahit_large.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_megahit_large.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_megahit_large.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
mh_lrg<-rbind(a,b,c)
mh_lrg$assembler_type<-rep('MEGAHIT large',nrow(mh_lrg))
#process_data(mh_lrg,'megahit_large')


a<-read.csv('coverage_0_1_output/gene_by_gene_megahit_sensitive.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_megahit_sensitive.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_megahit_sensitive.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
mh_sens<-rbind(a,b,c)
mh_sens$assembler_type<-rep('MEGAHIT sensitive',nrow(mh_sens))
#process_data(mh_sens,'megahit_sensitive')

geneDF<-rbind(ms,mh_df,mh_lrg,mh_sens)

geneDF$Avg_fold_contig<-as.numeric(geneDF$Avg_fold_contig)
geneDF$Avg_fold_gene<-as.numeric(geneDF$Avg_fold_gene)
geneDF$Gene.Length<-as.numeric(geneDF$Gene.Length)
geneDF$Contig.Length<-as.numeric(geneDF$Contig.Length)

geneDF<-geneDF[!is.na(geneDF$Avg_fold_contig),]

geneDF$group<-paste(geneDF$coverage,geneDF$assembler_type,sep=', ')
geneDF$group<-as.factor(geneDF$group)
geneDF$group<-factor(geneDF$group,levels(geneDF$group)[c(9:12,1:4,5:8)])#,c('0x-1x coverage, metaSPAdes','0x-1x coverage, MEGAHIT large','0x-1x coverage, MEGAHIT sensitive','0x-1x coverage, MEGAHIT default','0x-10x coverage, metaSPAdes','0x-10x coverage, MEGAHIT large','0x-10x coverage, MEGAHIT sensitive','0x-10x coverage, MEGAHIT default','10x-20x coverage, metaSPAdes','10x-20x coverage, MEGAHIT large','10x-20x coverage, MEGAHIT sensitive','10x-20x coverage, MEGAHIT default'))
geneDF$coverage <- relevel(as.factor(geneDF$coverage), ref = 2)

geneDF$False.Positive..method.1.<-as.factor(ifelse(as.numeric(geneDF$False.Positive..method.1.) > 0,'False Positive','True Positive'))
geneDF$False.Positive..method.1. <- factor(geneDF$False.Positive..method.1., levels = c("False Positive","True Positive"))

pdf('synthetic_false_positive_histograms_gene_length.pdf',width = 10, height = 20)
ggplot(geneDF, aes(Gene.Length,group,fill=geneDF$False.Positive..method.1.)) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Gene Length')
dev.off()

pdf('synthetic_false_positive_histograms_contig_length.pdf',width = 10, height = 20)
ggplot(geneDF, aes(Contig.Length,group,fill=geneDF$False.Positive..method.1.)) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Contig Length')
dev.off()

pdf('synthetic_false_positive_histograms_gene_coverage.pdf',width = 10, height = 20)
ggplot(geneDF, aes(log(Avg_fold_gene+1),group,fill=geneDF$False.Positive..method.1.)) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Gene Average Fold Coverage')
dev.off()

pdf('synthetic_false_positive_histograms_contig_coverage.pdf',width = 10, height = 20)
ggplot(geneDF, aes(log(Avg_fold_contig+1),group,fill=geneDF$False.Positive..method.1.)) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Contig Average Fold Coverage')
dev.off()

geneDF<-rbind(ms,mh_df,mh_lrg,mh_sens)
geneDF$Avg_fold_contig<-as.numeric(geneDF$Avg_fold_contig)
geneDF$Avg_fold_gene<-as.numeric(geneDF$Avg_fold_gene)
geneDF$Gene.Length<-as.numeric(geneDF$Gene.Length)
geneDF$Contig.Length<-as.numeric(geneDF$Contig.Length)

geneDF<-geneDF[!is.na(geneDF$Avg_fold_contig),]

geneDF$group<-paste(geneDF$coverage,geneDF$assembler_type,sep=', ')
geneDF$group<-as.factor(geneDF$group)
geneDF$group<-factor(geneDF$group,levels(geneDF$group)[c(9:12,1:4,5:8)])#,c('0x-1x coverage, metaSPAdes','0x-1x coverage, MEGAHIT large','0x-1x coverage, MEGAHIT sensitive','0x-1x coverage, MEGAHIT default','0x-10x coverage, metaSPAdes','0x-10x coverage, MEGAHIT large','0x-10x coverage, MEGAHIT sensitive','0x-10x coverage, MEGAHIT default','10x-20x coverage, metaSPAdes','10x-20x coverage, MEGAHIT large','10x-20x coverage, MEGAHIT sensitive','10x-20x coverage, MEGAHIT default'))
geneDF$coverage <- relevel(as.factor(geneDF$coverage), ref = 2)

#geneDF<-geneDF[geneDF$Avg_fold_gene<100,]
#geneDF<-geneDF[geneDF$Avg_fold_contig<100,]

fit1 <- glm(as.factor(False.Positive..method.1.) ~  I(Gene.Length/sd(Gene.Length)) + I(as.numeric(Avg_fold_gene)/sd(Avg_fold_gene))+ as.factor(assembler_type) + as.factor(coverage),data=geneDF,family="binomial")
fe1<-as.data.frame(summary(fit1)$coefficients)
fe1$odds_ratios<-exp(coef(fit1))

fit2 <- glm(as.factor(False.Positive..method.1.) ~ I(Contig.Length/sd(Contig.Length)) + I(as.numeric(Avg_fold_contig)/sd(Avg_fold_contig))+ as.factor(assembler_type) + as.factor(coverage),data=geneDF,family='binomial')
fe2<-as.data.frame(summary(fit2)$coefficients)
fe2$odds_ratios<-exp(coef(fit2))

write.csv(fe1,'gene_by_gene_regression_output_genes_100.csv')
write.csv(fe2,'gene_by_gene_regression_output_contigs.csv')

