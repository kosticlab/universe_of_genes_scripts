###run correlations and plotting for gene by gene data
##20190221
#Tierney

library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)
library(ggridges)


setwd('~/Dropbox (HMS)/orfletons/revisions/synthetic_data/')

process_data<-function(df,title){
  toTest<-c("Gene.Length","Contig.Length","Avg_fold_gene","Avg_fold_contig")
  plotList<-list()
  for(t in toTest)
    local({
      t<-t
      p<-ggplot(df,aes(x=df[,t],fill=as.factor(df$X50_perc_status)))+geom_histogram(bins=250)+xlab(gsub('[._]+',' ',toupper(t)))+ylab('Frequency')+ theme(legend.title = element_blank())
      plotList[[t]]<<-p
    })
  png(paste(gsub('.csv','',title),'_singleton_histograms.png',sep=''),width = 5000,height=2500,res=600)
  do.call("grid.arrange", c(plotList, ncol=2))
  dev.off()
  pdf(paste(gsub('.csv','',title),'_singleton_histograms.pdf',sep=''),width = 20,height=15)
  do.call("grid.arrange", c(plotList, ncol=2))
  dev.off()
}

a<-read.csv('coverage_0_1_output/gene_by_gene_with_singleton_data_metaspades.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_with_singleton_data_metaspades.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_with_singleton_data_metaspades.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
ms<-rbind(a,b,c)
ms$assembler_type<-rep('metaSPAdes',nrow(ms))
ms$X50_perc_status <- ifelse(ms$X50_perc_status >1 , 'Non-singleton','Singleton')
#process_data(ms,'metaSPAdes')

a<-read.csv('coverage_0_1_output/gene_by_gene_with_singleton_data_megahit_default.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_with_singleton_data_megahit_default.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_with_singleton_data_megahit_default.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
mh_df<-rbind(a,b,c)
mh_df$X50_perc_status <- ifelse(mh_df$X50_perc_status >1 , 'Non-singleton','Singleton')
mh_df$assembler_type<-rep('MEGAHIT default',nrow(mh_df))
#process_data(mh_df,'megahit_default')

a<-read.csv('coverage_0_1_output/gene_by_gene_with_singleton_data_megahit_large.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_with_singleton_data_megahit_large.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_with_singleton_data_megahit_large.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
mh_lrg<-rbind(a,b,c)
mh_lrg$X50_perc_status <- ifelse(mh_lrg$X50_perc_status >1 , 'Non-singleton','Singleton')
mh_lrg$assembler_type<-rep('MEGAHIT large',nrow(mh_lrg))
#process_data(mh_lrg,'megahit_large')

a<-read.csv('coverage_0_1_output/gene_by_gene_with_singleton_data_megahit_sensitive.csv',stringsAsFactors=FALSE)
b<-read.csv('coverage_0_10_output/gene_by_gene_with_singleton_data_megahit_sensitive.csv',stringsAsFactors=FALSE)
c<-read.csv('coverage_10_20_output/gene_by_gene_with_singleton_data_megahit_sensitive.csv',stringsAsFactors=FALSE)
a$coverage<-rep('0x-1x coverage',nrow(a))
b$coverage<-rep('0x-10x coverage',nrow(b))
c$coverage<-rep('10x-20x coverage',nrow(c))
mh_sens<-rbind(a,b,c)
mh_sens$X50_perc_status <- ifelse(mh_sens$X50_perc_status >1 , 'Non-singleton','Singleton')
mh_sens$assembler_type<-rep('MEGAHIT sensitive',nrow(mh_sens))
#process_data(mh_sens,'megahit_sensitive')

sfp_df<-rbind(ms,mh_df,mh_lrg,mh_sens)

sfp_df$Gene.Length<-as.numeric(sfp_df$Gene.Length)
sfp_df<-sfp_df[!is.na(sfp_df$Gene.Length),]


sfp_df$group<-paste(sfp_df$coverage,sfp_df$assembler_type,sep=', ')
sfp_df$group<-as.factor(sfp_df$group)
sfp_df$group<-factor(sfp_df$group,levels(sfp_df$group)[c(9:12,1:4,5:8)])#,c('0x-1x coverage, metaSPAdes','0x-1x coverage, MEGAHIT large','0x-1x coverage, MEGAHIT sensitive','0x-1x coverage, MEGAHIT default','0x-10x coverage, metaSPAdes','0x-10x coverage, MEGAHIT large','0x-10x coverage, MEGAHIT sensitive','0x-10x coverage, MEGAHIT default','10x-20x coverage, metaSPAdes','10x-20x coverage, MEGAHIT large','10x-20x coverage, MEGAHIT sensitive','10x-20x coverage, MEGAHIT default'))
sfp_df$coverage <- relevel(as.factor(sfp_df$coverage), ref = 2)

pdf('synthetic_singleton_histograms_gene_length.pdf',width = 10, height = 20)
ggplot(sfp_df, aes(Gene.Length,group,fill=as.factor(sfp_df$X50_perc_status))) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Gene Length')
dev.off()

pdf('synthetic_singleton_histograms_contig_length.pdf',width = 10, height = 20)
ggplot(sfp_df, aes(Contig.Length,group,fill=as.factor(sfp_df$X50_perc_status))) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Contig Length')
dev.off()

pdf('synthetic_singleton_histograms_gene_coverage.pdf',width = 10, height = 20)
ggplot(sfp_df, aes(log(Avg_fold_gene+1),group,fill=as.factor(sfp_df$X50_perc_status))) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Gene Average Fold Coverage')
dev.off()

pdf('synthetic_singleton_histograms_contig_coverage.pdf',width = 10, height = 20)
ggplot(sfp_df, aes(log(Avg_fold_contig+1),group,fill=as.factor(sfp_df$X50_perc_status))) + stat_density_ridges(alpha=.5)+ theme(legend.title = element_blank())+ylab('')+xlab('Contig Average Fold Coverage')
dev.off()

#sfp_df<-sfp_df[sfp_df$Avg_fold_gene<100,]
#sfp_df<-sfp_df[sfp_df$Avg_fold_contig<100,]

sfp_df$X50_perc_status <- ifelse(sfp_df$X50_perc_status=="Non-singleton", 0, 1)

fit1 <- glm(as.factor(X50_perc_status) ~ as.factor(False.Positive..method.1.) + I(Gene.Length/sd(Gene.Length)) + I(as.numeric(Avg_fold_gene)/sd(Avg_fold_gene)) + as.factor(assembler_type)+ as.factor(coverage),data=sfp_df,family='binomial')
fe1<-as.data.frame(summary(fit1)$coefficients)
fe1$odds_ratios<-exp(coef(fit1))

#fit1 <- glm(X50_perc_status ~ False.Positive..method.1.,data=sfp_df,family='binomial')
#fe1<-as.data.frame(summary(fit1)$coefficients)
#fe1$odds_ratios<-exp(coef(fit1))

fit2 <- glm(X50_perc_status ~ as.factor(False.Positive..method.1.)+ I(as.numeric(Contig.Length)/sd(Contig.Length)) + I(as.numeric(Avg_fold_contig)/sd(Avg_fold_contig)) + as.factor(assembler_type) + as.factor(coverage),data=sfp_df,family='binomial')
fe2<-as.data.frame(summary(fit2)$coefficients)
fe2$odds_ratios<-exp(coef(fit2))

write.csv(fe1,'gene_by_gene_regression_output_genes_singletons_withcutoff.csv')
write.csv(fe2,'gene_by_gene_regression_output_contigs_singletons_withcutoff.csv')



