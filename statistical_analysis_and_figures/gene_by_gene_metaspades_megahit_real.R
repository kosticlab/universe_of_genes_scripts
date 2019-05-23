###run correlations and plotting for gene by gene data
##20190221
#Tierney

library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)

setwd('~/Dropbox (HMS)/orfletons/revisions/megahit_metaspades/')

process_data<-function(df,title){
  toTest<-c("Gene.Length","Contig.Length","Avg_fold_gene","Avg_fold_contig")
  plotList<-list()
  for(t in toTest)
    local({
      t<-t
      p<-ggplot(df,aes(x=df[,t],fill=as.factor(df$X50_perc_status)))+geom_histogram(bins=250)+xlab(gsub('[._]+',' ',toupper(t)))+ylab('Frequency')+ ggtitle(title)+theme(legend.title = element_blank())
      plotList[[t]]<<-p
    })
  png(paste(gsub('.csv','',title),'_singleton_histograms_real_data.png',sep=''),width = 10000,height=5000)
  do.call("grid.arrange", c(plotList, ncol=3))
  dev.off()
}

ms<-read.csv('gene_by_gene_with_singleton_data_metaspades.csv',stringsAsFactors=FALSE)
ms$assembler_type<-rep('metaSPAdes',nrow(ms))
ms$X50_perc_status <- ifelse(ms$X50_perc_status >1 , 'Non-singleton','Singleton')
process_data(ms,'metaSPAdes')

mh<-read.csv('gene_by_gene_with_singleton_data_megahit_default.csv',stringsAsFactors=FALSE)
mh$assembler_type<-rep('megahit_default',nrow(mh))
mh$X50_perc_status <- ifelse(mh$X50_perc_status >1 , 'Non-singleton','Singleton')
process_data(mh,'MEGAHIT Default')

sfp_df<-rbind(ms,mh)
sfp_df$X50_perc_status <- ifelse(sfp_df$X50_perc_status=="Non-singleton", 0, 1)

#sfp_df<-sfp_df[sfp_df$Avg_fold_gene<100,]
#sfp_df<-sfp_df[sfp_df$Avg_fold_contig<100,]

fit1 <- glm(X50_perc_status ~ as.factor(assembler_type),data=sfp_df,family='binomial')
fe1<-as.data.frame(summary(fit1)$coefficients)
fe1$odds_ratios<-exp(coef(fit1))

#fit2 <- glm(X50_perc_status ~ I(as.numeric(Contig.Length)/sd(Contig.Length)) +  I(as.numeric(Avg_fold_contig)/sd(Avg_fold_contig))  + assembler_type,data=sfp_df,family='binomial')
#fe2<-as.data.frame(summary(fit2)$coefficients)
#fe2$odds_ratios<-exp(coef(fit2))

write.csv(fe1,'gene_by_gene_megahit_metaspades_regression_output_genes_singletons_nocutoff.csv')
#write.csv(fe2,'gene_by_gene_megahit_metaspades_regression_output_contigs_singletons_nocutoff.csv')



