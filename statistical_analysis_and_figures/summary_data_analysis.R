###run correlations on summary data
##20190221
#Tierney

library(ggplot2)
library(gridExtra)
library(reshape2)
library(broom)
library(ggpubr)
library(cowplot)

setwd('~/Dropbox (HMS)/orfletons/revisions/synthetic_data')
df<-read.csv('summary_data_coverage_analysis.csv')

covRanges<-unique(df$coverage_range)
covRanges<-covRanges[order(covRanges)]

data1<-df[df$coverage_range==covRanges[1],]
data2<-df[df$coverage_range==covRanges[2],]
data3<-df[df$coverage_range==covRanges[3],]

compute_summary_correlations<-function(data,covRange){
  #false positive rates
  msfp<-data$False.Positives..Compared.to.Ground.Truth..prokka_output_metaspades_contigs/data$Prokka.Gene.Count.prokka_output_metaspades_contigs
  mh_dffp<-data$False.Positives..Compared.to.Ground.Truth..prokka_output_megahit_contigs_default/data$Prokka.Gene.Count.prokka_output_megahit_contigs_default
  mh_sensfp<-data$False.Positives..Compared.to.Ground.Truth..prokka_output_megahit_contigs_sensitive/data$Prokka.Gene.Count.prokka_output_megahit_contigs_sensitive
  mh_largefp<-data$False.Positives..Compared.to.Ground.Truth..prokka_output_megahit_contigs_large/data$Prokka.Gene.Count.prokka_output_megahit_contigs_large
  
  fp_data<-data.frame(msfp,mh_dffp,mh_sensfp,mh_largefp)
  colnames(fp_data)<-c('MetaSPAdes','Megahit Default','Megahit Sensitive','Megahit Large')
  fp_data_long<-melt(fp_data)
  pdf(paste('false_positive_rate_boxplot_',covRange,'.pdf',sep=''))
  print(ggplot(fp_data_long,aes(y=value,x=variable))+geom_boxplot()+xlab('')+ylab('False Positive Rate')+ggtitle(paste('False Positive Rate by Assembly Conditions ','Coverage Range ',covRange,sep=''))+theme(plot.title = element_text(size = 8, face = "bold")))
  dev.off()
  fp_regression<-tidy(summary(lm(fp_data_long$value~fp_data_long$variable)))
  fp_regression$term<-c('INTERCEPT','Megahit Default','Megahit Sensitive','Megahit Long')
  write.csv(fp_regression,paste('false_positive_regression_by_variable_',covRange,'.csv',sep=''))

  ###correlations and plots
  #contig length
  a<-ggplot(data,aes(x=msfp,y=data$Mean.Contig.Length.prokka_output_metaspades_contigs)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Length')+ggtitle('False Positives vs. Mean Contig Length\n(MetaSPAdes)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  b<-ggplot(data,aes(x=mh_dffp,y=data$Mean.Contig.Length.prokka_output_megahit_contigs_default)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Length')+ggtitle('False Positives vs. Mean Contig Length\n(Megahit Default)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  c<-ggplot(data,aes(x=mh_sensfp,y=data$Mean.Contig.Length.prokka_output_megahit_contigs_sensitive)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Length')+ggtitle('False Positives vs. Mean Contig Length\n(Megahit Sensitive)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  d<-ggplot(data,aes(x=mh_largefp,y=data$Mean.Contig.Length.prokka_output_megahit_contigs_large)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Length')+ggtitle('False Positives vs. Mean Contig Length\n(Megahit Large)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  pdf(paste('contig_length_fp_cor_',covRange,'.pdf',sep=''))
  grid.arrange(a,b,c,d)
  dev.off()
  
  #gene length
  a<-ggplot(data,aes(x=msfp,y=data$Mean.Gene.Length.prokka_output_metaspades_contigs)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Length')+ggtitle('False Positives vs. Mean Gene Length\n(MetaSPAdes)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  b<-ggplot(data,aes(x=mh_dffp,y=data$Mean.Gene.Length.prokka_output_megahit_contigs_default)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Length')+ggtitle('False Positives vs. Mean Gene Length\n(Megahit Default)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  c<-ggplot(data,aes(x=mh_sensfp,y=data$Mean.Gene.Length.prokka_output_megahit_contigs_sensitive)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Length')+ggtitle('False Positives vs. Mean Gene Length\n(Megahit Sensitive)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  d<-ggplot(data,aes(x=mh_largefp,y=data$Mean.Gene.Length.prokka_output_megahit_contigs_large)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Length')+ggtitle('False Positives vs. Mean Gene Length\n(Megahit Large)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  pdf(paste('gene_length_fp_cor_',covRange,'.pdf',sep=''))
  grid.arrange(a,b,c,d)
  dev.off()
  
  #contig coverage (fold)
  a<-ggplot(data,aes(x=msfp,y=data$Contig.Average.Coverage.prokka_output_metaspades_contigs)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Coverage')+ggtitle('False Positives vs. Mean Contig Fold Coverage\n(MetaSPAdes)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  b<-ggplot(data,aes(x=mh_dffp,y=data$Contig.Average.Coverage.prokka_output_megahit_contigs_default)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Coverage')+ggtitle('False Positives vs. Mean Contig Fold Coverage\n(Megahit Default)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  c<-ggplot(data,aes(x=mh_sensfp,y=data$Contig.Average.Coverage.prokka_output_megahit_contigs_sensitive)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Coverage')+ggtitle('False Positives vs. Mean Contig Fold Coverage\n(Megahit Sensitive)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  d<-ggplot(data,aes(x=mh_largefp,y=data$Contig.Average.Coverage.prokka_output_megahit_contigs_large)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Contig Coverage')+ggtitle('False Positives vs. Mean Contig Fold Coverage\n(Megahit Large)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  pdf(paste('contig_cov_fp_cor_',covRange,'.pdf',sep=''))
  grid.arrange(a,b,c,d)
  dev.off()
  
  #gene coverage (fold)
  a<-ggplot(data,aes(x=msfp,y=data$Gene.Average.Coverage.prokka_output_metaspades_contigs)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Coverage')+ggtitle('False Positives vs. Mean Gene Fold Coverage\n(MetaSPAdes)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  b<-ggplot(data,aes(x=mh_dffp,y=data$Gene.Average.Coverage.prokka_output_megahit_contigs_default)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Coverage')+ggtitle('False Positives vs. Mean Gene Fold Coverage\n(Megahit Default)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  c<-ggplot(data,aes(x=mh_sensfp,y=data$Gene.Average.Coverage.prokka_output_megahit_contigs_sensitive)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Coverage')+ggtitle('False Positives vs. Mean Gene Fold Coverage\n(Megahit Sensitive)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  d<-ggplot(data,aes(x=mh_largefp,y=data$Gene.Average.Coverage.prokka_output_megahit_contigs_large)) + geom_point()+stat_cor(method = "pearson")+xlab('False Positive Rate')+ylab('Mean Gene Coverage')+ggtitle('False Positives vs. Mean Gene Fold Coverage\n(Megahit Large)')+ theme(plot.title = element_text(size = 8, face = "bold"))
  pdf(paste('gene_cov_fp_cor_',covRange,'.pdf',sep=''))
  grid.arrange(a,b,c,d)
  dev.off()
}

compute_summary_correlations(data1,covRanges[1])
compute_summary_correlations(data2,covRanges[2])
compute_summary_correlations(data3,covRanges[3])

###repeat for overall data
compute_summary_correlations(df,'0_30')


