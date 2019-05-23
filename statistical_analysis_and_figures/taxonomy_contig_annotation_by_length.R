#plot taxonomy contig relationship by different lengths

library(ggplot2)
library(ggpubr)
library(cowplot)

#get distribution by gene count and by length
setwd('~/Dropbox (HMS)/orfletons/revisions')

plotRel<-function(filename1,filename2,number,cutoff){
  contigtaxabin<-read.csv(filename1)
  colnames(contigtaxabin)<-c('Taxonomic_Annotation','Singleton_Counts','Non_Singleton_Counts')
  contigtaxabin$Normalized_Counts<-contigtaxabin$Non_Singleton_Counts/sum(contigtaxabin$Non_Singleton_Counts)
  contigtaxabin$Normalized_Singleton_Counts<-contigtaxabin$Singleton_Counts/sum(contigtaxabin$Singleton_Counts)
  contigtaxabin[,2]<-contigtaxabin[,2]+1
  contigtaxabin[,3]<-contigtaxabin[,3]+1
  
  contigtaxabin1<-read.csv(filename2)
  colnames(contigtaxabin1)<-c('Taxonomic_Annotation','Singleton_Counts','Non_Singleton_Counts')
  contigtaxabin1$Normalized_Counts<-contigtaxabin1$Non_Singleton_Counts/sum(contigtaxabin1$Non_Singleton_Counts)
  contigtaxabin1$Normalized_Singleton_Counts<-contigtaxabin1$Singleton_Counts/sum(contigtaxabin1$Singleton_Counts)
  contigtaxabin1[,2]<-contigtaxabin1[,2]+1
  contigtaxabin1[,3]<-contigtaxabin1[,3]+1
  ggplot() + geom_point(data=contigtaxabin, aes(x=log(Normalized_Counts), y=log(Normalized_Singleton_Counts),color='oral'))+ ylab('ln(Normalized non-singleton counts)')+xlab('ln(Normalized singleton counts)')+geom_point(data=contigtaxabin1, aes(x=log(Normalized_Counts), y=log(Normalized_Singleton_Counts),color='gut')) +labs(title=paste("Taxonomic bins of singleton and non-singleton contigs: ","min length ",number,"bp, min ",cutoff," per contig",sep='')) +geom_abline(slope=1, intercept=0)+ scale_x_continuous(limits=c(-20,0))+ scale_y_continuous(limits=c(-20,0))+theme(legend.position="bottom")+theme(legend.title=element_blank())
  ggsave(paste('./taxa_relationship_contig_length/taxa_relationship_combined_',number,'_',cutoff,'.pdf',sep=''),width = 10,height=10)
  
  foo<-rbind(contigtaxabin,contigtaxabin1)
  foo2=cor.test(foo$Normalized_Counts,foo$Normalized_Singleton_Counts)
  .GlobalEnv$corList[length(.GlobalEnv$corList)+1]<-as.numeric(foo2$estimate)
  .GlobalEnv$plist[length(.GlobalEnv$corList)+1]<-as.numeric(foo2$p.value)
  .GlobalEnv$rows[length(.GlobalEnv$corList)+1]<-paste("min length ",number,"bp, min ",cutoff," per contig",sep='')
  
}

corList<-list()
plist<-list()
rows<-list()


plotRel('gene_catalog_binned_contig_data_with_1_lengths_oral_combined_500.csv','gene_catalog_binned_contig_data_with_1_lengths_gut_combined_500.csv','500','2 genes')
plotRel('gene_catalog_binned_contig_data_with_1_lengths_oral_combined_1000.csv','gene_catalog_binned_contig_data_with_1_lengths_gut_combined_1000.csv','1000','2 genes')
plotRel('gene_catalog_binned_contig_data_with_1_lengths_oral_combined_1500.csv','gene_catalog_binned_contig_data_with_1_lengths_gut_combined_1500.csv','1500','2 genes')
plotRel('gene_catalog_binned_contig_data_with_1_lengths_oral_combined_2000.csv','gene_catalog_binned_contig_data_with_1_lengths_gut_combined_2000.csv','2000','2 genes')
plotRel('gene_catalog_binned_contig_data_with_2_lengths_oral_combined_500.csv','gene_catalog_binned_contig_data_with_2_lengths_gut_combined_500.csv','500','3 genes')
plotRel('gene_catalog_binned_contig_data_with_2_lengths_oral_combined_1000.csv','gene_catalog_binned_contig_data_with_2_lengths_gut_combined_1000.csv','1000','3 genes')
plotRel('gene_catalog_binned_contig_data_with_2_lengths_oral_combined_1500.csv','gene_catalog_binned_contig_data_with_2_lengths_gut_combined_1500.csv','1500','3 genes')
plotRel('gene_catalog_binned_contig_data_with_2_lengths_oral_combined_2000.csv','gene_catalog_binned_contig_data_with_2_lengths_gut_combined_2000.csv','2000','3 genes')

corOutput<-do.call(rbind, Map(data.frame, corellation=unlist(corList), pvalue=unlist(plist)))
row.names(corOutput)<-unlist(rows)
write.csv(corOutput,'./taxa_relationship_contig_length/correlation_output.csv')




