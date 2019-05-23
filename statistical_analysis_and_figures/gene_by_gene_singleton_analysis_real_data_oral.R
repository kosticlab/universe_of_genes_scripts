
setwd('~/Dropbox (HMS)/orfletons/revisions/nton_status_by_length/')

d<-read.csv('nton_by_gene_contig_length_with_coverage_data.csv',stringsAsFactors=FALSE)

d$GENE.NTON.STATUS <- ifelse(d$GENE.NTON.STATUS=="Non-singleton", 0, 1)

d$READ.COUNT.F1<-as.numeric(d$READ.COUNT.F1)
d$READ.COUNT.F2<-as.numeric(d$READ.COUNT.F2)
d$READ.COUNT.F2[is.na(d$READ.COUNT.F2)] <- 0
d$READ.COUNT.F1[is.na(d$READ.COUNT.F1)] <- 0
d$TOTAL.READS <- d$READ.COUNT.F1+d$READ.COUNT.F2
d$TOTAL.READS<-as.numeric(d$TOTAL.READS)/sd(as.numeric(d$TOTAL.READS))
d$CONTIG.LENGTH<-as.numeric(d$CONTIG.LENGTH)/sd(as.numeric(d$CONTIG.LENGTH))
d$CONTIG.LENGTH[is.na(d$CONTIG.LENGTH)] <- 0
d$FOLD.COVERAGE.CONTIG<-as.numeric(as.character(d$FOLD.COVERAGE.CONTIG))/sd(as.numeric(as.character(d$FOLD.COVERAGE.CONTIG)))
d$GENE.LENGTH <- d$GENE.LENGTH/sd(d$GENE.LENGTH)
d$GENE.LENGTH[is.na(d$GENE.LENGTH)] <- 0

fit1 <- glm(as.factor(GENE.NTON.STATUS) ~ CONTIG.LENGTH + FOLD.COVERAGE.CONTIG + TOTAL.READS,data=d,family='binomial')
fe1 <- as.data.frame(summary(fit1)$coefficients)
fe1$odds_ratios<-exp(coef(fit1))

write.csv(fe1,'gene_by_gene_regression_output_contigs_singletons_real_data_oral.csv')

fit2 <- glm(as.factor(GENE.NTON.STATUS) ~ GENE.LENGTH + TOTAL.READS,data=d,family='binomial')
fe2<-as.data.frame(summary(fit2)$coefficients)
fe2$odds_ratios<-exp(coef(fit2))

write.csv(fe2,'gene_by_gene_regression_output_genes_singletons_real_data_oral.csv')
