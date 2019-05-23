## Directory Overview
This directory contains the R scripts used to run statistical analyses and generate all the figures in the paper.

## Directory Contents
#### gene_by_gene_analysis.R 
Computes assocation between false positive genes and various covariates like gene/contig length/coverage and assembler type on synthetic data.

#### gene_by_gene_metaspades_megahit_real.R 
Computes assocation between singleton genes and assembler type on REAL assembled metagenomes.

#### gene_by_gene_singleton_analysis.R 	
Computes assocation between singleton genes and various covariates like gene/contig length/coverage and assembler type on synthetic data.

#### gene_by_gene_singleton_analysis_real_data_oral.R 	
Computes association between singleton genes and various covariates like gene/contig length/coverage and assembler type on real data.

#### orfleton_figures_both.Rmd 	
Generates majority of plots used in publication.

#### plot_nton_status_by_length.R 	
Plots n-ton status (if a gene is a singleton or non-singleton) by gene/contig length/coverage for real data.

#### plot_read_count_singleton_data.R 
Plots the association between depth of sequencing and number of singletons per sample on real data.

#### summary_data_analysis.R 	
Computes assocation between false positive genes and assembler type/parameters.

#### taxonomy_contig_annotation_by_length.R
Plots number of singleton/non-singleton contigs with different taxonomic annotations for contigs of different lengths.
