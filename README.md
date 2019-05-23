# Scripts for the "Universe of Genes" project

Code for running the analysis in our universe of genes project. Relevant scripts for different analytical steps are referenced in methods section of publication. Each directory has an individual README describing the function of each script.

## Directories

### contig_analysis

Contains scripts for extracting contigs from raw PROKKA output, binning them, and grouping them by if they contain only singletons, non-singletons, or a mixture of the two.

### database_construction

Tool for generating csv file (to be ingested into database) from gene catalog and its associated annotations.

### extrapolation

Scripts for extrapolating total gene content of the microbiome.

### gene_catalog_construction

Pipeline for assembling metagenomes, predicting genes, and constructing and processing gene catalog.

### statistical_analysis_and_figures

Statistical analyses and R-derived figures used in the paper.

### synthetic_data_benchmarking

Pipeline for generating and processing (assembling, etc) synthetic metagenomes.

## Project website:

https://microbial-genes.bio/

## Project abstract: 

Despite substantial interest in the species diversity of the human microbiome and its role in human disease, we have not quantitated the scale of human microbiome genetic diversity, an instrumental task for understanding human-microbe interactions. Here, to do so, we conducted a cross-study meta-analysis of metagenomes from two niches, the mouth and gut, amassing 3,655 samples from 13 studies. We found staggering genetic heterogeneity in our dataset, identifying, at the 95% identity level, a total of 45,666,334 non-redundant genes (23,961,508 in the oral, 22,254,436 in the gut). We found that 50% of the genes in both datasets were “singletons”, meaning they were unique to a single metagenomic sample. We identified that singletons were enriched for discrete functions (compared to non-singletons) and arose from sub-population specific and extremely rare microbial strains. Overall, these results serve as a potential explanation for the large, unexplained heterogeneity observed in microbiome-derived human phenotypes. 



