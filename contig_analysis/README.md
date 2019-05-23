## Directory Overview

This directory contains probably the most complicated portions of our analysis. In short, we did the following:

1) Parse the gff files containing information on predicted genes and the contigs they arose from for each of our samples. Then constructe a tab-separated file where each line corresponded to a different contig in our dataset. Each line contained information on 1) what genes/consensus were on the contig 2) if those genes were singletons or non-singletons 3) the functional/taxonomic annotations of those genes. Each of these groups of information (i.e. species for each gene) are separated by pipes (|). For example:


Contig ID | Consensus Genes | Raw Genes | Nton-status | Singletons IDs | Non-singleton IDs | Species IDs | Functions  
------------ | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
contig_id | consensus_gene1\|consensus_gene2 | raw_gene1\|raw_gene2 | singleton\|non-singleton | raw_gene1 | raw_gene2 | taxa1\|taxa2 | fxn1\|fxn2

2) Separate contigs into those that were mixtures of singletons and non-singletons or exclusively one or the other. (split_contigdb_by_nton_status.py)
3) Bin contigs according to the predominant taxonomy of their component genes. (bin_contigs_species.py)

