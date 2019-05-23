## Directory overview

The folder contains the scripts we used to generate and run initial processing on our gene catalogs. Most of the raw analysis for our project – save the contig-level work – can be found in this directory.

## Directory contents

#### full_contig_parsing_and_singleton_hunting_pipeline.py
Master script that runs pipeline to identify singletons, build gene catalog, filter short genes, taxonomically annotate genes, run iterative cdhit clustering (at progressively lower percent identities), build binary gene matrices, compute sorenseon dice index.

#### binary_abundance_matrix.py	
Generates a binary matrix from CD-HIT output data, where the input is a CD-HIT .clstr file and the output is an nXm matrix, where n is the number of consensus genes and m is the number of samples.

#### parse_minpath_output_to_cdhit_format.py	
Converts minpath output to CD-HIT .clstr format for building binary abundance matrix.

#### binary_abundance_spoofed_minpath.py	
Builds a binary matrix from minpath output converted to CD-HIT .clstr format with parse_minpath_output_to_cdhit_format.py. 

#### dna_2_protein.py
Translates in-frame dna sequences to protein sequences for a given fasta file.

#### iterative_cdhit.sh
Iteratively runs CD-HIT from 100% amino acid identity down to 50%.

#### parse_iterative_cdhit.py
Processes iterative CD-HIT output to compute gene cluster sizes at different percent identities. 

#### remove_short_genes_gc.py	
Removes genes shorter than 100 nucleotides from gene catalog. 

#### prep_cluster_file_for_filtering.sh	
First step in removing short genes (below 100bp) from CD-HIT .clstr file.

#### remove_short_genes_from_cluster.py
Second step in removing short genes (below 100bp) from CD-HIT .clstr file.

#### sorensen.cpp
Sorenson-Dice distance calculation for all genes in a binary abundance matrix. Needs to be compiled with openmp first. (g++ -fopenmp sorensen.cpp)
