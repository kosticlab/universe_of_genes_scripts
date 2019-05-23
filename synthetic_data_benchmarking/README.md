## Directory Overview

You can find here the scripts we used to run our synthetic data analysis found in the supplemental materials. In short, we generated synthetic metagenomes from fully sequenced complete genomes, assembled our metagenomes with a variety of parameters/assemblers, identified singletons and non-singletons, identified false and true positive genes, and searched for associations between gene n-ton status and a variety of parameteres (i.e. assembler type, depth of sequencing, etc).

## Directory Contents

#### download_homd_data.py 	
Acquire complete genomes from HOMD.

#### art_parameters_example.sh 
Sample parameters for synthetic metagenome construction.

#### run_synthetic_data_modeling.py 
Assemble synthetic metagenomes.

#### cat_summary_data.py 	
Combine summary data for multiple runs of synthetic metagenome processing into a single file.

#### compute_gene_contig_coverage.sh 	
Compute read fold coverage for all predicted genes/assembled contigs.

#### run_synthetic_cdhit_analysis.py 	
Cluster predicted from metagenomes genes and identify singletons/non-singletons.


