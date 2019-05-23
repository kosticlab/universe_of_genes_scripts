###fully process genecat and prep all output files for figure gen
##20181015
#Tierney

import os
from remove_short_genes_gc import remove_short
from remove_short_genes_from_cluster import remove_short_cluster
from dna_2_protein import dna_2_protein
from find_cluster_distribution_50perc import *
from merge_oral_nton_taxa_data import parse_contig_database
from binary_abundance_matrix import build_abundance_matrix,load_cluster_mapping,load_prokka_dict
from binary_abundance_spoofed_minpath import *
from parse_minpath_output_to_cdhit_format import prep_minpath_data
#arguments: location of all orfs, list of all gff files from prokka output, location of raw sequence files, nr database made with taxonmap, ncbi taxon nodes and taxid files
#minpath folder in working directory and properly functioning
#pandas, numpy installed in addition to base python packages
#compiled sorenson dice code in working directory
#cdhit and diamond installed in working directory
#prokka id mapping

###TODO: #Extract 95% cluster size distribution
		 #Figure out how to point to gff files

#run cdhit; 1 day
os.system('./cdhit/cd-hit-est -i allGenes -o gene_catalog_95perc -c 0.95 -n 10 -T 0 -aS .9 -s .9 -M 0')

#filter out short genes fewer than 100bp long; 30 minutes
remove_short('gene_catalog_95perc','gene_catalog_95perc_filtered_sequences','gene_catalog_95perc_removed_sequences','nr.dmnd','./diamond')
os.system('chmod +x prep_cluster_file_for_filtering.sh; ./prep_cluster_file_for_filtering.sh gene_catalog_95perc.clstr gene_catalog_95perc_removed_sequences')
remove_short_cluster('gene_catalog_95perc')

#find taxonomies and parse out to file; 24 hours
os.system('./diamond blastx -d ./nr.dmnd -q gene_catalog_95perc_filtered_sequences -o gene_catalog_taxa_mapping --taxonmap ./prot.accession2taxid.gz --taxonnodes ./nodes.dmp --outfmt 102 -k 1')

#use following line to generate taxonomy mapping file from ncbi data
#os.system("grep 'scientific name' names.dmp  | cut -f1,2 -d'|' | sed 's/|//g' | sed 's/\t\t/\t/g'  > ncbi_taxa_ID_mapping.tsv")

#translate gene catalog; 1 hour
dna_2_protein('gene_catalog_95perc_filtered_sequences','gene_catalog_filtered_sequences_translated_protein')

#run iterative cdhit; 1 day
os.system('chmod +x iterative_cdhit.sh; ./iterative_cdhit.sh ./cdhit/cd-hit')

#map genes between percent identities; 1 hour
run_iterative_processing('iterative_cdhit_output/')

#build binary abundance matrices for 95% and 50%; 1 day, 12 hours, respectively
build_abundance_matrix('iterative_cdhit_output/gene_catalog_95perc_filtered','prokka_id_mapping','95perc_orf_mapping.csv','gene_catalog_binary_abundance_95perc.csv')
build_abundance_matrix('iterative_cdhit_output/gene_catalog_filtered_sequences_translated_protein_50_perc','prokka_id_mapping','50perc_orf_mapping.csv','gene_catalog_binary_abundance_50perc.csv')
#build sorensen-dice indices for 95%, 50%, and minpath; will take in total about 2 days
os.system('./sorenson gene_catalog_binary_abundance_50perc.csv;  mv sorenson.csv gene_catalog_sd_50perc.csv')
os.system('./sorenson gene_catalog_binary_abundance_95perc.csv;  mv sorenson.csv gene_catalog_sd_95perc.csv')

######################################building and processing contig pipeline########################
#build contig database - 1 day
#os.system('chmod +x parallel_contig_parse.sh; ./parallel_contig_parse.sh')
#os.system('cat gene_catalog_contig_data_50perc_*_*.csv > gene_catalog_contig_data_50perc.csv')

#sort contigs into different groups and generate files for minpath and taxa/ecid R analysis; 10 minutes
#parse_contig_database()

#run minpath and build associated binary abundance matrix + sorensen matrix (1 minute)
#os.system('python MinPath/MinPath1.4.py -any gene_catalog_ecids_for_minpath -map ec2path -report gene_catalog_50perc.ec.report -details gene_catalog_50perc.ec.details')
#prep_minpath_data('gene_catalog_50perc.ec.report','gene_catalog_50perc.ec.details','gene_catalog_raw_orf_ecid_mapping.tsv')
#build_abundance_matrix_minpath('minpath_95','prokka_id_mapping','gene_catalog_binary_abundance_minpath_95.csv')
#os.system('./sorenson gene_catalog_binary_abundance_minpath_95.csv; mv sorenson.csv gene_catalog_sd_minpath_95.csv')
#build_abundance_matrix_minpath('minpath_50','prokka_id_mapping','gene_catalog_binary_abundance_minpath_50.csv')
#os.system('./sorenson gene_catalog_binary_abundance_minpath_50.csv; mv sorenson.csv gene_catalog_sd_minpath_50.csv')
