###run cdhit and create new analysis table
##20190221
#Tierney

import os
import pandas as pd

#to be run from within a particular coverage folder

def get_singletons(filename):
    singletonMap=[]
    samps=[]
    with open(filename) as f:
        for line in f:
        	if 'Cluster' in line:
        		if len(samps)>0:
        			singletonMap.append([congene,len(samps)])
        			samps=[]
        	else:
        		samps.append(line.rstrip().split('>')[1].split('_')[0])
        		if '*' in line:
        			congene=line.rstrip().split('>')[1].split('...')[0]
    singletonMap.append([congene,len(samps)])
    return singletonMap

def run_cdhit_and_merge_output(tag,gene_by_gene,prokka_folder):
    os.system('cat */%s/*faa > all_genes_for_cdhit_%s'%(prokka_folder,tag))
    os.system('cat */%s > gene_by_gene_%s.csv'%(gene_by_gene,tag))
    os.system('../cdhit/cd-hit -n 3 -i all_genes_for_cdhit_%s -T 0 -M 0 -s .9 -aS .9 -c .5 -o cdhit_output_50perc_%s'%(tag,tag))
    map50=pd.DataFrame(get_singletons('cdhit_output_50perc_%s.clstr'%tag))
    map50.columns=['name','50_perc_status']
    mh_df=pd.read_csv('gene_by_gene_%s.csv'%tag)
    merged=pd.merge(map50,mh_df,left_on='name',right_on='Gene ID',how='left')
    merged.to_csv('gene_by_gene_with_singleton_data_%s.csv'%tag)

run_cdhit_and_merge_output('megahit_default','prokka_output_megahit_contigs_default_gene_by_gene.csv','prokka_output_megahit_contigs_default')
run_cdhit_and_merge_output('megahit_sensitive','prokka_output_megahit_contigs_sensitive_gene_by_gene.csv','prokka_output_megahit_contigs_sensitive')
run_cdhit_and_merge_output('megahit_large','prokka_output_megahit_contigs_large_gene_by_gene.csv','prokka_output_megahit_contigs_large')
run_cdhit_and_merge_output('metaspades','prokka_output_metaspades_contigs_gene_by_gene.csv','prokka_output_metaspades_contigs')
