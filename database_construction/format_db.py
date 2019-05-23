###to format tsv files into initial db structure


import os
import pandas as pd
import numpy as np
from random import choice

#def create_random_seq(geneLength):
#    DNA=""
#    for count in range(geneLength):
#        DNA+=choice("CGTA")
#    return DNA

def load_gc(filename):
    gcDict={}
    with open(filename) as f:
        seqs=[]
        for line in f:
            if '>' in line:
                if len(seqs)>0:
                    sequence=''.join(seqs)
                    gcDict[geneId[1:]]=sequence
                    seqs=[]
                geneId=line.rstrip()
            else:
                seqs.append(line.rstrip())
    sequence=''.join(seqs)
    gcDict[geneId]=sequence
    return gcDict

def parse_cluster_file(clusterFile):
    gene_congene_mapping={}
    cluster_size_mapping={}
    for line in clusterFile:
        subGenes=[]
        if '>' in line:
            if len(subGenes)>1:
                for x in subGenes:
                    congene_gene_mapping[congene]=x

                    cluster_size_mapping[x]=len(subGenes)
                subGenes=[]
        else:
            if '*' in line:
                congene=line.rstrip().split('>')[1].split('...')[0]
            subGenes.append(line.rstrip().split('>')[1].split('...')[0])
    for x in subGenes:
        congene_gene_mapping[congene]=x
        try:
            gene_congene_mapping.setdefault(congene, gene_congene_mapping[congene]).append(x)
        except:
            gene_congene_mapping[x]=congene
    cluster_size_mapping[x]=len(subGenes)
    return [gene_congene_mapping,congene_gene_mapping,congene_cluster_size_mapping]

def parse_cluster_file(clusterFile,bodySite):
    cluster_size_mapping={}
    congene_gene_mapping={}
    bodySiteMapping={}
    subGenes=[]
    with open(clusterFile) as f:
        for line in f:
            if 'Cluster' in line:
                if len(subGenes)>0:
                    for x in subGenes:
                        cluster_size_mapping[x]=len(subGenes)
                        bodySiteMapping[x]=bodySite
                        try:
                            congene_gene_mapping.setdefault(congene, gene_congene_mapping[congene]).append(x)
                        except:
                            congene_gene_mapping[x]=congene
                            subGenes=[]
            else:
                if '*' in line:
                    congene=line.rstrip().split('>')[1].split('...')[0]
                subGenes.append(line.rstrip().split('>')[1].split('...')[0])
        for x in subGenes:
            cluster_size_mapping[x]=len(subGenes)
            bodySiteMapping[x]=bodySite
            try:
                congene_gene_mapping.setdefault(congene, gene_congene_mapping[congene]).append(x)
            except:
                congene_gene_mapping[x]=congene
        return [cluster_size_mapping,congene_gene_mapping,bodySiteMapping]

def get_body_site_nton(gene):
    count=0
    try:
        count+=clusterSizesOral[gene]
    except:
        pass
    try:
        count+=clusterSizesGut[gene]
    except:
        pass
    try:
        count+=clusterSizesAll[gene]-1
    except:
        pass
    bodySites=[]
    try:
        bodySites.append(geneMappingOral[gene])
    except:
        pass
    try:
        bodySites.append(geneMappingGut[gene])
    except:
        pass
    bodySites='|'.join(list(set(bodySites)))
    return [count,bodySites]


#load fasta file, mapping geneIDs to sequence
seqDict=load_gc('gene_catalog_oral_gut_95_nucl')

#get cluster sizes and body site mappings
clusterSizesOral,congeneGeneOral,geneMappingOral=parse_cluster_file('oral_1473_gene_catalogue_filtered_sequences.clstr','oral')
clusterSizesGut,congeneGeneGut,geneMappingGut=parse_cluster_file('gene_catalog_95perc_filtered.clstr','gut')
clusterSizesAll,congeneGeneAll,geneMappingAll=parse_cluster_file('gene_catalog_oral_gut_95_nucl.clstr','oral_gut')

#load taxa data
gene_taxon={}
with open('all_taxa_map') as f:
    for line in f:
        line=line.rstrip().split('\t')
        gene_taxon[line[0]]=line[1]

initialDB=[]
#parse  daa from annotation software
with open('tsv_output') as f:
  for line in f:
    line=line.rstrip().split('\t')
    try:
        test=seqDict[line[0]]
        ntonCount,bodySites=get_body_site_nton(line[0])
    except:
        continue
    if len(line)==3:
      line=line[:2]+['','']+[line[-1]]
      line.append(gene_taxon[line[0]])
      line.append(ntonCount)
      line.append(bodySites)
      line.append(len(seqDict[line[0]]))
      line.append(seqDict[line[0]])
      initialDB.append(line)
    elif len(line)==4:
      if len(line[2].split('.'))==4:
          line=line[:2]+['']+line[-2:]
          line.append(gene_taxon[line[0]])
          line.append(ntonCount)
          line.append(bodySites)
          line.append(len(seqDict[line[0]]))
          line.append(seqDict[line[0]])
          initialDB.append(line)
      else:
          initialDB.append(line)     
          line=line[:3]+['']+[line[-1]]
          line.append(gene_taxon[line[0]])
          line.append(ntonCount)
          line.append(bodySites)
          line.append(len(seqDict[line[0]]))
          line.append(seqDict[line[0]])
          initialDB.append(line)
    elif len(line)==5:
      line.append(gene_taxon[line[0]])
      line.append(ntonCount)
      line.append(bodySites)
      line.append(len(seqDict[line[0]]))
      line.append(seqDict[line[0]])
      initialDB.append(line)

initialDB=pd.DataFrame(initialDB)

initialDB.columns=['GENE ID','CDS','ANNOTATION','ECID','FULL GENE NAME','NCBI TAXON ID','NUMBER_OF_GENES_IN_CLUSTER','BODY_SITES','GENE_LENGTH','SEQUENCE']
initialDB=initialDB.drop(initialDB.columns[1],axis=1)

ecids=list(initialDB.ECID)
for i,x in enumerate(ecids):
    if '.' not in x or x!='':
        ecids[i]=''

initialDB.loc[:,'ECID']=ecids

initialDB=initialDB.dropna()

#write to file
initialDB.to_csv('gene_database_oral_gut_full.csv',sep=',')
