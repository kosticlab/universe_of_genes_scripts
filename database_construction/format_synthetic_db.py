###to format tsv files into initial db structure


import os
import pandas as pd
import numpy as np
from random import choice

def create_random_seq(geneLength):
    DNA=""
    for count in range(geneLength):
        DNA+=choice("CGTA")
    return DNA

#os.system("cat */*tsv > all_prokka_tsv_output ; grep 'CDS' all_prokka_tsv_output > foo ; mv foo all_prokka_tsv_output")
initialDB=[]
#parse  daa from annotation software
with open('all_prokka_tsv_output') as f:
  for line in f:
    line=line.rstrip().split('\t')
    if len(line)==3:
      line=line[:2]+['','']+[line[-1]]
      initialDB.append(line)
    elif len(line)==4:
      line=line[:3]+['']+[line[-1]]
      initialDB.append(line)
    elif len(line)==5:
      initialDB.append(line)

#convert list of lists to dataframe for ease of use
initialDB=pd.DataFrame(initialDB)
initialDB.columns=['GENE ID','CDS','ANNOTATION','ECID','GENE NAME']
#drop irrelevant column
initialDB.drop('CDS',axis=1,inplace=True)
#get length of dataframe for generation of additional columns
rowLength=len(initialDB.index)

create synthetic body site data
bs=['oral']*(int(rowLength)/2) + ['gut']*int(rowLength)
initialDB.loc[:,'BODY SITE']=bs[:rowLength]

create synthetic nton status data
ntonstat=np.random.randint(1,2000,size=rowLength)
initialDB.loc[:,'NUMBER OF GENES IN CLUSTER']=ntonstat

create length data
lengths=np.random.randint(100,2000,size=rowLength)
initialDB.loc[:,'GENE LENGTH']=lengths

create consensus sequence data
seqs=[]
for x in lengths:
    seqs.append(create_random_seq(x))

initialDB.loc[:,'CONSENSUS SEQUENCE']=seqs

create study id data
sid=['study 1']*(int(rowLength)/4) + ['study 2']*(int(rowLength)/4)+ ['study 3']*(int(rowLength)/4)+ ['study 4']*(int(rowLength)/3)
initialDB.loc[:,'STUDY ID']=sid[:rowLength]

#write to file
initialDB.to_csv('synthetic_database_table_for_el.tsv',sep='\t')
