#!/bin/python

###metagenomic phylogeny analysis
##20180522
#Tierney

#after running cdhit at different percent identities (can use iterative cdhit script), identify number of singletons in each file

import os
import sys
import pandas as pd
from collections import Counter
import subprocess

#identify singleton genes in catalog cluster file
def find_singletons(f):
	ids = []
	output=[]
	cluster_catalogue = open(cdhitoutdir+'/'+f, "r")
	for line in cluster_catalogue:
		if line[0] != '>':
			ids.append(line.split('>')[1].split('...')[0].split("_")[0])
		else:
			if len(ids)>0:
				output.append(ids)
				ids=[]
			else:
				continue
	output.append(ids)
	output=[len(set(x)) for x in output]
	return output

cdhitoutdir='cdhitout'

#for each cdhit run, find and count singleton vs non-singleton genes
#only look at up to 9-tons, then group everything together
files=os.listdir(cdhitoutdir)
files.sort()
files=list(reversed(files))
counts=[]
for i,f in enumerate(files):
	print f
	out=find_singletons(f)
	out=Counter(out)
	out=pd.DataFrame.from_dict(out,orient='index')
	out.columns=['number_%s'%f]
	if i==0:
		out.loc[:,'CLUSTER_SIZE']=[int(x) for x in out.index]
		out2=out.iloc[0:9,:]
		bottom=out.iloc[9:,:].sum()
		bottom.iloc[1]='10+'
		out2=out2.append(bottom,ignore_index=True)
		out2=out2.reindex_axis(sorted(out2.columns), axis=1)
		final=out2
		continue
	else:
		out2=out.iloc[0:9,:]
		bottom=out.iloc[9:,:].sum()
		out2=out2.append(bottom,ignore_index=True)
		final=pd.concat([final,out2.iloc[:,0]],axis=1)
		continue

final.to_csv('oral_iterative_cdhit_out.csv')
