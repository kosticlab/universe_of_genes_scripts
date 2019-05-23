#!/usr/bin/python
def map_prokka_names_to_file_minpath():
	prokkaOutput=os.listdir('/mnt/oral/contigparse/aminoithitcontigparse/oralmicrobiome')
	prokkaIDs=[]
	d={}
	for each0 in prokkaOutput:
		files=os.listdir('/mnt/oral/contigparse/aminoithitcontigparse/oralmicrobiome/'+each0)
		files=[x for x in files if '.ffn' in x]
		for each in files:
			with open('/mnt/oral/contigparse/aminoithitcontigparse/oralmicrobiome/'+each0+'/'+each) as f:
				first_line = f.readline()
				d[first_line.split(' ')[0][1:].split('_')[0]]=each0.split('_prokka_out')[0]
				print d[first_line.split(' ')[0][1:].split('_')[0]]
	return d

def load_prokka_dict_minpath(nameMap):
	lines={}
	with open(nameMap) as f:
		for line in f:
			lines[line.rstrip().split('\t')[0]]=line.rstrip().split('\t')[1]
	return lines

def build_abundance_matrix_minpath(sequence,nameMap,outName): #opens and parses cluster file, finds genes in each cluster and maps them to samples, iteratively writing each line of binary matrix
	import os
	from collections import Counter
	import pandas as pd
	import subprocess
	print 'Loading prokka name mapping...'
	d = load_prokka_dict_minpath(nameMap)
	allvals=[]
	clusters=[]
	f=open('%s.clstr'%sequence)
	print 'Parsing cluster file...'
	count=-1
	for cluster in f:
		if 'Cluster' in cluster:
			if len(clusters)>0:
				count+=1
				#####modification for functionality at 50% ID level
				clusters=list(set([d[x.split('_')[0]] for x in clusters if x!='']))
				freq=Counter(clusters)
				df=pd.DataFrame.from_dict(freq, orient='index')
				df[df>1]=1
				df.columns=[name]
				df=df.reindex(d.values()).transpose()
				df=df.fillna(0)
				clusters=[]
				if count==0:
					with open(outName,'a') as f:
						df.to_csv(f,header=True)
				else:
					with open(outName,'a') as f:
						df.to_csv(f,header=False)
			name=cluster.rstrip().split(' ')[-1]
		else:
			if '*' in cluster:
				congene=cluster.split(',')[1].split('...')[0][2:]
			cluster=cluster.split(',')[1].split('...')[0][2:]
			####modification for funcionality at 50% ID level
			#clusters.append(d[cluster])
			clusters.append(cluster)
	if len(clusters)>0:
		count+=1
		#####modification for functionality at 50% ID level
		clusters=list(set([d[x.split('_')[0]] for x in clusters if x!='']))
		freq=Counter(clusters)
		df=pd.DataFrame.from_dict(freq, orient='index')
		df[df>1]=1
		df.columns=[name]
		df=df.reindex(d.values()).transpose()
		df=df.fillna(0)
		if count==0:
			with open(outName,'a') as f:
				df.to_csv(f,header=True)
		else:
			with open(outName,'a') as f:
				df.to_csv(f,header=False)


#build_abundance_matrix_minpath('minpath','prokka_id_mapping','gene_catalog_binary_abundance_minpath.csv')
