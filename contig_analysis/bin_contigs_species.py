###bin contigs from the same species
##20180915
#Tierney

#find contigs mapping to each species given a custom contig database
#you can select number of genes required to present for each contig to be considered

import pandas as pd
import sys
from collections import Counter

def bin_species(contigData,taxaMap,mingenecontent,mixed):
	classDict1={}
	classDictCount1={}
	classDict2={}
	classDictCount2={}
	classDict1_diff={}
	classDictCount1_diff={}
	classDict2_diff={}
	classDictCount2_diff={}
	unused=[]
	unused1=[]
	unused2=[]
	with open(contigData) as f:
		for i,line in enumerate(f):
			species=[]
			speciesNonSing=[]
			speciesSing=[]
			diff=False
			print i
			line=line.rstrip().split(',')
			line=[x.split('|') for x in line]
			if len(line[1])>mingenecontent:
				if mixed==True:
					singledexes=[]
					for i,val in enumerate(line[3]):
						if val=='Singleton' or val == 'singleton':
							singledexes.append(i)
					species=list(line[-1])
					speciesSing=[x for i,x in enumerate(species) if i in singledexes]
					speciesNonSing=[x for i,x in enumerate(species) if i not in singledexes]
					speciesNonSing=[x for x in speciesNonSing if x!='' or x!='0' or x!='2']
					speciesSing=[x for x in speciesSing if x!='' or x!='0' or x!='2']
					species=[x for x in species if x!='' or x!='0' or x!='2']
				else:
					species=list(set(line[-1]))
					species=[x for x in species if x!='']
				speciesCount=pd.DataFrame.from_dict(Counter(species),orient='index')
				if len(speciesCount.index)==0:
					continue
				speciesCount=speciesCount/float(speciesCount.sum()[0])
				maxVal=speciesCount.iloc[:,0].max()
				if maxVal>=.75:
					maxIdx=speciesCount.iloc[:,0].idxmax()
					species=[maxIdx]
					try:
						speciesName=taxaMap[species[0]]
					except:
						speciesName=species[0]
					if mixed==True and speciesName!='Bacteria':
						for i,x in enumerate(speciesSing):
							try:
								speciesSing[i]=taxaMap[x].split(' ')[0]
							except:
								continue
						if speciesName.split(' ')[0] not in speciesSing and len(speciesSing)>0:
							diff=True
					try:
						#classDict[line[0][0]]=speciesName
						if len(line[1])>2 and diff==False:
							classDict2[line[0][0]]=speciesName
							classDictCount2.setdefault(speciesName, classDictCount[speciesName]).append(line[0][0])
						if len(line[1])>2 and diff==True:
							classDict2_diff[line[0][0]]=speciesName
							classDictCount2_diff.setdefault(speciesName, classDictCount2_diff[speciesName]).append(line[0][0])
					except:
						#classDictCount[speciesName]=[line[0][0]]
						if len(line[1])>2 and diff==False:
							classDictCount2[speciesName]=[line[0][0]]
						if len(line[1])>2 and diff==True:
							classDictCount2_diff[speciesName]=[line[0][0]]
					try:
						if len(line[1])>1 and diff==False:
							classDict1[line[0][0]]=speciesName
							classDictCount1.setdefault(speciesName, classDictCount1[speciesName]).append(line[0][0])
						if len(line[1])>1 and diff==True:
							classDict1_diff[line[0][0]]=speciesName
							classDictCount1_diff.setdefault(speciesName, classDictCount1_diff[speciesName]).append(line[0][0])
					except:
						if len(line[1])>1 and diff==False:
							classDictCount1[speciesName]=[line[0][0]]
						if len(line[1])>1 and diff==True:
							classDictCount1_diff[speciesName]=[line[0][0]]
				else:
				#	unused.append(line)
					if len(line[1])>1:
						unused1.append(line)
					if len(line[1])>2:
						unused2.append(line)
		#for key in classDictCount.keys():
		for key in classDictCount1.keys():
			classDictCount1[key]=len(classDictCount1[key])
		for key in classDictCount2.keys():
			classDictCount2[key]=len(classDictCount2[key])
		for key in classDictCount1_diff.keys():
			classDictCount1_diff[key]=len(classDictCount1_diff[key])
		for key in classDictCount2_diff.keys():
			classDictCount2_diff[key]=len(classDictCount2_diff[key])
	return classDict1,classDictCount1,classDict1_diff,classDictCount1_diff,unused1,classDict2,classDictCount2,classDict2_diff,classDictCount2_diff,unused2

def run_parser(contigData,taxaData,number,nton):
	mixture=False
	if 'mixed' in nton:
		mixture=True
	binnedContigData1,binnedContigDataCount1,binnedContigData1_diff,binnedContigDataCount1_diff,notUsed1,binnedContigData2,binnedContigDataCount2,binnedContigData2_diff,binnedContigDataCount2_diff,notUsed2=bin_species(contigData,taxaData,number,mixture)
	#convert to dataframes and send to tsvs
	binnedContigData1=pd.DataFrame.from_dict(binnedContigData1,orient='index')
	binnedContigData1.to_csv('gene_catalog_binned_contig_data_%s_1.csv'%(nton),sep='\t')
	binnedContigDataCount1=pd.DataFrame.from_dict(binnedContigDataCount1,orient='index')
	binnedContigDataCount1.to_csv('gene_catalog_binned_contig_data_%s_1_count.csv'%(nton),sep='\t')
	binnedContigData2=pd.DataFrame.from_dict(binnedContigData2,orient='index')
	binnedContigData2.to_csv('gene_catalog_binned_contig_data_%s_2.csv'%(nton),sep='\t')
	binnedContigDataCount2=pd.DataFrame.from_dict(binnedContigDataCount2,orient='index')
	binnedContigDataCount2.to_csv('gene_catalog_binned_contig_data_%s_2_count.csv'%(nton),sep='\t')
	binnedContigData1_diff=pd.DataFrame.from_dict(binnedContigData1_diff,orient='index')
	binnedContigData1_diff.to_csv('gene_catalog_binned_contig_data_%s_1_diff.csv'%(nton),sep='\t')
	binnedContigDataCount1_diff=pd.DataFrame.from_dict(binnedContigDataCount1_diff,orient='index')
	binnedContigDataCount1_diff.to_csv('gene_catalog_binned_contig_data_%s_1_count_diff.csv'%(nton),sep='\t')
	binnedContigData2_diff=pd.DataFrame.from_dict(binnedContigData2_diff,orient='index')
	binnedContigData2_diff.to_csv('gene_catalog_binned_contig_data_%s_2_diff.csv'%(nton),sep='\t')
	binnedContigDataCount2_diff=pd.DataFrame.from_dict(binnedContigDataCount2_diff,orient='index')
	binnedContigDataCount2_diff.to_csv('gene_catalog_binned_contig_data_%s_2_count_diff.csv'%(nton),sep='\t')
	if nton=='mixed_singleton':
		with open('gene_catalog_binned_contig_data_%s_%s_unused_contigs.csv'%(nton,1),'w') as w:
			for line in notUsed1:
				line=['|'.join(x) for x in line]
				line=','.join(line)+'\n'
				w.write(line)
		with open('gene_catalog_binned_contig_data_%s_%s_unused_contigs.csv'%(nton,2),'w') as w:
			for line in notUsed2:
				line=['|'.join(x) for x in line]
				line=','.join(line)+'\n'
				w.write(line)


contigsNonSingle='gene_catalog_contig_data_nonsingletons_50perc.csv'
contigsMixed='gene_catalog_contig_data_mixed_50perc.csv'
contigsSingle='gene_catalog_contig_data_singletons_50perc.csv'

#contigsNonSingle='oral_contig_data_nonsingletons_50perc.csv'
#contigsMixed='oral_contig_data_mixed_50perc.csv'
#contigsSingle='oral_contig_data_singletons_50perc.csv'


#load taxonomy data
taxaNameMap={}
with open('ncbi_taxa_ID_mapping.tsv') as f:
  for line in f:
	taxaNameMap[line.split('\t')[0]]=line.rstrip().split('\t')[1]

#run_parser(contigsMixedSingle,taxaNameMap,0,'mixed_singleton')

run_parser(contigsSingle,taxaNameMap,0,'singleton')
run_parser(contigsMixed,taxaNameMap,0,'mixed')
run_parser(contigsNonSingle,taxaNameMap,0,'non_singleton')
