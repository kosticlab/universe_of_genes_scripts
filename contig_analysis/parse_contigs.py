###build contig database
##20180911
#Tierney

from collections import Counter
import pandas as pd
import sys
import os
import multiprocessing as mp

###create tsv of information on all contigs in dataset
#parse prokka gff files, extract taxonomic+functional+singleton/nonsingleton information mapping gene ids to contigs
def parse_gff(foldername):
	files=os.listdir(foldername)
	filename=[x for x in files if 'gff' in x][0]
	fxnDict={}
	contigData=[]
	genes=[]
	taxa=[]
	singletons=[]
	genes2=[]
	done=[]
	ecids=[]
	nonsingletons=[]
	oldcontig=''
	ntonStatuses=[]
	with open(foldername+'/'+filename) as f:
	  print 'Parsing '+foldername+'/'+filename+'...'
	  for i,line in enumerate(f):
		#confirm this is not a header or tail line in the gff file
		line=line.rstrip().split('\t')
		#make sure there is data in the line
		if len(line)<8:
		  continue
		if len(line)>1:
		  if 'minced' in str(line):
			  continue
		  #get contig id
		  contig=line[0]
		  #create unique contig ID from sample ID and contig identifier
		  contigName=[line[0]+'_'+line[8].split('ID=')[1].split('_')[0]]
		  #if the contig id is different from before, reinitialize lists and store data
		  if contig!=oldcontig and i!=0:
			for i,g in enumerate(genes):
				assert(len(ecids)==len(genes))
				if ecids[i]!='':
					try:
						fxnDict.setdefault(g, fxnDict[g]).append(ecids[i])
					except:
						fxnDict[g]=ecids[i]
				try:
					taxa.append(taxaDict[g])
				except:
					taxa.append('')
					continue
			if len(genes)>0:
				contigData.append([contigName,genes,genes2,ntonStatuses,singletons,nonsingletons,ecids,taxa])
			#build functional annotation dictionary
			genes=[]
			genes2=[]
			singletons=[]
			nonsingletons=[]
			taxa=[]
			ecids=[]
			ntonStatuses=[]
		  #get the gene described on this line of the gff file
		  tempGene=line[8].split('ID=')[1].split(';')[0].split('_gene')[0]
		  #if we have seen these gene before but NOT been able to assign an EC ID, attempt to do so and continue to next iteration
		  if tempGene in done:
			if 'Prodigal' in str(line) and 'eC_number' in str(line):
			  if ecids[-1]=='':
				ecids[-1]=line[8].split('eC_number=')[1].split(';')[0]
			continue
		  #find consensus gene
		  try:
			genes.append(geneMapping[tempGene])
			genes2.append(tempGene)
		  except:
			continue
		  #find singleton status
		  nton=ntonStatus[geneMapping[tempGene]]
		  ntonStatuses.append(nton)
		  if nton=='singleton':
			singletons.append(geneMapping[tempGene])
		  if nton=='non-singleton':
			nonsingletons.append(geneMapping[tempGene])
		  #look for EC ID
		  if 'Prodigal' in str(line) and 'eC_number' in str(line):
			ecids.append(line[8].split('eC_number=')[1].split(';')[0])
		  if 'Prodigal' in str(line) and 'eC_number' not in str(line):
			ecids.append('')
		  if 'prokka' in str(line):
			ecids.append('')
		  if 'Aragorn' in str(line):
			ecids.append('')
		  oldcontig=line[0]
		  done.append(tempGene)
	  if len(genes)>0:
		contigData.append([contigName,genes,genes2,ntonStatuses,singletons,nonsingletons,ecids,taxa])
	  return contigData

geneMapFile=sys.argv[1] #comma separated file containing two columns: left is each raw gene ID, right is appropriate consensus gene
ntonMapFile=sys.argv[2] #comma separated file containing two columns: left is each raw gene ID, right is singleton status (singleton or non-singleton)
fileToProcess=sys.argv[3] #file containing location of all prokka output folders (each of which, in turn, should have a gff file)
taxaMapping=sys.argv[4] #comma separated file containing two columns: left is each raw gene ID, right is ncbi taxonomic annotation
bottom=int(sys.argv[5]) #for parallel processing: gff file number to start at
top=int(sys.argv[6]) #for parallel processing: gff file number to end at

#for bottom and top arguments – to process all at once, set bottom at 0 and top at total number of gff files

#identify paths to prokka output
folders=[]
with open(fileToProcess) as f:
	for line in f:
		folders.append(line.rstrip())

#load taxonomic data
taxaDict={}
with open(taxadata) as f:
	for line in f:
		taxaDict[line.rstrip().split(',')[0]]=line.rstrip().split(',')[1]

#load list of files to work with

#load nton status data
ntonStatus={}
with open(ntonMapFile) as f:
  for line in f:
    ntonStatus[line.split(',')[0]]=line.rstrip().split(',')[1]

#load the geneMapping
geneMapping={}
with open(geneMapFile) as f:
  for line in f:
    geneMapping[line.split(',')[0]]=line.rstrip().split(',')[1]

#iterate through files, extract necessary information
#parallel version for speed

pool = mp.Pool(16)
allData = pool.map(parse_gff, folders[bottom:top])

#unlist
allGffData = [item for sublist in allData for item in sublist]

#store full gff output and compress it to make extra space
with open('gene_catalog_contig_data_50perc_%s_%s.csv'%(bottom,top),'w') as w:
	for line in allGffData:
		line=['|'.join(x) for x in line]
w.write(','.join(line)+'\n')

#will need to cat output from parallel processing together
