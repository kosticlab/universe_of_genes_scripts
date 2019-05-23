###parse minpath results to cdhit cluster format
##20181015
#Tierney

#load path to kegg mappings
def prep_minpath_data(report,details,map):
	import pandas as pd
	path_kegg={}
	with open(report) as f:
		for line in f:
			if 'minpath 0' in line:
				continue
			line=line.rstrip().split(' ')
			key=line[1]
			val=line[-1]
			try:
				path_kegg.setdefault(key, path_kegg[key]).append(val)
			except:
				path_kegg[key]=[val]
	#load kegg to ecid mappings
	kegg_ecid={}
	with open(details) as f:
		for line in f:
			if 'fam-found' in line:
				key=line.split('# ')[1].rstrip()
				continue
			else:
				val=line.split('# ')[1].rstrip()
			try:
				kegg_ecid.setdefault(key, kegg_ecid[key]).append(val)
			except:
				kegg_ecid[key]=[val]
	#load minpath to gene mappings
	ecid_gene={}
	with open(map) as f:
		for line in f:
			line=line.rstrip().split('\t')
			if len(line)==1:
				key='0'
				val=line[0]
			else:
				key=line[1]
				val=line[0]
			try:
				ecid_gene.setdefault(key, ecid_gene[key]).append(val)
			except:
				ecid_gene[key]=[val]
	#map genes to paths
	gene_path={}
	genePathList=[]
	bad=[]
	###WEIRD; THERE ARE SOME ECIDS (9 OUT OF ABOUT 1000) IN THE MINPATH OUTPUT THAT ARE NOT ACTUALLY IN MY DATASET
	for minpath in path_kegg.keys():
		keggs=path_kegg[minpath]
		for kegg in keggs:
			ecids=kegg_ecid[kegg]
			for ecid in ecids:
				try:
					genes=ecid_gene[ecid]
				except:
					bad.append(ecid)
					continue
				for gene in genes:
					try:
						gene_path.setdefault(minpath, gene_path[minpath]).append(gene)
					except:
						gene_path[minpath]=[gene]
					genePathList.append([gene,ecid,kegg,minpath])

	with open('gene_ecid_minpath_mapping.tsv','w') as w:
		for line in genePathList:
			w.write('\t'.join(line)+'\n')

	#load singleton mapping data
	geneNtonMap={}
	with open('50perc_orf_nton_status.csv') as f:
	  for line in f:
	    geneNtonMap[line.split(',')[0]]=line.rstrip().split(',')[1]

	map50={}
	with open('50perc_orf_mapping.csv') as f:
		for line in f:
			line=line.rstrip().split(',')
			map50[line[0]]=line[1]

	forCountTable=[]
	for x in genePathList:
		forCountTable.append([x[-1],geneNtonMap[x[0]]])

	#group together and collapse
	forCountTable=pd.DataFrame(forCountTable)
	forCountTable=forCountTable.pivot_table(index=forCountTable.columns[1],columns=forCountTable.columns[0],aggfunc=lambda x: len(x)).transpose()
	forCountTable.fillna(0,inplace=True)
	forCountTable.columns.name = None
	forCountTable=forCountTable.sort_values(by=forCountTable.columns[1],ascending=False)

	#write to file
	forCountTable.to_csv('gene_catalog_minpath_nton_counts_50perc.tsv',sep='\t')

	#spoof cluster file
	with open('minpath_50.clstr','w') as w:
		for value in gene_path.keys():
			w.write('>Cluster '+value+'\n')
			genes=gene_path[value]
			for i,gene in enumerate(genes):
				if i==0:
					w.write('0      aa, >'+map50[gene]+'... *'+'\n')
				else:
					w.write('0      aa, >'+map50[gene]+'...'+'\n')

	########repeat for 95% identity
	geneNtonMap={}
	with open('95perc_orf_nton_status.csv') as f:
	  for line in f:
	    geneNtonMap[line.split(',')[0]]=line.rstrip().split(',')[1]

	map95={}
	with open('95perc_orf_mapping.csv') as f:
		for line in f:
			line=line.rstrip().split(',')
			map95[line[0]]=line[1]

	forCountTable=[]
	for x in genePathList:
		forCountTable.append([x[-1],geneNtonMap[x[0]]])

	#group together and collapse
	forCountTable=pd.DataFrame(forCountTable)
	forCountTable=forCountTable.pivot_table(index=forCountTable.columns[1],columns=forCountTable.columns[0],aggfunc=lambda x: len(x)).transpose()
	forCountTable.fillna(0,inplace=True)
	forCountTable.columns.name = None
	forCountTable=forCountTable.sort_values(by=forCountTable.columns[1],ascending=False)

	#write to file
	forCountTable.to_csv('gene_catalog_minpath_nton_counts_95perc.tsv',sep='\t')

	#spoof cluster file
	with open('minpath_95.clstr','w') as w:
		for value in gene_path.keys():
			w.write('>Cluster '+value+'\n')
			genes=gene_path[value]
			for i,gene in enumerate(genes):
				if i==0:
					w.write('0      aa, >'+map95[gene]+'... *'+'\n')
				else:
					w.write('0      aa, >'+map95[gene]+'...'+'\n')

prep_minpath_data('gene_catalog_50perc.ec.report','gene_catalog_50perc.ec.details','gene_catalog_raw_orf_ecid_mapping.tsv')
