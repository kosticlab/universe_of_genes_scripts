###filter contigdb into contigs of different types
##20181015
#Tierney

def process(inputFile):
	allGffData=[]
	print 'Loading data...'
	with open(inputFile) as f:
		for line in f:
			line=line.rstrip().split(',')
			allGffData.append(line)
	print 'Sorting contigs...'
	singletigs,nonsingletigs,mixed=find_nton_contig_type(allGffData)

	print 'Writing nonsingleton data to file...'
	with open('gene_catalog_contig_data_nonsingletons_50perc.csv','w') as w:
		for line in nonsingletigs:
			w.write(','.join(line)+'\n')

	print 'Writing singleton data to file...'
	with open('gene_catalog_contig_data_mixed_50perc.csv','w') as w:
		for line in mixed:
			w.write(','.join(line)+'\n')

	print 'Writing mixed data to file...'
	with open('gene_catalog_contig_data_singletons_50perc.csv','w') as w:
		for line in singletigs:
			w.write(','.join(line)+'\n')

	os.system('cat gene_catalog_contig_data_singletons_50perc.csv gene_catalog_contig_data_mixed_50perc.csv > gene_catalog_contig_data_singletons_mixed_50perc.csv')

def find_nton_contig_type(contigData):
	s=[]
	ns=[]
	m=[]
	for l in contigData:
		if list(set(l[3].split('|')))[0]=='Singleton' and len(list(set(l[3].split('|'))))==1:
			s.append(l)
		if list(set(l[3].split('|')))[0]=='Non-singleton' and len(list(set(l[3].split('|'))))==1:
			ns.append(l)
		if len(list(set(l[3].split('|'))))==2:
			m.append(l)
	return [s,ns,m]

process('gene_catalog_contig_data_50perc.csv')












