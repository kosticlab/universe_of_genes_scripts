#open and parse sequences + names + gene lengths from gene catalog
def remove_short(inputName,outputNameFilt,outputNameRem,diamondIndex,diamond):
	import sys
	import numpy as np
	import os
	names=[]
	seqs=[]
	sequences=[]
	f=open(inputName)
	for line in f:
		if '>' in line:
			names.append(line.rstrip())
			if len(seqs)>0:
				sequences.append(''.join(seqs))
				seqs=[]
		else:
			seqs.append(line.rstrip())
	sequences.append(''.join(seqs))
	lengths=np.array([len(x) for x in sequences])
	print 'GENE MEAN LENGTH PRIOR TO REMOVAL OF SHORT SEQUENCES:'+str(np.mean(lengths))
	#find long sequences
	saved=[]
	saved_names=[]
	for i,x in enumerate(sequences):
		if len(x)>=100:
			saved_names.append(names[i])
			saved.append(x)
	#find short sequences
	removed=[]
	removed_names=[]
	for i,x in enumerate(sequences):
		if len(x)<100:
			removed_names.append(names[i])
			removed.append(x)
	removed_names=[x.split(' ')[0] for x in removed_names]
	#temporarily write short sequences
	w=open('removed_sequences_temp','w')
	for i,line in enumerate(removed):
		w.write(removed_names[i]+'\n')
		w.write(line+'\n')
	w.close()
	#align short sequences back to refseq
	os.system("%s blastx --query 'removed_sequences_temp' --id '95' -d '%s' -k '1' -o 'removed_sequences_temp_aligned'"%(diamond,diamondIndex))
	#find aligned sequences
	validNames=[]
	f=open('removed_sequences_temp_aligned')
	for line in f:
		validNames.append('>'+line.split('\t')[0])
	#remove temp tiles
	os.system('rm removed_sequences_temp')
	os.system('rm removed_sequences_temp_aligned')
	validNames=list(set(validNames))
	setvalid=set(validNames)
	validNamesIndices=[removed_names.index(x) for x in validNames]
	#identify sequences that aligned
	validSeq=[]
	for x in validNamesIndices:
		validSeq.append(removed[x])
	setValidNamesIndices=set(validNamesIndices)
	#append aligned short sequences to sequences to save list
	for i,name in enumerate(validNames):
		saved_names.append(name)
		saved.append(validSeq[i])
	saved_names=[x.split(' ')[0] for x in saved_names]
	#filter removed genes list by aligned short sequences
	removed_names=[x for i,x in enumerate(removed_names) if i not in setValidNamesIndices]
	removed=[x for i,x in enumerate(removed) if i not in setValidNamesIndices]
	print 'GENE MEAN LENGTH AFTER REMOVAL OF SHORT SEQUENCES:'+str(np.mean([len(x) for x in saved]))
	#write saved and removed genes to file
	w=open(outputNameRem,'w')
	for i,line in enumerate(removed):
		w.write(removed_names[i]+'\n')
		w.write(line+'\n')
	w.close()
	w=open(outputNameFilt,'w')
	for i,line2 in enumerate(saved):
		w.write(saved_names[i]+'\n')
		w.write(line2+'\n')
	w.close()
