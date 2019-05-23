#remove lines that are too short from abundance matrix

import os
import sys

#find and load gene ID's to remove
os.system("grep '>' oral_1476_gene_catalogue_removed_sequences | cut -d' ' -f1 > temp2")

short=[]
f=open('temp2')
for line in f:
	short.append(line.rstrip()[1:])

short=set(short)

f.close()

#iterate through abundance matrix
#if a gene ID is in the list of those to remove, don't write it to the file
f=open(sys.argv[1])
w=open(sys.argv[1]+'_filtered','w')
for line in f:
	if line.split(',')[0] in short:
		continue
	else:
		w.write(line)

f.close()
w.close()

os.system('rm temp2')
