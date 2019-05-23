
def remove_short_cluster(inputFile):
	from time import time
	import sys

	print "==========================================================================================="
	print "reading list of genes to remove..."
	start = time()
	with open("unaligned_sort.txt") as f:
	    to_remove = f.readlines()
	end = time()
	print "finished in {} seconds".format(end-start)
	print "==========================================================================================="

	print "removing whitespaces and newlines..."
	start = time()
	to_remove = [x.strip() for x in to_remove]
	end = time()
	print "finished in {} seconds".format(end-start)
	print "==========================================================================================="

	i = 0
	skip = 1
	done = 0
	length = len(to_remove)
	print length
	output = open(inputFile+"_filtered.clstr", "w")
	print "iterating through cluster file and removing..."
	start = time()
	with open(inputFile+".clstr") as f:
	    for line in f:
	        if not done:
	            if line[0] == ">":
	                if not skip:
	                    output.write(cluster)
	                if i == length:
	                    done = 1
	                    output.write(line)
	                    next
	                cluster = line
	                skip = 0
	            elif not skip:
	                if to_remove[i] == line.split(" >")[1].split("...")[0]:
	                    i += 1
	                    skip = 1
	                else:
	                    cluster += line
	        else:
	            output.write(line)

	if not skip:
	    output.write(cluster)
	output.close()
	end = time()
	print "finished in {} seconds".format(end-start)
	print "==========================================================================================="
