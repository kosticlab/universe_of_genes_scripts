###run synthetic data modeling
##20190719
#Tierney

import os
import random
import subprocess
import numpy as np
import pandas as pd
import sys
import time
###outside requirements:
#bbmap in working directory
#diamond in working directory
#megahit in working directory
#prokka in working directory
#metaspades in working directory
#folder called refgenomes containing readdata and orfs folders with appropriate information inside (eg raw reads generated with ART with corresponding orf files)
#compute_gene_contig_coverage.sh script, executable and in working directory
#seqtk in path

def create_synthetic_metagenome(filesToUse,orfList,seqFiles,iteration): #CHECKED
	print 'Creating synthetic metagenomic data...'
	os.system('rm synthetic_metagenome* known_genes tempFile_* 2> error.log')
	for f in filesToUse:
		coverage=random.uniform(covMin,covMax)
		coverage=round(coverage,2)
		sampleFiles=[x for x in seqFiles if f in x]
		lowCov=False
		# if low coverage tempfile present (coverage less than one), cat with other data and clean up extra files
		orfFile=[x for x in orfList if f in x][0]
		os.system('cat ./refgenomes/orfs/%s >> known_genes'%orfFile)
		if coverage<1:
			lowCov=True
			totalLength=subprocess.check_output("grep '@' ./refgenomes/readdata/%s | wc -l"%sampleFiles[0],shell=True).rstrip()
			x=int(round(float(totalLength)*coverage))
			seed=random.randint(0, 100)
			os.system("seqtk sample -s%s ./refgenomes/readdata/%s %s > tempFile_1_%s "%(seed,sampleFiles[0],x,sampleFiles[0]))
			os.system("seqtk sample -s%s ./refgenomes/readdata/%s %s > tempFile_2_%s "%(seed,sampleFiles[1],x,sampleFiles[1]))
			continue
		else:
			coverage=round(coverage)
		for i in range(0,int(coverage)):
			#create metagenome with coverage estimate
			os.system('cat ./refgenomes/readdata/%s >> synthetic_metagenome_1.fq'%sampleFiles[0])
			os.system('cat ./refgenomes/readdata/%s >> synthetic_metagenome_2.fq'%sampleFiles[1])
	if lowCov==True:
		os.system('cat tempFile_1_* synthetic_metagenome_1.fq > foo ; mv foo synthetic_metagenome_1.fq ; rm tempFile_1_*')
		os.system('cat tempFile_2_* synthetic_metagenome_2.fq > foo ; mv foo synthetic_metagenome_2.fq ; rm tempFile_2_*')
	#get all genes
	os.system('rm tempFile_* 2>> error.log')

def assemble(m1,m2): # CHECKED
	print 'Assembling with megahit (default)...'
	os.system('./megahit/megahit --mem-flag 2 -1 ./%s -2 ./%s -o temp_mh_out_default 2>> error.log'%(m1,m2))
	print 'Assembling with megahit (sensitive)...'
	os.system('./megahit/megahit --presets meta-sensitive --mem-flag 2 -1 ./%s -2 ./%s -o temp_mh_out_sensitive 2>> error.log'%(m1,m2))
	print 'Assembling with megahit (large)...'
	os.system('./megahit/megahit --presets meta-large --mem-flag 2 -1 ./%s -2 ./%s -o temp_mh_out_large 2>> error.log'%(m1,m2))
	print 'Assembling with metaspades...'
	os.system('python SPAdes-3.13.0-Linux/bin/metaspades.py --only-assembler -1 ./%s -2 ./%s -o temp_ms_out 2>> error.log'%(m1,m2))
	os.system('mkdir assembly_output ; mv temp_mh_out_default/final.contigs.fa assembly_output/megahit_contigs_default ; mv temp_mh_out_sensitive/final.contigs.fa assembly_output/megahit_contigs_sensitive ; mv temp_mh_out_large/final.contigs.fa assembly_output/megahit_contigs_large ; mv temp_ms_out/contigs.fasta assembly_output/metaspades_contigs ; rm -r temp_mh_out* temp_ms_out')

def run_prokka(): #CHECKED
	samples=os.listdir('assembly_output')
	for f in samples:
		print 'Running prokka for %s...'%f
		os.system("./prokka/bin/prokka --outdir prokka_output_%s --addgenes --metagenome --cpus 0 --mincontiglen 1 assembly_output/%s 2> prokka.log"%(f,f))
	os.system('rm -rf assembly_output')

def false_positive_check_1(prokkaOutput): #CHECKED
	print 'Finding false positive genes (method 1)...'
	os.system('./diamond makedb --in %s --db tempdb'%prokkaOutput)
	os.system('./diamond blastx -q ./known_genes --db ./tempdb -o tempalignments --max-target-seqs 1 --id .95')
	output=[]
	with open('tempalignments') as f:
		for line in f:
			output.append(line.rstrip().split('\t')[1])
	output=set(output)
	output2=[]
	with open(prokkaOutput) as f:
		for line in f:
			if line[0]=='>':
				output2.append(line.rstrip()[1:].split(' ')[0])
	output3=list(set(output2)-set(output))
	falsepositive=len(output3)
	os.system('rm tempalignments tempdb*')
	return [falsepositive,output3]

def find_singleton_count(clusterFile): #CHECKED
	count=0
	singletonCount=0
	output=[]
	with open(clusterFile) as f:
		for line in f:
			if 'Cluster' in line:
				if count==1:
					singletonCount+=1
					output.append(geneid)
				count=0
			else:
				geneid=line.split('>')[1].split('...')[0]
				count+=1
	if count==1:
		singletonCount+=1
		output.append(line.split('>')[1].split('...')[0])
		count=0
	return [singletonCount,set(output)]

def false_positive_check_2(prokkaOutput): #CHECKED
	print 'Finding false positive genes (method 2)...'
	os.system('./diamond blastx --db uniref100.dmnd -q %s  -k 1 -o tempalignments --un diamond_unaligned 2> diamond.log'%(prokkaOutput))
	falsepositive=subprocess.check_output("grep '>' diamond_unaligned  | wc -l",shell=True).rstrip()
	output=[]
	with open('diamond_unaligned') as f:
		for line in f:
			if '>' in line:
				output.append(line[1:].rstrip())
	os.system('rm diamond_unaligned tempalignments diamond.log')
	return [falsepositive,output]

def map_gene_contig(gfffile): #CHECKED
	mapping={}
	with open(gfffile) as f:
		for line in f:
				if 'ID' in line:
					id=line.split('ID=')[1].split(';')[0].split('_gene')[0].split(' ')[0]
					contig=line.split('\t')[0]
					mapping[id]=contig
	return mapping

def parse_fasta_lengths_to_dict(fasta): #CHECKED
	#d={}
	d2={}
	lines=[]
	with open(fasta) as f:
		for i,line in enumerate(f):
			if '>' in line:
				if i == 0:
					idval = line.rstrip().replace('>','').split(' ')[0]
					continue
				seq=''.join(lines)
				#d[id]=seq
				d2[idval]=len(seq)
				lines=[]
				idval = line.rstrip().replace('>','').split(' ')[0]
			else:
				lines.append(line.rstrip())
		#d[id]=''.join(lines)
		d2[idval]=len(''.join(lines))
	return d2

def get_contig_lengths(contigFile): #CHECKED
	print 'Getting contig lengths...'
	lengthDict=parse_fasta_lengths_to_dict(contigFile)
	meanLength=np.mean(lengthDict.values())
	lengthStd=np.std(lengthDict.values())
	return [meanLength,lengthStd,lengthDict]

def get_gene_lengths_and_contig_mapping(fasta,contigLengthDict,geneContigMapping):
	print 'Getting gene lengths...'
	lengthDict = parse_fasta_lengths_to_dict(fasta)
	output=[]
	for val in lengthDict.keys():
		contigId=geneContigMapping[val]
		contigLength=contigLengthDict[contigId]
		geneLength=lengthDict[val]
		output.append([val,geneLength,contigId,contigLength])
	meanLength=np.mean([x for x in lengthDict.values()])
	geneStd=np.std([x for x in lengthDict.values()])
	return [output,meanLength,geneStd]

def process_output(prokkaOutput,covMin,covMax,iteration):
	print 'Processing prokka output...'
	summaryOutput=[]
	summaryOutputNames=[]
	summaryOutput.append(subprocess.check_output("grep '@' synthetic_metagenome_1.fq | wc -l",shell=True).rstrip())
	summaryOutputNames.append('Raw Reads in PE 1')
	summaryOutput.append(subprocess.check_output("grep '@' synthetic_metagenome_2.fq | wc -l",shell=True).rstrip())
	summaryOutputNames.append('Raw Reads in PE 2')
	summaryOutput.append(subprocess.check_output("grep '>' known_genes | wc -l",shell=True).rstrip())
	summaryOutputNames.append('Ground Truth Gene Count')
	# count genes in genomes (to later compare to prokka output)
	for p in prokkaOutput:
		print 'Processing %s...'%p
		p_sub=os.listdir(p)
		p2=[x for x in p_sub if 'ffn' in x][0]
		#count genes in prokka output
		summaryOutput.append(subprocess.check_output("grep '>' %s/%s | wc -l"%(p,p2),shell=True).rstrip())
		summaryOutputNames.append('Prokka Gene Count %s'%p)
		#compare prokka output to ground truth using cd-hit clustering
		fp1,fp1Ids=false_positive_check_1(p+'/'+p2.split('.')[0]+'.faa')
		summaryOutput.append(fp1)
		summaryOutputNames.append('False Positives (Compared to Ground Truth) %s'%p)
		#compare prokka output to nr to find genes never identified before
	#	fp2,fp2Ids=false_positive_check_2(p+'/'+p2)
	#	summaryOutput.append(fp2)
	#	summaryOutputNames.append('False Positives (Compared to NR database) %s'%p)
		#record length of each contig
		contigMean,contigStd,contigLengthDict=get_contig_lengths(p+'/'+p2.split('.')[0]+'.fna')
		summaryOutput.append(contigMean)
		summaryOutputNames.append('Mean Contig Length %s'%p)
		summaryOutput.append(contigStd)
		summaryOutputNames.append('Contig Length Standard Deviation %s'%p)
		geneContigMapping=map_gene_contig(p+'/'+p2.split('.')[0]+'.gff')
		#find gene lengths, map genes to contigs
		geneLevelData,geneMeanLength,geneStd=get_gene_lengths_and_contig_mapping(p+'/'+p2,contigLengthDict,geneContigMapping)
		summaryOutput.append(geneMeanLength)
		summaryOutputNames.append('Mean Gene Length %s'%p)
		summaryOutput.append(geneStd)
		summaryOutputNames.append('Gene Length Standard Deviation %s'%p)
		#merge gene data with false positive data
		for i,x in enumerate(geneLevelData):
			geneId=x[0]
			if geneId in fp1Ids:
				x.append(1)
			if geneId not in fp1Ids:
				x.append(0)
		#	if geneId in fp2Ids:
		#		x.append(1)
		#	if geneId not in fp2Ids:
		#		x.append(0)
			geneLevelData[i]=x
		output=pd.DataFrame(geneLevelData)
		output.columns=['Gene ID','Gene Length','Contig ID','Contig Length','False Positive (method 1)']#,'False Positive (method 2)']
		# compute coverage of each gene and of each contig
		print 'Calculating contig-by-contig coverage...'
		os.system('./compute_gene_contig_coverage.sh %s synthetic_metagenome_1.fq synthetic_metagenome_2.fq'%(p))
		#merge in coverage summary/specific data
		contigData=pd.read_csv('%s/%s_contig_cov/cov.txt'%(p,p),sep='\t',index_col=0)
		contigData=contigData[['Avg_fold','Covered_percent']]
		contigData.columns=['Avg_fold_contig','Covered_percent_contig']
		geneData=pd.read_csv('%s/%s_gene_cov/cov.txt'%(p,p),sep='\t',index_col=0)
		geneData=geneData[['Avg_fold','Covered_percent']]
		geneData.columns=['Avg_fold_gene','Covered_percent_gene']
		geneData.index=[x.split(' ')[0] for x in geneData.index]
		print 'Cleaning and writing output to file...'
		output=output.merge(contigData,left_on=output.columns[2],right_index=True)
		output=output.merge(geneData,left_on=output.columns[0],right_index=True)
		# save output
		output.to_csv('coverage_%s_%s/%s/%s_gene_by_gene.csv'%(covMin,covMax,iteration,p))
		with open('%s/%s_contig_cov/coverage_summary.txt'%(p,p)) as f:
			for line in f:
				if 'Percent mapped' in line or 'Average coverage' in line or 'Percent scaffolds with any coverage' in line:
					summaryOutput.append(line.split(' ')[-1].strip())
		summaryOutputNames.append('Contig Percent Mapped %s'%p)
		summaryOutputNames.append('Contig Average Coverage %s'%p)
		summaryOutputNames.append('Percent Contigs with Any Coverage %s'%p)
		with open('%s/%s_gene_cov/coverage_summary.txt'%(p,p)) as f:
			for line in f:
				if 'Percent mapped' in line or 'Average coverage' in line or 'Percent scaffolds with any coverage' in line:
					summaryOutput.append(line.split(' ')[-1].strip())
		summaryOutputNames.append('Gene Percent Mapped %s'%p)
		summaryOutputNames.append('Gene Average Coverage %s'%p)
		summaryOutputNames.append('Percent Genes with Any Coverage %s'%p)
	summaryInfoAll=pd.DataFrame(summaryOutput).transpose()
	summaryInfoAll.columns=summaryOutputNames
	summaryInfoAll.to_csv('coverage_%s_%s/%s/summary_data.csv'%(covMin,covMax,iteration))
	os.system('mv prokka_output* coverage_%s_%s/%s/ ; mv synthetic_metagenome* coverage_%s_%s/%s/ ; mv known_genes coverage_%s_%s/%s/'%(covMin,covMax,iteration,covMin,covMax,iteration,covMin,covMax,iteration))

def main():
	for i in range(10):
		start=time.time()
		iteration=i
		os.system('mkdir coverage_%s_%s 2> error.log'%(covMin,covMax))
		os.system('mkdir coverage_%s_%s/%s'%(covMin,covMax,iteration))
		# pick 95 genomes randomly
		fileNames=os.listdir('./refgenomes/readdata')
		fileNames=list(set(['_'.join(x.split('_fullSeq')[:-1]) for x in fileNames]))
		random.shuffle(fileNames)
		filesToUse=fileNames[:95]
		seqFiles=os.listdir('./refgenomes/readdata')
		orfList=os.listdir('./refgenomes/orfs')
		# generate a random number for coverage value and, based on that number, combine genomes in correct proportions to create synthetic metagenome
		create_synthetic_metagenome(filesToUse,orfList,seqFiles,iteration)
		# Run assembler on ART output
		assemble('synthetic_metagenome_1.fq','synthetic_metagenome_2.fq')
		# Run prokka out assembled contigs
		run_prokka()
		#find prokka output
		prokkaOutput=os.listdir('.')
		prokkaOutput=[x for x in prokkaOutput if 'prokka_output' in x]
		process_output(prokkaOutput,covMin,covMax,iteration)
		end=time.time()
		print 'Iteration %s took %s seconds.'%(str(i),str(round(end-start,2)))

if __name__ == '__main__':
	covMin,covMax=int(sys.argv[1]),int(sys.argv[2])
	main()



#iterations: 10 "normal" distribution of coverage (10 to 30X), 10 medium coverage (1 to 20X), 10 low coverage (0 to 10X), 10 extremely low coverage (0 to 1X) , 10 broad coverage (0, 30X)
