#!/bin/bash

#run this script for each sample

$1 = file 1 of paired ended fastq files
$2 = file 2 of paired ended fastq files
$3 = base sample name (i.e. SRR1234567)

megahit --mem-flag 2 -1 $1 -2 $2 -o $3_mh_out
prokka --outdir $3_prokka_out --addgenes --notbl2asn --metagenome --cpus 0 --mincontiglen 1 $3_mh_out/final.contigs.fa 
