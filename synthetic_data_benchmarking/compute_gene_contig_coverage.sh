#!/bin/bash
echo 'Computing gene and contig coverage...'

id=$1

mkdir -p ${id}/${id}_contig_cov
mkdir -p ${id}/${id}_gene_cov

echo "raw reads in pair_1" >> ${id}/reads_count.txt && grep "@" $2 | wc -l >> ${id}/reads_count.txt
echo "raw reads in pair_2" >> ${id}/reads_count.txt && grep "@" $3 | wc -l >> ${id}/reads_count.txt
#Align reads with bbwrap.sh to contigs
#input = prokka output fna/ffn file (contig/gene), raw sequences pair 1, raw sequences pair 2
#output = aln.sam.gz, log file bbwrap.log
bbmap/bbwrap.sh ref=$(echo "${id}"/*.fna) in=$2 in2=$3 \
out=${id}/"${id}_contig_cov"/aln.sam.gz kfilter=22 subfilter=15 maxindel=80 2> ${id}/"${id}_contig_cov"/bbwrap.log
rm -r ./ref/
bbmap/bbwrap.sh ref=$(echo "${id}"/*.ffn) in=$2 in2=$3 \
out=${id}/"${id}_gene_cov"/aln.sam.gz kfilter=22 subfilter=15 maxindel=80 2> ${id}/"${id}_gene_cov"/bbwrap.log
rm -r ./ref/

#Output per contig coverage to cov.txt with pileup.sh
#input = aln.sam.gz
#output = coverage table cov.txt, log file coverage_summary.txt
bbmap/pileup.sh in=${id}/"${id}_contig_cov"/aln.sam.gz out=${id}/"${id}_contig_cov"/cov.txt 2>${id}/"${id}_contig_cov"/coverage_summary.txt
bbmap/pileup.sh in=${id}/"${id}_gene_cov"/aln.sam.gz out=${id}/"${id}_gene_cov"/cov.txt 2>${id}/"${id}_gene_cov"/coverage_summary.txt
