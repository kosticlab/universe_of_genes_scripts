#!/bin/bash

$1 -n "5" -i "gene_catalog_filtered_sequences_translated_protein" -T 0 -M 0 -o "gene_catalog_filtered_sequences_translated_protein_100_perc" -c 1 -s .9 -aS .9

name="gene_catalog_filtered_sequences_translated_protein"
name1="gene_catalog_filtered_sequences_translated_protein_100_perc"
word=5
for i in $(seq 100 -5 50);

do

if [ "$i" -ge 70 -a "$i" -le 100 ];
then word=5
fi
if [ "$i" -ge 60 -a "$i" -lt 70 ];
then word=4
fi
if [ "$i" -ge 50 -a "$i" -lt 60 ];
then word=3
fi
if [ "$i" -ge 40 -a "$i" -lt 50 ];
then word=2
fi

$1 -n "$word" -i "$name1" -T 0 -M 0 -o "$name"_"$i"_perc -c ."$i" -s .9 -aS .9

name1="$name"_"$i"_perc

done

mkdir iterative_cdhit_output
mv gene_catalog_filtered_sequences_translated_protein*perc* iterative_cdhit_output
mv *95*filtered* iterative_cdhit_output
