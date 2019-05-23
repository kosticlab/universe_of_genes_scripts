#prep sequences for cluster filtering

#argument1: cluster file, argument2: list of genes to remove

grep '>' $2 > genes_to_remove

sed -i 's/>//g' genes_to_remove

grep "*" $1 | cut -d ">" -f 2 | cut -d "." -f 1 > consensus_genes.txt

awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' genes_to_remove consensus_genes.txt > unaligned_sort.txt
