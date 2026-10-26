#!/bin/bash

cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/data/dmr_all.sorted.withID.tsv .
# `-c` Reporting the number of DMRs overlapping per each feature
for GROUP in f_vs_i f_vs_p f_vs_y i_vs_y p_vs_i y_vs_p; do
	for FEATURE in *.PparFemVer2024.transcriptome.6col.bed;
	do
	bedtools intersect -c \
	    -a $FEATURE \
	    -b $GROUP.sorted.default_DMR_combined.bed > DMRsOvelapping.$GROUP.$FEATURE.tmp
	done
done 

# Filter by the parent genes [thanks to StackOverflow](https://stackoverflow.com/questions/23465478/summing-up-a-column-of-numbers-for-each-unique-item-in-another-column) 
# New 3 column file: [chr geneID CountOfDMRsOverlapping] 
for FILE in DMRsOvelapping*;
do 
awk '$4!=prev && NR>1 { print $1 "\t" prev "\t" sum; sum = 0;} {prev = $4; sum += $7} END {print $1 "\t" prev "\t" sum}' $FILE > countPerGene.$FILE
done 


# This is for filtering only the DMRs that overlap with Features in the + strand. however this approach is ignoring the logic of the simetric conversion, therefore should be avoided
# for FILE in DMRsOvelapping*;
# do 
# awk '$4!=prev && NR>1 && $6 == "+" { print $1 "\t" prev "\t" sum; sum = 0;} {prev = $4; sum += $7} END {print $1 "\t" prev "\t" sum}' $FILE > countPerGene.$FILE
# done 


# Report Per chromosome, how many genes are overlaping with at least one DMR
# New 2 column file: [chr CountOfDMRsOverlappingPerGene] 
for FILE in countPerGene.*;
do
awk '$3!="0" {print $1}' $FILE | uniq -c | awk '{print $2  "\t" $1}' > perChr.$FILE
done


# Make reference count of genes with a given feature
for FILE in CDS* d2000* exon* five_prime_UTR_exon* five_prime_UTR* geneBody_exon* geneBody_intron* geneBody* gene* intron* mRNA* three_prime_UTR_exon* three_prime_UTR* u2000* UTR_exon* UTR*;
do
awk -F '[\t]' '{print $1, $4}' $FILE | uniq | awk '{count[$1]++} END {for (chr in count) print chr "\t" count[chr]}' | sort -k1,1V > perChr.countOfGenes.reference.$FILE 
done
