#!/bin/bash
# `-wb` Write the original entry in B for each overlap in A (Genes)
for GROUP in f_vs_i f_vs_p f_vs_y i_vs_y p_vs_i y_vs_p; do
	for FEATURE in *.PparFemVer2024.transcriptome.6col.bed;
	do
	bedtools intersect -wb \
	    -a $FEATURE \
	    -b $GROUP.sorted.default_DMR_combined.bed > wb.DMRs.$GROUP.$FEATURE.tmp
	done
done 
