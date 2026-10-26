#!/bin/bash
# export PATH=/programs/bedtools2-2.29.2/bin:$PATH 
# `-wb` Write the original entry in B for each overlap in A (Genes)
for GROUP in f_vs_i f_vs_p f_vs_y i_vs_y p_vs_i y_vs_p; do
	for FEATURE in CDS d2000 exon five_prime_UTR_exon five_prime_UTR geneBody_exon geneBody_intron geneBody gene intron mRNA three_prime_UTR_exon three_prime_UTR u2000 UTR_exon UTR;
	do
	bedtools intersect -wb \
	    -a $FEATURE.PparFemVer2024.transcriptome.6col.bed \
	    -b $GROUP.sorted.default_DMR_combined.bed > wb.DMRs.$GROUP.$FEATURE.tmp
	done
done 
