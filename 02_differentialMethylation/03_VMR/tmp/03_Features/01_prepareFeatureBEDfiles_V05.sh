#!/bin/bash
# DATE:  	  20250825 
# AUTHOR: 	MZF
# SCRIPT: 	01_prepareFeatureBEDfiles_05.sh
#
# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL: 	Extract default features and custom features from gff3 file 

# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 1) Make a bed file form the gtf and format it 
# Add bedops to the path
export PATH=/programs/bedops-2.4.35/bin:$PATH
# Make a bed file from a gft file (this converts the different numbering system and breaks it down in columns that can be easily manipulated by awk)
gff2bed < P_parae_Female.gtf > PparFemVer2024.transcriptome.bed


# Reformat the file so it now contains all the same columns but in a more useful order
# PparFemVer2024.transcriptome.6col.bed --> Chr, Start, End, Feature, Strand, GeneID
awk -F'[\t|;]' '{print $1 "\t" $2 "\t" $3 "\t" $8 "\t" $6 "\t" $10}' PparFemVer2024.transcriptome.bed | awk -F'[\t|"]' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}'  > PparFemVer2024.transcriptome.6col.bed

# Check the kinds of annotations available in the bed file
# Generate an output file with a summary of the features available in the transcriptome
awk '{print $4}' PparFemVer2024.transcriptome.6col.bed | sort | uniq -c | awk '{print $2 "\t" $1}' > featureSummary.txt &&
echo feature$'\t'count > header.txt &&
cat header.txt featureSummary.txt > featureSummary.PparFemVer2024.transcriptome.6col.txt.tmp &&
rm header.txt && rm featureSummary.txt

echo "[✓] 1: Made reference file & and listed available features "

# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 2) Make individual bed files for each feature available
# FEATURE.PparFemVer2024.transcriptome.6col.bed --> Chr, Start, End, Feature, Strand, GeneID
for FEATURE in CDS exon five_prime_UTR gene mRNA three_prime_UTR; do
  awk -v FEATURE="$FEATURE" -F'\t' \
      '$4==FEATURE {print $1,$2,$3,$6,$4,$5}' OFS="\t" \
      PparFemVer2024.transcriptome.6col.bed \
  > ${FEATURE}.PparFemVer2024.transcriptome.6col.bed
done

echo "[✓] 2: Made individual bed files for each feature available"

# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 3) De-duplicate the base feature BEDs (safe)
for F in gene exon CDS five_prime_UTR three_prime_UTR; do
  bedtools sort -i ${F}.PparFemVer2024.transcriptome.6col.bed \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' \
  | sort -u -k1,1V -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 \
  > dedup.${F}.6col.bed
done

# UTR = 5' ∪ 3' (deduped)
# De-duplicates identical rows and keeps one copy. 
# The sort keys are:
# -k1,1V = “version” sort on chr (so Parae_2 < Parae_10, not lexicographic)
# -k2,2n = numeric start
# -k3,3n = numeric end
# -k4,4 = name (gene ID)
# -k5,5 = score/label (I am using the feature label here in some steps)
# -k6,6 = strand
# Using -u with this full key means any perfect duplicates (same coords, gene, label, strand) collapse to one line. This is the safe fix for cases where multiple transcripts or earlier steps create exact duplicates. (Importantly, this is not uniq -u, which would drop all duplicates entirely.)

cat dedup.five_prime_UTR.6col.bed dedup.three_prime_UTR.6col.bed \
| bedtools sort -i - \
| sort -u -k1,1V -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 \
> dedup.UTR.6col.bed

echo "[✓] 3: Deduplicated available features "

# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 4) Make `geneBody` and Intron using **gene-aware** subtraction
# Gene-aware intron = gene − exon (same gene only)

# Tag chrom with "|geneID" for both genes and exons, then sort (bedtools needs sorted)
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,0,$6}' dedup.gene.6col.bed \
  | bedtools sort -i - > tag.gene.bed

awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,0,$6}' dedup.exon.6col.bed \
  | bedtools sort -i - > tag.exon.bed

# Subtract on the tagged space (so only same-gene exons subtract), then untag chrom
bedtools subtract -s -a tag.gene.bed -b tag.exon.bed \
| awk 'BEGIN{OFS="\t"}{split($1,a,"|"); print a[1],$2,$3,$4,"intron",$6}' \
| bedtools sort -i - \
> intron.PparFemVer2024.transcriptome.6col.bed


# To restrict subtraction to UTRs from the **same gene**, we temporarily tag the chromosome with the gene ID
# in both files so only same-gene records can subtract each other.

# Tag chrom with |geneID, SORT, subtract, then untag
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,0,$6}' dedup.gene.6col.bed \
  | bedtools sort -i - > tag.gene.bed
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,0,$6}' dedup.UTR.6col.bed  \
  | bedtools sort -i - > tag.UTR.bed

bedtools subtract -s -a tag.gene.bed -b tag.UTR.bed \
| awk 'BEGIN{OFS="\t"}{split($1,a,"|"); print a[1],$2,$3,$4,"geneBody",$6}' \
| bedtools sort -i - \
> geneBody.PparFemVer2024.transcriptome.6col.bed

echo "[✓] 4: Calculated geneBodies and Introns using gene awarnes"

# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 5) Compute the intragenic tracks
# geneBody_exon = CDS ∩ geneBody  (safe guard)
bedtools intersect -s \
  -a dedup.CDS.6col.bed \
  -b geneBody.PparFemVer2024.transcriptome.6col.bed -wa \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,"geneBody_exon",$6}' \
> geneBody_exon.PparFemVer2024.transcriptome.6col.bed

# geneBody_intron = geneBody − CDS  (GENE-AWARE)
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,0,$6}' geneBody.PparFemVer2024.transcriptome.6col.bed \
  | bedtools sort -i - > tag.geneBody.bed
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,0,$6}' dedup.CDS.6col.bed \
  | bedtools sort -i - > tag.CDS.bed
bedtools subtract -s -a tag.geneBody.bed -b tag.CDS.bed \
| awk 'BEGIN{OFS="\t"}{split($1,a,"|"); print a[1],$2,$3,$4,"geneBody_intron",$6}' \
| bedtools sort -i - \
> geneBody_intron.PparFemVer2024.transcriptome.6col.bed


echo "[✓] 5: Computed intragenic tracks"


# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 6) Compute fixed regions before and after the gene Body

# --- u2000 ---
# build
awk -F'\t' '$5=="geneBody" && $6=="+" {print $1, $2-2000, $2, $4, "u2000", $6}
            $5=="geneBody" && $6=="-" {print $1, $3,      $3+2000, $4, "u2000", $6}' \
  OFS="\t" geneBody.PparFemVer2024.transcriptome.6col.bed \
> u2000.PparFemVer2024.transcriptome.6col.bed

# clip starts, then sort
awk -F'\t' '$2>=0' u2000.PparFemVer2024.transcriptome.6col.bed > tmp.u2000 && mv tmp.u2000 u2000.PparFemVer2024.transcriptome.6col.bed
bedtools sort -i u2000.PparFemVer2024.transcriptome.6col.bed > tmp.u2000 && mv tmp.u2000 u2000.PparFemVer2024.transcriptome.6col.bed

# --- d2000 ---
# build
awk -F'\t' '$5=="geneBody" && $6=="-" {print $1, $2-2000, $2, $4, "d2000", $6}
            $5=="geneBody" && $6=="+" {print $1, $3,      $3+2000, $4, "d2000", $6}' \
  OFS="\t" geneBody.PparFemVer2024.transcriptome.6col.bed \
> d2000.PparFemVer2024.transcriptome.6col.bed

# clip starts, then sort
awk -F'\t' '$2>=0' d2000.PparFemVer2024.transcriptome.6col.bed > tmp.d2000 && mv tmp.d2000 d2000.PparFemVer2024.transcriptome.6col.bed
bedtools sort -i d2000.PparFemVer2024.transcriptome.6col.bed > tmp.d2000 && mv tmp.d2000 d2000.PparFemVer2024.transcriptome.6col.bed

echo "[✓] 6: Computed additional features (fixed length UTRs)"


# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 7) Make individual files for 5' and 3' UTR_exon
# 7) UTR_exon (gene-aware) + split into 5′/3′

# Tag by geneID and sort (bedtools likes sorted input)
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,$5,$6}' dedup.UTR.6col.bed  \
  | bedtools sort -i - > tag.utr.bed
awk 'BEGIN{OFS="\t"}{print $1"|"$4,$2,$3,$4,$5,$6}' dedup.exon.6col.bed \
  | bedtools sort -i - > tag.exon.bed

# Intersect within the same gene (strand-aware), then untag and relabel
bedtools intersect -s -wa -a tag.utr.bed -b tag.exon.bed \
| awk 'BEGIN{OFS="\t"}{
         split($1,a,"|"); chr=a[1];
         # $5 is either "five_prime_UTR" or "three_prime_UTR"
         print chr,$2,$3,$4,$5"_exon",$6
      }' \
| bedtools sort -i - \
| sort -u -k1,1V -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 \
> UTR_exon.PparFemVer2024.transcriptome.6col.bed

# Split into 5′ and 3′ (already sorted)
awk -F'\t' '$5=="five_prime_UTR_exon"'  UTR_exon.PparFemVer2024.transcriptome.6col.bed \
> five_prime_UTR_exon.PparFemVer2024.transcriptome.6col.bed

awk -F'\t' '$5=="three_prime_UTR_exon"' UTR_exon.PparFemVer2024.transcriptome.6col.bed \
> three_prime_UTR_exon.PparFemVer2024.transcriptome.6col.bed


echo "[✓] 7: Split UTR regions 5′/3 and calculate exons within them"


# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 8) Filter files to keep only unique instances of the features
# Base tracks (use dedup)
for F in CDS exon five_prime_UTR three_prime_UTR UTR gene mRNA; do
  SRC="dedup.${F}.6col.bed"
  [ -s "$SRC" ] || SRC="${F}.PparFemVer2024.transcriptome.6col.bed"
  bedtools sort -i "$SRC" | sort -u -k1,1V -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 \
  | awk 'BEGIN{OFS="\t"}{print $0, $3-$2}' > "${F}.feature.bed"
done

# Derived tracks (use freshly built files)
for F in intron geneBody geneBody_exon geneBody_intron \
         UTR_exon five_prime_UTR_exon three_prime_UTR_exon u2000 d2000; do
  bedtools sort -i "${F}.PparFemVer2024.transcriptome.6col.bed" 2>/dev/null \
  | sort -u -k1,1V -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 \
  | awk 'BEGIN{OFS="\t"}{print $0, $3-$2}' > "${F}.feature.bed" || true
done

echo "[✓] 8: Filter files to keep only unique instances of the features "

# ~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~%
# 9) clean files
echo "[...] cleaning files"
rm *.PparFemVer2024.transcriptome.6col.bed
rm dedup.*
rm tag.*
rm *.tmp

echo "[✓] Feature tracks ready! "
