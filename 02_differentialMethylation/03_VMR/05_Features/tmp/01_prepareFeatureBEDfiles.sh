#!/bin/bash
#SBATCH --job-name=features
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# 1) Copy reference .GFF3 file from Ehren's folder
cp /local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/P_parae_female_gene_models.gff3 ./PparFemVer2024.transcriptome.gff3
# 1.1) Copy the fererence .fasta file
cp /local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/PparFemMitoLambdaPuc19Ver2024.fasta ./PparFemMitoLambdaPuc19Ver2024.fasta

# 2) Make an index using the transcriptome
# Check AGAT gff2gtf index functionality
singularity run /programs/agat-1.2.0/agat.sif agat --tools
# Run the following command to make the transcriptome index
singularity run --bind $PWD --pwd $PWD /programs/agat-1.2.0/agat.sif agat_convert_sp_gff2gtf.pl -gff PparFemVer2024.transcriptome.gff3 -o PparFemVer2024.transcriptome.gtf

# 3) Add bedops to the path
export PATH=/programs/bedops-2.4.35/bin:$PATH
# 2) Make a bed file from a gff3 file (this converts the different numbering system and breaks it down in columns that can be easily manipulated by awk)
gff2bed < PparFemVer2024.transcriptome.gtf >PparFemVer2024.transcriptome.bed


# 4) Reformat the file so it now contains all the same columns but in a more useful order
# PparFemVer2024.transcriptome.6col.bed --> Chr, Start, End, Feature, Strand, GeneID
awk -F'[\t|;]' '{print $1 "\t" $2 "\t" $3 "\t" $8 "\t" $6 "\t" $10}' PparFemVer2024.transcriptome.bed > PparFemVer2024.transcriptome.6col.bed


# 5) Check the kinds of annotations available in the bed file
# Generate an output file with a summary of the features available in the transcriptome
awk '{print $4}' PparFemVer2024.transcriptome.6col.bed | sort | uniq -c | awk '{print $2 "\t" $1}' > featureSummary.txt &&
echo feature$'\t'count > header.txt &&
cat header.txt featureSummary.txt > featureSummary.PparFemVer2024.transcriptome.6col.txt &&
rm header.txt && rm featureSummary.txt


# This is how to filter by strand 
# awk -F'[\t|"]' '$4 =="gene" && $5=="+" {print $1 "\t" $3 "\t" $3+2000 "\t" $4 "\t" $5 "\t" $7}' PparFemVer2024.transcriptome.6col.bed > d2000.PparFemVer2024.transcriptome.6col.bed


# 6) make individual bed files for each feature available
# FEATURE.PparFemVer2024.transcriptome.6col.bed --> Chr, Start, End, Feature, Strand, GeneID
for FEATURE in CDS exon five_prime_UTR gene mRNA three_prime_UTR;
	do
	awk -v "FEATURE=$FEATURE" -F'[\t|"]' '$4 ==FEATURE {print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $4 "\t" $5}' PparFemVer2024.transcriptome.6col.bed > $FEATURE.PparFemVer2024.transcriptome.6col.bed
	done

# 7) Make new files for new features
# 7.1) UTRs
cat five_prime_UTR.PparFemVer2024.transcriptome.6col.bed three_prime_UTR.PparFemVer2024.transcriptome.6col.bed | sort -k1,1V -k2,2n > UTR.PparFemVer2024.transcriptome.6col.bed

# add BedTools to the path
export PATH=/programs/bedtools2-2.29.2/bin:$PATH 

# 7.2) intron <== gene - exon
bedtools subtract -s \
	    -a gene.PparFemVer2024.transcriptome.6col.bed \
	    -b exon.PparFemVer2024.transcriptome.6col.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "intron" "\t" $6}' > intron.PparFemVer2024.transcriptome.6col.bed
	    

# To check for the source of the count difference between exon and intron count 
# awk '{print $4}' exon.PparFemVer2024.transcriptome.6col.bed | sort -k1,1V -k2,2n| uniq -c | wc -l 
# awk '{print $4}' exon.PparFemVer2024.transcriptome.6col.bed | sort -k1,1V -k2,2n| uniq -c | awk '$1=="1" {print $0}' | wc -l
# awk '{print $4}' intron.PparFemVer2024.transcriptome.6col.bed | sort -k1,1V -k2,2n | uniq -c | wc -l

# 7.3) geneBody <== gene - UTRs
bedtools subtract -s \
	    -a gene.PparFemVer2024.transcriptome.6col.bed \
	    -b UTR.PparFemVer2024.transcriptome.6col.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "geneBody" "\t" $6}' > geneBody.PparFemVer2024.transcriptome.6col.bed

# 7.4) 2000bp before and after each geneBody

# 7.4.1.1)  u2000 <== Upstream 2000pb from genebody 
awk -F'[\t]' '$5 =="geneBody" && $6=="+" {print $1 "\t" $2-2000 "\t" $2 "\t" $4 "\t" "u2000" "\t" $6}' geneBody.PparFemVer2024.transcriptome.6col.bed > u2000.PparFemVer2024.transcriptome.6col.bed &&
awk -F'[\t]' '$5 =="geneBody" && $6=="-" {print $1 "\t" $3 "\t" $3+2000 "\t" $4 "\t" "u2000" "\t" $6}' geneBody.PparFemVer2024.transcriptome.6col.bed >> u2000.PparFemVer2024.transcriptome.6col.bed &&
sort -k1,1V -k2,2n u2000.PparFemVer2024.transcriptome.6col.bed > sorted.u2000.PparFemVer2024.transcriptome.6col.bed &&
mv sorted.u2000.PparFemVer2024.transcriptome.6col.bed u2000.PparFemVer2024.transcriptome.6col.bed

# 7.4.1.2) Make sure there are not out of-frame coordinates
awk -F'[\t]' '$2 > 0 {print $0}' u2000.PparFemVer2024.transcriptome.6col.bed > temp.u2000 &&
mv temp.u2000 u2000.PparFemVer2024.transcriptome.6col.bed

# 7.4.2)   d2000 <== Downstream 2000pb from genebody 
awk -F'[\t]' '$5 =="geneBody" && $6=="-" {print $1 "\t" $2-2000 "\t" $2 "\t" $4 "\t" "d2000" "\t" $6}' geneBody.PparFemVer2024.transcriptome.6col.bed > d2000.PparFemVer2024.transcriptome.6col.bed &&
awk -F'[\t]' '$5 =="geneBody" && $6=="+" {print $1 "\t" $3 "\t" $3+2000 "\t" $4 "\t" "d2000" "\t" $6}' geneBody.PparFemVer2024.transcriptome.6col.bed >> d2000.PparFemVer2024.transcriptome.6col.bed &&
sort -k1,1V -k2,2n d2000.PparFemVer2024.transcriptome.6col.bed > sorted.d2000.PparFemVer2024.transcriptome.6col.bed &&
mv sorted.d2000.PparFemVer2024.transcriptome.6col.bed d2000.PparFemVer2024.transcriptome.6col.bed

# 7.5) GeneBody-intron <== geneBody - CDS
bedtools subtract -s \
	    -a gene.PparFemVer2024.transcriptome.6col.bed \
	    -b CDS.PparFemVer2024.transcriptome.6col.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "geneBody_intron" "\t" $6}' > geneBody_intron.PparFemVer2024.transcriptome.6col.bed

# 7.6) GeneBody-exon <== CDS
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "geneBody_exon" "\t" $6}' CDS.PparFemVer2024.transcriptome.6col.bed > geneBody_exon.PparFemVer2024.transcriptome.6col.bed

# 7.7) UTR-exon <== UTRs ∩ exon
bedtools intersect -s \
    -a UTR.PparFemVer2024.transcriptome.6col.bed \
    -b exon.PparFemVer2024.transcriptome.6col.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "_exon" "\t" $6}' > UTR_exon.PparFemVer2024.transcriptome.6col.bed

# 7.7.1) Make individual files for 5' and 3' UTR_exon
awk -F'[\t]' '$5 == "five_prime_UTR_exon" {print $0}' UTR_exon.PparFemVer2024.transcriptome.6col.bed > five_prime_UTR_exon.PparFemVer2024.transcriptome.6col.bed && 
awk -F'[\t]' '$5 == "three_prime_UTR_exon" {print $0}' UTR_exon.PparFemVer2024.transcriptome.6col.bed > three_prime_UTR_exon.PparFemVer2024.transcriptome.6col.bed
