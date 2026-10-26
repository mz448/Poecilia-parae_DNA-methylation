#!/bin/bash
# Working folder
WORKDIR="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/data/features"

cd $WORKDIR

# Copy the GFF3 file from Ehren's folder
# cp /local/storage/Projects/EJB/P_parae_genome/P_parae_female_gene_models.gff3 /local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/

# copy the transcriptome into the working folder
# Notice that the transcriptome file was renamed to PparFemVer2024.transcriptome.gff3
cp /local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/P_parae_female_gene_models.gff3 $WORKDIR/PparFemVer2024.transcriptome.gff3
