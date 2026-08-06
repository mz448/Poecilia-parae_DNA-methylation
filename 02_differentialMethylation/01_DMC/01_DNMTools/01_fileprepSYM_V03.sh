#!/bin/bash
#SBATCH --job-name=fileprepSYM
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=30
#SBATCH --mem=300G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./%x_%j.out
#SBATCH -e ./%x_%j.err
# NOTE: %x_%j is going to be the name of the job

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# AUTHOR: MZF
# DATE: 20260718
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTES: 
#   This version has commented progress for clear progress documentation via 
#   logs/fileprepSYM.out file
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

#set DNMTools environment 
export LD_LIBRARY_PATH=/programs/htslib-1.16/lib 
export PATH=/programs/dnmtools-1.2.2/bin:$PATH 

date
echo 'Preparing CpG reports for DNMTools analysis...'

SCRIPTS="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/scripts/01_DNMTools"
WORKDIR="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/data"
cd "$WORKDIR" || exit 1


# Runing the custom python script against the CpGreports to generate .meth files
echo '↳ generating .meth files....'

for INPUT in *.CpG_report.txt;
	do
	OUTPUT="${INPUT}.meth"
	python3 ${SCRIPTS}/CpG_report-to-.meth_V2.py "$INPUT" "$OUTPUT"
	echo "DONE: ${OUTPUT} [✓]"
	done &&


# Removing flt() in $6 from python-crafted .meth files
for FILE in *.meth;
	do
	awk -F'[\t|.]' '{print $1 "\t" $2 "\t" $3 "\t" "CpG" "\t" $5 "." $6 "\t" $7}' $FILE > forSym.$FILE
	done &&

echo '=======Files converted to .meth format [✓]======='


#Run command
# dnmtools

# Running DNMTools SYM
echo 'Running DNMTools SYM...'

for FILE in forSym.ppar*;
do
    echo "Processing: ${FILE} ..."
    dnmtools sym -o "Sym.${FILE}" "$FILE"
    echo "DONE: Sym.${FILE} [✓]"
done &&

echo '=======DNMTools SYM completed [✓]======='


# Make the proportion tables
echo 'Making Proportion Tables....'
dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth -o proportion-table_f-vs-p.txt &&
echo "DONE: f-vs-p [✓]"
dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth  Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth -o proportion-table_f-vs-y.txt &&
echo "DONE: f-vs-y [✓]"
dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth -o proportion-table_f-vs-i.txt &&
echo "DONE: f-vs-i [✓]"
dnmtools merge -t -ignore -radmeth Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth -o proportion-table_p-vs-y.txt &&
echo "DONE: p-vs-y [✓]"
dnmtools merge -t -ignore -radmeth Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth -o proportion-table_p-vs-i.txt &&
echo "DONE: p-vs-i [✓]"
dnmtools merge -t -ignore -radmeth Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth -o proportion-table_y-vs-i.txt
echo "DONE: y-vs-i [✓]"
echo '=======Generated proportion tables [✓]======='

echo 'Symetric .meth files are ready for DNMTools analysis [✓]'
date
