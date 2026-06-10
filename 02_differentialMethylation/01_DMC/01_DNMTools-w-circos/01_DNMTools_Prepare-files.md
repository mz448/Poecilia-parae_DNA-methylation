____
**Creation Date:**  20240717 
**Modification date**: 2420240717 16:14:30
**Author(s):** Max 
_____
**GOAL :**
Get counts of methylations per **individual** & **per position** using the BISMARK alignment
____
# 0) Path from BISMARK alignments!
```bash
# Female
cd /local/storage/Projects/ppar_emseq/analysis/007_emseq_pparae_muscle/nfcoreMethylseq_Bismark_3/bismark/coverage2cytosine/reports/curated_reports
# Immaculata
cd /local/storage/Projects/ppar_emseq/analysis/008_emseq_pparae_muscle/Immac/nfcoreMethylseq_Bismark/bismark/coverage2cytosine/reports
# Parae
cd /local/storage/Projects/ppar_emseq/analysis/008_emseq_pparae_muscle/Parae/nfcoreMethylseq_Bismark/bismark/coverage2cytosine/reports
# Yellow
cd /local/storage/Projects/ppar_emseq/analysis/008_emseq_pparae_muscle/Yellow/nfcoreMethylseq_Bismark/bismark/coverage2cytosine/reports

```

# 1) Prepare files for DNMTools by crafting count files `.meth` from Bismark output
```bash

#set DNMTools environment 
export LD_LIBRARY_PATH=/programs/htslib-1.16/lib 
export PATH=/programs/dnmtools-1.2.2/bin:$PATH 

# Make a new directory inside the Bismark report folder `bismark/coverage2cytosine/reports/curated_reports/`
mkdir DNMTools &&

# Unzip samples
for ZIP in *.gz; do
gunzip $ZIP 
done &&

# Copy files from the parent folder
cp *CpG_report.txt DNMTools/. &&

# Get into the processing folder
cd  DNMTools &&

cp /local/storage/Projects/ppar_emseq/analysis/008_emseq_pparae_muscle/Immac/nfcoreMethylseq_Bismark/bismark/coverage2cytosine/reports/DNMTools/CpG_report-to-.meth_V2.py . &&

# Run a python script against the CpGreports to generate .meth files
for INPUT in *.CpG_report.txt; do 
OUTPUT="${INPUT}.meth" 
python3 CpG_report-to-.meth_V2.py "$INPUT" "$OUTPUT"
done &&

# Removing flt() in $6 from python-crafted .meth files
for FILE in *.meth; do 
awk -F'[\t|.]' '{print $1 "\t" $2 "\t" $3 "\t" "CpG" "\t" $5 "." $6 "\t" $7}' $FILE > forSym.$FILE
done &&

#Run command
# dnmtools

# Running DNMTools SYM
for FILE in  forSym.ppar* ; do
	dnmtools sym -o Sym.$FILE $FILE
	done &&

dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth  Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth -o proportion-table_f-vs-y.txt

```
# 2) Make proportion tables
```
# Make the proportion tables
dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth -o proportion-table_f-vs-p.txt &&

dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth  Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth -o proportion-table_f-vs-y.txt &&

dnmtools merge -t -ignore -radmeth Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth -o proportion-table_f-vs-i.txt &&

dnmtools merge -t -ignore -radmeth Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth -o proportion-table_p-vs-y.txt &&

dnmtools merge -t -ignore -radmeth  Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth -o proportion-table_p-vs-i.txt &&

dnmtools merge -t -ignore -radmeth Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth -o proportion-table_y-vs-i.txt
```


## Make the design-matrixes
>[!Tip]
>Only include the name of the files in the matrix, **exclude** the `.meth`
>Otherwise, DNMTools crashes
>This is because the Proportion-table does not use the the `.meth` suffix for referencing the files, but their names only.


```bash
# make files manually
nano FILE_NAME.txt

# design-matrix_f-vs-p.txt
base case 
Sym.forSym.pparfmem001.CpG_report.txt 1 0 
Sym.forSym.pparfmem002.CpG_report.txt 1 0 
Sym.forSym.pparfmem003.CpG_report.txt 1 0 
Sym.forSym.pparpmem004.CpG_report.txt 1 1
Sym.forSym.pparpmem005.CpG_report.txt 1 1
Sym.forSym.pparpmem006.CpG_report.txt 1 1

# design-matrix_f-vs-y.txt
base case 
Sym.forSym.pparfmem001.CpG_report.txt 1 0 
Sym.forSym.pparfmem002.CpG_report.txt 1 0 
Sym.forSym.pparfmem003.CpG_report.txt 1 0 
Sym.forSym.pparymem007.CpG_report.txt 1 1 
Sym.forSym.pparymem008.CpG_report.txt 1 1 
Sym.forSym.pparymem009.CpG_report.txt 1 1

# design-matrix_f-vs-i.txt
base case 
Sym.forSym.pparfmem001.CpG_report.txt 1 0 
Sym.forSym.pparfmem002.CpG_report.txt 1 0 
Sym.forSym.pparfmem003.CpG_report.txt 1 0 
Sym.forSym.pparimem010.CpG_report.txt 1 1 
Sym.forSym.pparimem011.CpG_report.txt 1 1 
Sym.forSym.pparimem012.CpG_report.txt 1 1

# ---------------
# design-matrix_p-vs-y.txt
base case
Sym.forSym.pparpmem004.CpG_report.txt 1 0 
Sym.forSym.pparpmem005.CpG_report.txt 1 0 
Sym.forSym.pparpmem006.CpG_report.txt 1 0
Sym.forSym.pparymem007.CpG_report.txt 1 1 
Sym.forSym.pparymem008.CpG_report.txt 1 1 
Sym.forSym.pparymem009.CpG_report.txt 1 1

# design-matrix_p-vs-i.txt
base case 
Sym.forSym.pparpmem004.CpG_report.txt 1 0
Sym.forSym.pparpmem005.CpG_report.txt 1 0
Sym.forSym.pparpmem006.CpG_report.txt 1 0
Sym.forSym.pparimem010.CpG_report.txt 1 1 
Sym.forSym.pparimem011.CpG_report.txt 1 1 
Sym.forSym.pparimem012.CpG_report.txt 1 1

# design-matrix_y-vs-i.txt
base case 
Sym.forSym.pparymem007.CpG_report.txt 1 0 
Sym.forSym.pparymem008.CpG_report.txt 1 0 
Sym.forSym.pparymem009.CpG_report.txt 1 0
Sym.forSym.pparimem010.CpG_report.txt 1 1 
Sym.forSym.pparimem011.CpG_report.txt 1 1 
Sym.forSym.pparimem012.CpG_report.txt 1 1

```

# Run DNMTools using `parallel_dnmtools.sh`:
> [!note]- parallel_dnmtools.sh
> ![[parallel_dnmtools.sh]]

