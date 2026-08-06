#!/bin/bash
#SBATCH --job-name=dnmtools_parallel
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=60
#SBATCH --mem=800G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# AUTHOR: MZF
# DATE: 20260718
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTES: 
#   Hardcoded path to preset the run via job manager
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
set -euo pipefail

# Add DNMtools to the path
export LD_LIBRARY_PATH=/programs/htslib-1.16/lib 
export PATH=/programs/dnmtools-1.2.2/bin:$PATH 
export PATH=/programs/parallel/bin:$PATH

WORKDIR="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/data"
cd "$WORKDIR" || exit 1

# 1) Radmeth
echo "Starting Radmeth..."
radmeth_commands=(
    "dnmtools radmeth -factor case design-matrix_f-vs-p.txt proportion-table_f-vs-p.txt >radmeth_f-vs-p.bed"
    "dnmtools radmeth -factor case design-matrix_f-vs-y.txt proportion-table_f-vs-y.txt >radmeth_f-vs-y.bed"
    "dnmtools radmeth -factor case design-matrix_f-vs-i.txt proportion-table_f-vs-i.txt >radmeth_f-vs-i.bed"
    "dnmtools radmeth -factor case design-matrix_p-vs-i.txt proportion-table_p-vs-i.txt >radmeth_p-vs-i.bed"
    "dnmtools radmeth -factor case design-matrix_p-vs-y.txt proportion-table_p-vs-y.txt >radmeth_p-vs-y.bed"
    "dnmtools radmeth -factor case design-matrix_y-vs-i.txt proportion-table_y-vs-i.txt >radmeth_y-vs-i.bed"
)

printf "%s\n" "${radmeth_commands[@]}" | parallel --jobs 6

# 2) Radadjust
echo "Starting Radadjust..."
radadjust_commands=(
    "dnmtools radadjust -bins 1:200:1 -v radmeth_f-vs-p.bed > radadjust_f-vs-p.bed"
    "dnmtools radadjust -bins 1:200:1 -v radmeth_f-vs-y.bed > radadjust_f-vs-y.bed"
    "dnmtools radadjust -bins 1:200:1 -v radmeth_f-vs-i.bed > radadjust_f-vs-i.bed"
    "dnmtools radadjust -bins 1:200:1 -v radmeth_p-vs-i.bed > radadjust_p-vs-i.bed"
    "dnmtools radadjust -bins 1:200:1 -v radmeth_p-vs-y.bed > radadjust_p-vs-y.bed"
    "dnmtools radadjust -bins 1:200:1 -v radmeth_y-vs-i.bed > radadjust_y-vs-i.bed"
)

printf "%s\n" "${radadjust_commands[@]}" | parallel --jobs 6

# 3) Filtering for significant positions only
echo "Filtering significant positions..."
filter_commands=(
    "awk '\$7 <= 0.01' radadjust_f-vs-p.bed > significant.radadjust_f-vs-p.bed"
    "awk '\$7 <= 0.01' radadjust_f-vs-y.bed > significant.radadjust_f-vs-y.bed"
    "awk '\$7 <= 0.01' radadjust_f-vs-i.bed > significant.radadjust_f-vs-i.bed"
    "awk '\$7 <= 0.01' radadjust_p-vs-i.bed > significant.radadjust_p-vs-i.bed"
    "awk '\$7 <= 0.01' radadjust_p-vs-y.bed > significant.radadjust_p-vs-y.bed"
    "awk '\$7 <= 0.01' radadjust_y-vs-i.bed > significant.radadjust_y-vs-i.bed"
)

printf "%s\n" "${filter_commands[@]}" | parallel --jobs 6


# 4) Calculate single-CpG posteiror methylation probability
echo "Calculating Posterior methylation probabilities per Phenotype..."
filter_commands=(
    "dnmtools hmr-rep -post-meth postProb_meth.f.tsv -p params.f.txt -s 123 Sym.forSym.pparfmem001.CpG_report.txt.meth Sym.forSym.pparfmem002.CpG_report.txt.meth Sym.forSym.pparfmem003.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.y.tsv -p params.y.txt -s 123 Sym.forSym.pparymem007.CpG_report.txt.meth Sym.forSym.pparymem008.CpG_report.txt.meth Sym.forSym.pparymem009.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.p.tsv -p params.p.txt -s 123 Sym.forSym.pparpmem004.CpG_report.txt.meth Sym.forSym.pparpmem005.CpG_report.txt.meth Sym.forSym.pparpmem006.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.i.tsv -p params.i.txt -s 123 Sym.forSym.pparimem010.CpG_report.txt.meth Sym.forSym.pparimem011.CpG_report.txt.meth Sym.forSym.pparimem012.CpG_report.txt.meth"
)
printf "%s\n" "${filter_commands[@]}" | parallel --jobs 8

echo "Calculating Posterior methylation probabilities Individually..."
filter_commands=(
    "dnmtools hmr-rep -post-meth postProb_meth.fmem001.tsv -p params.fmem001.txt -s 123 Sym.forSym.pparfmem001.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.fmem002.tsv -p params.fmem002.txt -s 123 Sym.forSym.pparfmem002.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.fmem003.tsv -p params.fmem003.txt -s 123 Sym.forSym.pparfmem003.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.pmem004.tsv -p params.pmem004.txt -s 123 Sym.forSym.pparpmem004.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.pmem005.tsv -p params.pmem005.txt -s 123 Sym.forSym.pparpmem005.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.pmem006.tsv -p params.pmem006.txt -s 123 Sym.forSym.pparpmem006.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.ymem007.tsv -p params.ymem007.txt -s 123 Sym.forSym.pparymem007.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.ymem008.tsv -p params.ymem008.txt -s 123 Sym.forSym.pparymem008.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.ymem009.tsv -p params.ymem009.txt -s 123 Sym.forSym.pparymem009.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.imem010.tsv -p params.imem010.txt -s 123 Sym.forSym.pparimem010.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.imem011.tsv -p params.imem011.txt -s 123 Sym.forSym.pparimem011.CpG_report.txt.meth"
    "dnmtools hmr-rep -post-meth postProb_meth.imem012.tsv -p params.imem012.txt -s 123 Sym.forSym.pparimem012.CpG_report.txt.meth"
)
printf "%s\n" "${filter_commands[@]}" | parallel --jobs 12

echo "All tasks completed."
