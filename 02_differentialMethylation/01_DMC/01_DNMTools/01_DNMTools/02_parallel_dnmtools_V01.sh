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

echo "All tasks completed."
