# DATE: 20260803
# AUTHOR: MZF
# Script: 04_compilePermanova_V01.sh
# Version 01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Compile PERMANOVA and PERMDISPERSION results into a single file
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
set -euo pipefail

: > PERM.analysis.compiled.tsv

for FILE in */PERM*; do
	echo "$FILE" >> PERM.analysis.compiled.tsv
    colcomma "$FILE" | heads >> PERM.analysis.compiled.tsv
	echo $'\n' >> PERM.analysis.compiled.tsv
done

