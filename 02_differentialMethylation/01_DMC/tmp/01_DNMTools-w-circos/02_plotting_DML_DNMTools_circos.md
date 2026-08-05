____
**DATE:**  20250102 
**AUTHOR:** MZF
_____
**GOAL :**
Make Circos plots
____
# Path
```bash
cd /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_methylation/01_DNMTools/plotting_DML_circos
```
------
# 1) Make karyotype
```bash
cd Projects/ppar_emseq/data/009_emseq_pparae_muscle/genome/

# Find all scaffold names and lengths (and remove the header)
/programs/seqkit-0.15.0/seqkit fx2tab --length --name --header-line PparFemMitoLambdaPuc19Ver2024.fasta | tail -n +2 > scaff_lengths.txt

head -n 23 scaff_lengths.txt >head.txt &&
tail -n 3 scaff_lengths.txt > tail.txt &&
cat head.txt tail.txt > karyotype.PparFemMitoLambdaPuc19Ver2024.txt &&
rm scaff_lengths.txt head.txt tail.txt

# Print in proper format
# awk '{print "chr" "\t" "-" "\t" $1 "\t" $1 "\t" "0" "\t" $2 "\t" $1}' scaff_lengths.txt > karyotype.parae.txt

# copy into circos folder
cp /local/storage/Projects/ppar_emseq/data/009_emseq_pparae_muscle/genome/karyotype.PparFemMitoLambdaPuc19Ver2024.txt /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_methylation/01_DNMTools/plotting_DML_circos


```

# 2)  Make binned files (100kb bins)

#### NOTE: Reference file comes from radmeth: 
```bash
# radmeth structure
ch position strand CpG Pvalue reads-Control methylated-Control reads-Case methylated-Case
# Check radmeth lenght (is the same, that is why I use only one file as reference)
wc -l radmeth*
  15942610 radmeth_f-vs-i.bed
  15942610 radmeth_f-vs-p.bed
  15942610 radmeth_f-vs-y.bed
  15942610 radmeth_p-vs-i.bed
  15942610 radmeth_p-vs-y.bed
  15942610 radmeth_y-vs-i.bed
```

## 2.1) Provide simpler count files with chromosome and position only
```bash
# 1) Generate a reference file that will be used to normalize counts. This file contains all the possible sites that can be methylated (all CpGs)
awk '{print $1 "\t" $2}' ../radmeth_f-vs-i.bed > all.CpGpositions.txt &&
# 2) Generate count files  with only the first two columns representing the sites that are methilated in each comparison
awk '{print $1 "\t" $2}' ../significant.radadjust_f-vs-i.bed > CpGpositions_significant.radadjust_f-vs-i.txt &&
awk '{print $1 "\t" $2}' ../significant.radadjust_f-vs-p.bed > CpGpositions_significant.radadjust_f-vs-p.txt &&
awk '{print $1 "\t" $2}' ../significant.radadjust_f-vs-y.bed > CpGpositions_significant.radadjust_f-vs-y.txt &&
awk '{print $1 "\t" $2}' ../significant.radadjust_p-vs-i.bed > CpGpositions_significant.radadjust_p-vs-i.txt &&
awk '{print $1 "\t" $2}' ../significant.radadjust_p-vs-y.bed > CpGpositions_significant.radadjust_p-vs-y.txt &&
awk '{print $1 "\t" $2}' ../significant.radadjust_y-vs-i.bed > CpGpositions_significant.radadjust_y-vs-i.txt

```

## 2.2) Calculate the counts per bin using a custom python script
>[!note]- Python scripts
> ![[counts-for-histogram.py]]
> 
> ![[counts-for-histogram_02.py]]

```bash
# make sure it is executable
chmod +x counts-for-histogram.py
# Run counts-for-histogram.py in a for-loop
# for FILE in CpGpositions*; do ./counts-for-histogram.py karyotype* $FILE output.$FILE.txt done

# Run counts-for-histogram_02.py in a for-loop
for FILE in CpGpositions*; do 
./counts-for-histogram_02.py karyotype.PparFemMitoLambdaPuc19Ver2024.txt  $FILE output.$FILE 
done
# Make the reference histogram using the All.CpGpositions.txt file as the input
./counts-for-histogram_02.py karyotype.PparFemMitoLambdaPuc19Ver2024.txt  All.CpGpositions.txt reference.All.CpGpositions.txt 

```

# 3) Normalize the histograms, comparing each counts histogram against the reference

>[!note]- python scripts
> ![[histogram_normalizer_06.py]]

```bash
# run the normalization script using a for loop
for FILE in output*; do
./histogram_normalizer_06.py reference.All.CpGpositions.txt  $FILE normalized_$FILE.txt
done
```



# 4) Generate CIRCOS karyotype
```bash
# trim the karyotype to include chromosomes only
head -n 23 karyotype.PparFemMitoLambdaPuc19Ver2024.txt > karyotype.PparFemVer2024.txt

# Make the formatted karyotype using awk
awk '{
  # Define the chromosome prefix
  chr_name = $1;
  size = $2;
  
  # Extract numeric part of the chromosome name (after "Parae_")
  match(chr_name, /Parae_([0-9]+)/, arr);
  chr_number = arr[1];
  
  # Print the output in the required format
  start = 0;
  end = size;
  print "chr\t-\t" chr_name "\t" chr_number "\t" start "\t" end "\tchr" chr_number;
}' karyotype.PparFemVer2024.txt | sort -k6,6nr > karyotype_circos.PparFemVer2024.txt

```

# 5) Define Circos files
>[!note]- CircosFiles_morph
> ![[CircosFiles_morph]]

>[!note]- CircosFiles_sex
> ![[CircosFiles_sex]]

# 6) Run circos plots

```bash
#Starting the circos environment
export PATH=/local/storage/Environments/Anaconda3/bin:$PATH
source activate Circos
# Run the config file
circos -conf config.morph.conf
circos -conf config.sex.conf
```

# 6.1) Check the scale of plots
```bash
cd /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_methylation/01_DNMTools/plotting_DML_circos
# Identify each track's maximum value
for FILE in normalized*; do
ls $FILE && awk '{print $4}' $FILE | sort | tail -n 1
done
```
>[!summary]- Scales
>normalized_output.CpGpositions_significant.radadjust_f-vs-i.txt.txt
0.165349
normalized_output.CpGpositions_significant.radadjust_f-vs-p.txt.txt
0.168643
normalized_output.CpGpositions_significant.radadjust_f-vs-y.txt.txt
0.200264
normalized_output.CpGpositions_significant.radadjust_p-vs-i.txt.txt
0.075316
normalized_output.CpGpositions_significant.radadjust_p-vs-y.txt.txt
0.082359
normalized_output.CpGpositions_significant.radadjust_y-vs-i.txt.txt
0.056439

# 7) Do additional plotting using R