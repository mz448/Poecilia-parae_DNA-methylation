```bash
nano histogram_normalizer_06.py
chmod +x histogram_normalizer_06.py
```

```python
#!/usr/bin/env python3
"""
DATE: 20250103
AUTHOR: MAXIMILIANO ZULUAGA FORERO
VERSION: 06
DESCRIPTION:
This script compares two histogram files (tab-delimited `.txt` files) to generate a normalized histogram. 
Each input file must have four columns:
1. Chromosome name
2. Start position of the bin
3. End position of the bin
4. Count of sites in the bin

The script calculates the ratio of counts in the "counts file" to those in the "reference file" for rows where the first three columns match. 
The output is saved as a new `.txt` file with the same structure, but the fourth column contains the normalized counts.

USAGE:
Run this script with three arguments:
1. Reference histogram file
2. Counts histogram file
3. Output file name

Example for bash for-loop:
    for FILE in Histogram*; do
        ./normalize_histograms_ref_counts.py reference_file.txt $FILE normalized_$FILE.txt
    done
"""

import sys

def load_histogram_data(file_path):
    """Reads a histogram file and returns data as a dictionary."""
    histogram_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            chromosome, start_position, end_position, bin_count = (
                columns[0],
                int(columns[1]),
                int(columns[2]),
                int(columns[3]),
            )
            histogram_data[(chromosome, start_position, end_position)] = bin_count
    return histogram_data

def main():
    if len(sys.argv) != 4:
        print("Usage: normalize_histograms_ref_counts.py <reference_file> <counts_file> <output_file>")
        sys.exit(1)

    reference_file_path = sys.argv[1]
    counts_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    print(f"Starting normalization: {counts_file_path} (counts) normalized to {reference_file_path} (reference)")

    reference_histogram = load_histogram_data(reference_file_path)
    counts_histogram = load_histogram_data(counts_file_path)

    with open(output_file_path, 'w') as output_file:
        for bin_key in counts_histogram:
            if bin_key in reference_histogram:
                counts_value = counts_histogram[bin_key]
                reference_value = reference_histogram[bin_key]

                # Calculate the ratio, avoiding division by zero
                normalized_ratio = 0 if reference_value == 0 else counts_value / reference_value

                chromosome, start_position, end_position = bin_key
                output_file.write(f"{chromosome}\t{start_position}\t{end_position}\t{normalized_ratio:.6f}\n")

    print(f"Normalization complete. Results saved to {output_file_path}")

if __name__ == "__main__":
    main()

```