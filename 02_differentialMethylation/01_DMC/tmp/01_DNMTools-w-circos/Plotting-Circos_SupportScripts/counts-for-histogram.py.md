```python
#!/usr/bin/env python3
"""
DATE: 20250103
AUTHOR: MAXIMILIANO ZULUAGA FORERO
VERSION: 00

DESCRIPTION:
This script generates histograms of counts in fixed-size genomic bins (default: 100kb) based on a karyotype file 
and a file containing genomic sites of interest. The karyotype file provides chromosome names and sizes, while 
the sites file lists the positions of interest within those chromosomes. The script calculates the number of sites 
falling into each bin and outputs the results to a specified file.

USAGE:
The script accepts three arguments:
1. Karyotype file: A tab-separated file with two columns: chromosome name and chromosome length.
2. Sites file: A tab-separated file with two columns: chromosome name and site position.
3. Output file: The path to save the resulting histogram data.

The script is designed to be used in a bash loop, enabling the processing of multiple datasets 
by iteratively passing the input and output file names.

Example:
python script.py karyotype.txt sites.txt output.txt
"""

import sys

def calculate_intervals(chromosome_length, interval_size=100000):
    """Generates fixed-size intervals for a given chromosome length."""
    intervals = []
    start = 1  # Intervals are 1-based
    while start <= chromosome_length:
        end = min(start + interval_size - 1, chromosome_length)
        intervals.append((start, end))
        start = end + 1
    return intervals

def count_sites_in_intervals(sites, intervals):
    """Counts how many sites fall within each interval."""
    interval_counts = [0] * len(intervals)
    for site in sites:
        for i, (start, end) in enumerate(intervals):
            if start <= site <= end:
                interval_counts[i] += 1
                break
    return interval_counts

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <karyotype_file> <sites_file> <output_file>")
        sys.exit(1)

    chromosome_file_path = sys.argv[1]
    sites_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    # Parse karyotype file
    intervals_per_chromosome = {}
    with open(chromosome_file_path, 'r') as chromosome_file:
        for line in chromosome_file:
            chromosome, length = line.strip().split('\t')
            intervals_per_chromosome[chromosome] = calculate_intervals(int(length))

    # Parse sites file
    sites = []
    with open(sites_file_path, 'r') as sites_file:
        for line in sites_file:
            chromosome, position = line.strip().split('\t')
            sites.append((chromosome, int(position)))

    # Compute histogram
    output_lines = []
    for chromosome, intervals in intervals_per_chromosome.items():
        interval_counts = count_sites_in_intervals(
            [site[1] for site in sites if site[0] == chromosome], intervals
        )
        for (start, end), count in zip(intervals, interval_counts):
            output_lines.append(f"{chromosome}\t{start}\t{end}\t{count}")

    # Write output
    with open(output_file_path, 'w') as output_file:
        output_file.write('\n'.join(output_lines))

if __name__ == "__main__":
    main()


```

