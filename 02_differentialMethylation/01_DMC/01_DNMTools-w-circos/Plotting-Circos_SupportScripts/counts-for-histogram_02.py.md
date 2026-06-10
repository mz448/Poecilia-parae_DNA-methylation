
```python 
#!/usr/bin/env python3

"""
DATE: 20250103
AUTHOR: MAXIMILIANO ZULUAGA FORERO
VERSION: 02

DESCRIPTION:
This script calculates histograms of 100kb bins for genomic data. It takes a karyotype file with chromosome names and sizes,
and a sites file with genomic positions, to compute the number of sites in each bin. The output is a text file containing
chromosome, start, end, and site count for each bin.

Usage:
    python counts-for-histogram.py <karyotype_file> <sites_file> <output_file>

Parameters:
    <karyotype_file> : Path to the karyotype file (tab-separated with chromosome names and sizes).
    <sites_file>     : Path to the sites file (tab-separated with chromosome names and positions).
    <output_file>    : Path to the output file for histogram data.

Version 02: The script now includes verbose messages for progress updates.
"""

def calculate_intervals(chromosome_length, interval_size=100000):
    """
    Calculates non-overlapping intervals (bins) for a given chromosome length.
    """
    intervals = []
    start = 1  # Intervals are 1-based
    while start <= chromosome_length:
        end = min(start + interval_size - 1, chromosome_length)
        intervals.append((start, end))
        start = end + 1
    return intervals

def count_sites_in_intervals(sites, intervals):
    """
    Counts the number of sites that fall within each interval.
    """
    interval_counts = [0] * len(intervals)
    for site in sites:
        for i, (start, end) in enumerate(intervals):
            if start <= site <= end:
                interval_counts[i] += 1
                break
    return interval_counts

def verbose_message(message):
    """
    Prints a verbose message to indicate script progress.
    """
    print(f"[INFO] {message}")

def main():
    import sys

    if len(sys.argv) != 4:
        print("Usage: python counts-for-histogram.py <karyotype_file> <sites_file> <output_file>")
        sys.exit(1)

    chromosome_file_path = sys.argv[1]
    sites_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    verbose_message(f"Starting processing: {sites_file_path}")

    # Read karyotype data and calculate intervals
    verbose_message("Reading karyotype file...")
    intervals_per_chromosome = {}
    with open(chromosome_file_path, 'r') as chromosome_file:
        for line in chromosome_file:
            chromosome, length = line.strip().split('\t')
            intervals_per_chromosome[chromosome] = calculate_intervals(int(length))
    verbose_message("Karyotype file processing complete.")

    # Read sites data
    verbose_message("Reading sites file...")
    sites = []
    with open(sites_file_path, 'r') as sites_file:
        for line in sites_file:
            chromosome, position = line.strip().split('\t')
            sites.append((chromosome, int(position)))
    verbose_message("Sites file processing complete.")

    # Process each chromosome and calculate counts
    verbose_message("Calculating histogram bins and counts...")
    output_lines = []
    for chromosome, intervals in intervals_per_chromosome.items():
        interval_counts = count_sites_in_intervals([site[1] for site in sites if site[0] == chromosome], intervals)
        for (start, end), count in zip(intervals, interval_counts):
            output_lines.append(f"{chromosome}\t{start}\t{end}\t{count}")
    verbose_message("Histogram calculation complete.")

    # Write output to file
    verbose_message(f"Writing results to output file: {output_file_path}")
    with open(output_file_path, 'w') as output_file:
        output_file.write('\n'.join(output_lines))
    verbose_message("Processing complete.")

if __name__ == "__main__":
    main()

```

```bash
# Provide executable privileges
chmod +x counts-for-histogram_02.py
```