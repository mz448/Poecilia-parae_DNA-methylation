#!/usr/bin/env python3

"""
DATE: 20250109
AUTHOR: MAXIMILIANO ZULUAGA FORERO
VERSION: 02

DESCRIPTION:
This script calculates the total count of genomic sites per chromosome using a karyotype file and a sites file. 
The output is a text file containing chromosome names and site counts.

Usage:
    python counts-per-chromosome.py <karyotype_file> <sites_file> <output_file>

Parameters:
    <karyotype_file> : Path to the karyotype file (tab-separated with chromosome names and sizes).
    <sites_file>     : Path to the sites file (tab-separated with chromosome names and positions).
    <output_file>    : Path to the output file for chromosome counts.

Version 02: 
The script calculates total site counts per chromosome instead of histogram bins. This script is based in the counts-for-histogram_02.py script
V03:
  Works with files with more than 2 coliumns. Since only counts rows. 
"""
def verbose_message(message):
    """
    Prints a verbose message to indicate script progress.
    """
    print(f"[INFO] {message}")

def main():
    import sys

    if len(sys.argv) != 4:
        print("Usage: python counts-per-chromosome.py <karyotype_file> <sites_file> <output_file>")
        sys.exit(1)

    chromosome_file_path = sys.argv[1]
    sites_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    verbose_message(f"Starting processing: {sites_file_path}")

    # Read karyotype data to get chromosome names
    verbose_message("Reading karyotype file...")
    chromosomes = set()
    with open(chromosome_file_path, 'r') as chromosome_file:
        for line in chromosome_file:
            chromosome, _ = line.strip().split('\t')
            chromosomes.add(chromosome)
    verbose_message("Karyotype file processing complete.")

    # Initialize site counts per chromosome
    site_counts = {chromosome: 0 for chromosome in chromosomes}

    # Read sites data and count sites per chromosome
    verbose_message("Reading and processing sites file...")
    # with open(sites_file_path, 'r') as sites_file:
    #     for line in sites_file:
    #         chromosome, _ = line.strip().split('\t')
    #         if chromosome in site_counts:
    #             site_counts[chromosome] += 1
    
    with open(sites_file_path, 'r') as sites_file:
      for line in sites_file:
        chromosome = line.rstrip('\n').split('\t')[0]
        if chromosome in site_counts:
          site_counts[chromosome] += 1
    
    verbose_message("Sites file processing complete.")

    # Sort site counts by chromosome name
    sorted_site_counts = sorted(site_counts.items())

    # Write output to file
    verbose_message(f"Writing results to output file: {output_file_path}")
    with open(output_file_path, 'w') as output_file:
        for chromosome, count in sorted_site_counts:
            output_file.write(f"{chromosome}\t{count}\n")
    verbose_message("Processing complete.")

if __name__ == "__main__":
    main()
