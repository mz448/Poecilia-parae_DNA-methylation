#!/usr/bin/env python3

"""
DATE: 20250109
AUTHOR: MAXIMILIANO ZULUAGA FORERO
VERSION: 01

DESCRIPTION:
This script normalizes genomic site counts per chromosome by comparing a counts file to a reference file. 
The output is a text file containing chromosome names and normalized counts.

Usage:
    python normalize-chromosome-counts.py <reference_file> <counts_file> <output_file>

Parameters:
    <reference_file> : Path to the reference file (tab-separated with chromosome names and site counts).
    <counts_file>    : Path to the counts file (tab-separated with chromosome names and site counts).
    <output_file>    : Path to the output file for normalized counts.

Version 01: The script is an adapted version of histogram_normalizer_06.py to work with files only containing chromosome names and total site counts.
"""

def verbose_message(message):
    """
    Prints a verbose message to indicate script progress.
    """
    print(f"[INFO] {message}")

def read_counts_file(file_path):
    """
    Reads a file with chromosome names and counts, returning a dictionary.
    """
    counts = {}
    with open(file_path, 'r') as file:
        for line in file:
            chromosome, count = line.strip().split('\t')
            counts[chromosome] = int(count)
    return counts

def main():
    import sys

    if len(sys.argv) != 4:
        print("Usage: python normalize-chromosome-counts.py <reference_file> <counts_file> <output_file>")
        sys.exit(1)

    reference_file_path = sys.argv[1]
    counts_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    verbose_message(f"Starting normalization: {counts_file_path} (counts) normalized to {reference_file_path} (reference)")

    # Read reference and counts files
    verbose_message("Reading reference file...")
    reference_counts = read_counts_file(reference_file_path)
    verbose_message("Reference file read successfully.")

    verbose_message("Reading counts file...")
    counts = read_counts_file(counts_file_path)
    verbose_message("Counts file read successfully.")

    # Normalize counts
    verbose_message("Calculating normalized counts...")
    normalized_counts = {}
    for chromosome, count in counts.items():
        reference_count = reference_counts.get(chromosome, 0)
        normalized_counts[chromosome] = 0 if reference_count == 0 else count / reference_count

    # Sort normalized counts by chromosome name
    sorted_normalized_counts = sorted(normalized_counts.items())

    # Write output to file
    verbose_message(f"Writing results to output file: {output_file_path}")
    with open(output_file_path, 'w') as output_file:
        for chromosome, normalized_count in sorted_normalized_counts:
            output_file.write(f"{chromosome}\t{normalized_count:.6f}\n")
    verbose_message("Processing complete.")

if __name__ == "__main__":
    main()
