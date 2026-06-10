"""
DATE: 20230717
AUTHOR: Maximiliano Zuluaga-Forero
contact: mz448@cornell.edu
Cornell University, Sandkam Lab
Version: V02 (from process_for_meth-file_V1)

AIM: 
This code process a given list of files in the 'CpG_report' format given by -cytosisne_report in Bismark and produces formated '.meth' files compatible with DNMTools
"""
import argparse
import os

# Function to change file extension
def change_extension(filename, new_extension):
    base = os.path.splitext(filename)[0]
    return base + new_extension

# Function to calculate column 5
def calculate_column_5(column4, column5):
    if column4 + column5 == 0:
        return 0
    return column4 / (column4 + column5)

# Function to process each file
def process_file(original_file, output_file):
    if os.path.exists(output_file):
        print(f"Output file '{output_file}' already exists. Skipping...")
        return
    
    try:
        with open(original_file, 'r', encoding='latin-1') as input_file, open(output_file, 'w') as output:
            for line in input_file:
                columns = line.strip().split('\t')
                column1 = columns[0]
                column2 = columns[1]
                column3 = columns[2]
                column4 = columns[6]
                column5 = calculate_column_5(float(columns[3]), float(columns[4]))
                column6 = float(columns[3]) + float(columns[4])
                
                # Write the new line to the output file
                output.write(f"{column1}\t{column2}\t{column3}\t{column4}\t{column5}\t{column6}\n")
        
        print(f"Processed file '{original_file}'. Output saved as '{output_file}'.\n [DONE!]")
    except (FileNotFoundError, IsADirectoryError) as e:
        print(f"Error processing file '{original_file}': {str(e)}")

# Main function to parse arguments and call the processing function
def main():
    parser = argparse.ArgumentParser(description="Process CpG report files into .meth format.")
    parser.add_argument('input_file', help="Path to the input file.")
    parser.add_argument('output_file', help="Path to the output file.")
    
    args = parser.parse_args()
    
    process_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
