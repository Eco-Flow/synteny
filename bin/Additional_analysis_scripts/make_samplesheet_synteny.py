#!/usr/bin/env python3
import os
import csv

# Get the current working directory
cwd = os.getcwd()

# Get a list of all files in the current directory
files = os.listdir(cwd)

# Initialize dictionaries to store paths for each species
fa_files = {}
gff3_files = {}

# Separate .fa.gz and .gff3.gz files into different dictionaries based on species name
for file in files:
    if file.endswith('.fa.gz'):
        species_name = file.split('-')[0]
        fa_files[species_name] = os.path.join(cwd, file)
    elif file.endswith('.gff3.gz'):
        species_name = file.split('-')[0]
        gff3_files[species_name] = os.path.join(cwd, file)

# Open the output CSV file
with open('input.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Write the header (optional)
    # csvwriter.writerow(['Species', 'FA_File', 'GFF3_File'])
    
    # Find pairs and write to CSV
    for species_name in fa_files:
        if species_name in gff3_files:
            csvwriter.writerow([species_name, fa_files[species_name], gff3_files[species_name]])

print("input.csv file has been generated.")
