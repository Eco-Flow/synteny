import pandas as pd
import sys
import os

# python code written by Chris Wyatt to extract the junction locations (chromosome:start,end) from a file 

# Ensure proper usage
if len(sys.argv) != 4:
    print("Usage: python script.py <bed_file> <main_file> <output_file>")
    sys.exit(1)

# Assign command-line arguments to variables
bed_file = sys.argv[1]
main_file = sys.argv[2]
output_file = sys.argv[3]

# Extract species names from the main file name
filename = os.path.basename(main_file)
try:
    species_names = filename.split('.')[0:2]  # Assumes first two components are species names
    species1, species2 = species_names
except ValueError:
    print("Error: Could not parse species names from the file name.")
    sys.exit(1)

# Read BED file and create an ordered list of genes with coordinates
bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene', 'score', 'strand'])

# Create a mapping of each gene to its preceding gene and coordinates
gene_predecessor = {}
gene_coordinates = {}

for i in range(1, len(bed_df)):
    gene = bed_df['gene'].iloc[i]
    predecessor_gene = bed_df['gene'].iloc[i - 1]
    gene_predecessor[gene] = predecessor_gene
    gene_coordinates[gene] = (bed_df['chromosome'].iloc[i], bed_df['start'].iloc[i], bed_df['end'].iloc[i])
    gene_coordinates[predecessor_gene] = (bed_df['chromosome'].iloc[i - 1], bed_df['start'].iloc[i - 1], bed_df['end'].iloc[i - 1])

# Process the main file
with open(main_file, 'r') as main_f, open(output_file, 'w') as out_f:
    header = main_f.readline()  # Skip header
    # Write header for the output file
    out_f.write("event_type\tcurrent_gene\tpredecessor_gene\t"
                "current_chromosome\tcurrent_start\tcurrent_end\t"
                "predecessor_chromosome\tpredecessor_start\tpredecessor_end\t"
                "species1\tspecies2\n")
    
    for line in main_f:
        columns = line.strip().split('\t')
        event_type = columns[10]      # Column 8: type (e.g., INVER, TRANS, OTHER)
        current_gene = columns[9]    # Column 10: current gene name

        # Check if the gene exists in the BED file list and get the predecessor
        if current_gene in gene_predecessor:
            predecessor_gene = gene_predecessor[current_gene]
            # Retrieve coordinates for both genes
            current_gene_coords = gene_coordinates[current_gene]
            predecessor_gene_coords = gene_coordinates[predecessor_gene]
            
            # Write to the output file
            out_f.write(f"{event_type}\t{current_gene}\t{predecessor_gene}\t"
                        f"{current_gene_coords[0]}\t{current_gene_coords[1]}\t{current_gene_coords[2]}\t"
                        f"{predecessor_gene_coords[0]}\t{predecessor_gene_coords[1]}\t{predecessor_gene_coords[2]}\t"
                        f"{species1}\t{species2}\n")
        else:
            # Handle cases where there is no predecessor
            out_f.write(f"{event_type}\t{current_gene}\tNo_Predecessor_Found\t"
                        "N/A\tN/A\tN/A\tN/A\tN/A\tN/A\t"
                        f"{species1}\t{species2}\n")

print(f"Output saved to {output_file}")
