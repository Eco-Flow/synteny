use strict;
use warnings;

# Check for correct number of arguments
if (@ARGV != 1) {
    die "Usage: perl script.pl <input_file>\n";
}

# Get the input file name from command-line arguments
my $input_file = $ARGV[0];

# Extract the genus_species from the input file name
my ($genus_species) = $input_file =~ /^([^\.]+)/;

# Open the input file
open(my $in, '<', $input_file) or die "Cannot open $input_file: $!";

# Skip the header line
my $header = <$in>;

# Process the input file
my %output_files;
while (<$in>) {
    chomp;
    my @fields = split /\t/, $_;

    # Extract relevant columns
    my $final_classification = $fields[28]; # 29th column, 0-based index
    my $sp1_gene_b4          = $fields[4];  # 5th column
    my $sp1_gene_after       = $fields[8];  # 9th column

    # Open the output file for the classification if not already opened
    if (!exists $output_files{$final_classification}) {
        my $output_file_name = "${genus_species}.$final_classification.txt";
        open($output_files{$final_classification}, '>', $output_file_name)
          or die "Cannot open $output_file_name: $!";
    }

    # Write each gene name on a new line in the corresponding file
    my $fh = $output_files{$final_classification};
    print $fh "$sp1_gene_b4\n";
    print $fh "$sp1_gene_after\n";
}

# Close all output files
close($_) for values %output_files;

# Close the input file
close($in);

print "Processing complete. Output files created for each classification.\n";
