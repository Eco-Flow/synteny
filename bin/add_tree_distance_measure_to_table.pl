#!/usr/bin/perl
use strict;
use warnings;

# Check if the required arguments are provided
if (@ARGV != 2) {
    die "Usage: $0 <pairs_file> <matrix_file>\n";
}

# Get the file names from command line arguments
my ($pairs_file, $matrix_file) = @ARGV;

# Open the matrix file
open(my $matrix_fh, '<', $matrix_file) or die "Could not open '$matrix_file' $!";

# Read the header of the matrix file
my $header = <$matrix_fh>;
chomp $header;
my @matrix_species = split(/\t/, $header);

# Store the matrix values in a hash
my %matrix;
while (my $line = <$matrix_fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    my $species1 = shift @fields;  # Take out the first species
    for (my $i = 0; $i < @fields; $i++) {
        $matrix{"$species1.$matrix_species[$i]"} = $fields[$i];
        $matrix{"$matrix_species[$i].$species1"} = $fields[$i];  # Account for reverse order
    }
}
close $matrix_fh;

# Open the pairs file
open(my $pairs_fh, '<', $pairs_file) or die "Could not open '$pairs_file' $!";

# Read and process the pairs file
my $pairs_header = <$pairs_fh>;
chomp $pairs_header;
print "$pairs_header\tMatrix_Value\n";  # Add a new column in the output

while (my $line = <$pairs_fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    my $pair = $fields[0];  # Assuming species pair is in the first column
    my $value = exists $matrix{$pair} ? $matrix{$pair} : 'NA';  # Get value or 'NA' if not found
    print "$line\t$value\n";
}

close $pairs_fh;