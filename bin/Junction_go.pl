#!/usr/bin/perl
use strict;
use warnings;

#For each species in your Go folder work out species name, then store in a hash with species name -> path to file
my @gos=`ls Go/*`;
my %go_key;
foreach my $sp (@gos){
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
}

# Ensure the correct number of arguments are provided
if (@ARGV != 2) {
    die "Usage: $0 <percentage> <input_file>\n";
}

my $percentage = $ARGV[0];
my $input_file = $ARGV[1];

# Read the input file and store the values
open my $fh, '<', $input_file or die "Could not open file '$input_file': $!";

# Find out species we are using:
my @name=split("\.", $input_file);
my $species=$name[0];


my @data;
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split /\t/, $line;
    # Skip if the value is 'NA'
    next if $fields[2] eq 'NA';
    push @data, [$fields[1], $fields[2]];  # Store gene and value
}
close $fh;

# Sort the data by the third column (numeric comparison)
@data = sort { $a->[1] <=> $b->[1] } @data;

# Calculate the number of elements in the top/bottom X%
my $num_elements = int(($percentage / 100) * scalar(@data));

# Get the cutoff values for top and bottom X%
my $bottom_cutoff = $data[$num_elements - 1][1];
my $top_cutoff = $data[-$num_elements][1];

# Prepare output filenames
my $top_output_file = "top_${percentage}_percent.txt";
my $bottom_output_file = "bottom_${percentage}_percent.txt";

# Open files for writing
open my $top_fh, '>', $top_output_file or die "Could not open file '$top_output_file': $!";
open my $bottom_fh, '>', $bottom_output_file or die "Could not open file '$bottom_output_file': $!";

# Filter and write genes to the output files
foreach my $entry (@data) {
    my ($gene, $value) = @$entry;
    if ($value >= $top_cutoff) {
        print $top_fh "$gene\n";
    }
    if ($value <= $bottom_cutoff) {
        print $bottom_fh "$gene\n";
    }
}

# Close the file handles
close $top_fh;
close $bottom_fh;

print "Results written to $top_output_file and $bottom_output_file\n";

#Run GO enrichment analysis on lists
`ChopGO_VTS2_v12.pl -i $top_output_file --GO_file $go_key{$species} -bg $species\.$percentage\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $bottom_output_file --GO_file $go_key{$species} -bg $species\.$percentage\.bg.txt`;