#!/usr/bin/perl
use strict;
use warnings;

#For each species in your Go folder work out species name, then store in a hash with species name -> path to file
my @gos=`ls Go/*`;
my %go_key;
foreach my $sp (@gos){
    chomp $sp;
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
}

# Ensure the correct number of arguments are provided
if (@ARGV != 2) {
    die "Usage: $0 <percentage> <input_file>\n";
}

my $percentage = $ARGV[0];
my $species = $ARGV[1];

# Store the total scores and counts for each key
my %scores;
my %counts;

# Read all input files ending with '_gene_scores.txt'
my @files = glob("*.gene_scores.txt");

foreach my $file (@files) {
    open my $fh, '<', $file or die "Could not open file '$file' $!";
    print "Combining $file\n";
    while (<$fh>) {
        chomp;
        my ($col1, $col2, $score) = split /\t/;

        # Create a unique key for each row based on the first two columns
        my $key = "$col1\t$col2";

        # Accumulate scores and increment count for each key
        $scores{$key}+= $score;
        $counts{$key}++;
    }

    close $fh;
}

# Write the averaged scores to "Summary.tsv"
open my $out_fh, '>', "Summary.tsv" or die "Could not open file 'Summary.tsv' $!";

foreach my $key (keys %scores) {
    # Calculate average score
    my $average = $scores{$key} / $counts{$key};

    # Write to output file
    print $out_fh "$key\t$average\n";
}

close $out_fh;

print "Averaged scores written to 'Summary.tsv'.\n";

# Find out species we are using:
chomp $species;

#Save output of genes tested:
my $back="$species\.$percentage\.bg.txt";
open my $backsave, '>', $back or die "Could not open file '$back': $!";


#Read in the new merged table;
my $sumtab="Summary.tsv";
open my $fh2, '<', $sumtab or die "Could not open file '$sumtab' $!";

my @data;
while (my $line = <$fh2>) {
    chomp $line;
    my @fields = split /\t/, $line;
    # Skip if the value is 'NA'
    next if $fields[2] eq 'NA';
    
    if($fields[1] =~ m/\:/){
                my @sp1=split(/\:/, $fields[1]);
                $fields[1]=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($fields[1] =~ m/rna-/){
                my @sp1=split(/\-/, $fields[1]);
                $fields[1]=$sp1[1];
            }
    if($fields[1] =~ m/-/){
        $fields[1]=~ s/\-/\_/g;
    }
    push @data, [$fields[1], $fields[2]];  # Store gene and value
    print $backsave "$fields[1]\n";
}
close $fh2;
close $backsave;

# Sort the data by the third column (numeric comparison)
@data = sort { $a->[1] <=> $b->[1] } @data;

# Calculate the number of elements in the top/bottom X%
my $num_elements = int(($percentage / 100) * scalar(@data));

print "$num_elements\n";

#Check whether to plot GO or not (we need some values in the score files)
if ($num_elements == 0){
    print "There were no elements to sort into top and bottom genes. Likely this is because there were no junctions to check the distance to them\n";
}
else{
    # Get the cutoff values for top and bottom X%
    my $bottom_cutoff = $data[$num_elements - 1][1];
    my $top_cutoff = $data[-$num_elements][1];

    # Prepare output filenames
    my $top_output_file = "$species.${percentage}\.top.txt";
    my $bottom_output_file = "$species.${percentage}\.bottom.txt";

    # Open files for writing
    open my $top_fh, '>', $top_output_file or die "Could not open file '$top_output_file': $!";
    open my $bottom_fh, '>', $bottom_output_file or die "Could not open file '$bottom_output_file': $!";

    # Filter and write genes to the output files
    foreach my $entry (@data) {
        my ($gene, $value) = @$entry;
        if ($value >= $top_cutoff) {
            #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
            if($gene =~ m/\:/){
                my @sp1=split(/\:/, $gene);
                $gene=$sp1[1];
            }
            # Rename weird NCBI id rna- prefix
            if($gene =~ m/rna-/){
                my @sp1=split(/\-/, $gene);
                $gene=$sp1[1];
            }
            if($gene =~ m/-/){
                $gene=~ s/\-/\_/g;
            } 
            print $top_fh "$gene\n";
        }
        if ($value <= $bottom_cutoff) {
            #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
            if($gene =~ m/\:/){
                my @sp1=split(/\:/, $gene);
                $gene=$sp1[1];
            }
            # Rename weird NCBI id rna- prefix
            if($gene =~ m/rna-/){
                my @sp1=split(/\-/, $gene);
                $gene=$sp1[1];
            }
            if($gene =~ m/-/){
                $gene=~ s/\-/\_/g;
            }
            print $bottom_fh "$gene\n";
        }
    }

    # Close the file handles
    close $top_fh;
    close $bottom_fh;

    print "Results written to $top_output_file and $bottom_output_file\n";

    #Run GO enrichment analysis on lists

    print "ChopGO_VTS2_v12.pl -i $top_output_file --GO_file $go_key{$species} -bg $species\.$percentage\.bg.txt \n";
    print "ChopGO_VTS2_v12.pl -i $bottom_output_file --GO_file $go_key{$species} -bg $species\.$percentage\.bg.txt \n";
    `ChopGO_VTS2_v12.pl -i $top_output_file --GO_file $go_key{$species} -bg $species\.$percentage\.bg.txt`; 
    `ChopGO_VTS2_v12.pl -i $bottom_output_file --GO_file $go_key{$species} -bg $species\.$percentage\.bg.txt`;
}

