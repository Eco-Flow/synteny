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
my $input_file = $ARGV[1];

# Read the input file and store the values
open my $fh, '<', $input_file or die "Could not open file '$input_file': $!";

# Find out species we are using:
my @name=split(/\./, $input_file);
my $species=$name[0];
chomp $species;

#Save output of genes tested:
my $back="Background.txt";
open my $backsave, '>', $back or die "Could not open file '$back': $!";


my @data;
while (my $line = <$fh>) {
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
close $fh;
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

    print "ChopGO_VTS2_v12.pl -i $top_output_file --GO_file $go_key{$species} -bg Background.txt \n";
    print "ChopGO_VTS2_v12.pl -i $bottom_output_file --GO_file $go_key{$species} -bg Background.txt \n";
    `ChopGO_VTS2_v12.pl -i $top_output_file --GO_file $go_key{$species} -bg Background.txt`; 
    `ChopGO_VTS2_v12.pl -i $bottom_output_file --GO_file $go_key{$species} -bg Background.txt`;
}

