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
    die "Usage: $0 <distance_in_genes> <input_file>\n";
}

my $distance = $ARGV[0];
my $species = $ARGV[1];

# Store the total scores and counts for each key
my %scores;

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

        # Add new score if less than prior:
        if ($scores{$key}){
            if ($scores{$key} eq "NA"){
                $scores{$key} = $score;
            }
            else{
                if ($scores{$key} >= $score){
                    $scores{$key} = $score;
                }
            }
            
        }
        else{
            $scores{$key} = $score;
        }
        
        
    }

    close $fh;
}

# Write the lowest scores to "Summary.tsv"
open my $out_fh, '>', "Summary.tsv" or die "Could not open file 'Summary.tsv' $!";

#Save output of genes tested:
my $back="$species\.$distance\.bg.txt";
open my $backsave, '>', $back or die "Could not open file '$back': $!";

my $closest_output_file = "$species.${distance}\.closest.txt";
open my $closest_fh, '>', $closest_output_file or die "Could not open file '$closest_output_file': $!";

# Count of distance scores from break, e.g. <1 scores or <3 scores
my $distance_count=0;
my @closest_cutoff;
foreach my $key (keys %scores) {
    # Write to output file
    print $out_fh "$key\t$scores{$key}\n";
    if ($scores{$key} eq "NA"){
        #Do nothing.
    }
    else{
        if ($scores{$key} <= $distance){
            $distance_count++;
            push (@closest_cutoff, $key);
            my @splh=split("\t", $key);
            if($splh[1] =~ m/\:/){
                my @sp1=split(/\:/, $splh[1]);
                $splh[1]=$sp1[1];
            }
            # Rename weird NCBI id rna- prefix
            if($splh[1] =~ m/rna-/){
                        my @sp1=split(/\-/, $splh[1]);
                        $splh[1]=$sp1[1];
                    }
            if($splh[1] =~ m/-/){
                $splh[1] =~ s/\-/\_/g;
            }
            print $closest_fh "$splh[1]\n";
        }
        my @splh=split("\t", $key);
        if($splh[1] =~ m/\:/){
            my @sp1=split(/\:/, $splh[1]);
            $splh[1]=$sp1[1];
        }
        # Rename weird NCBI id rna- prefix
        if($splh[1] =~ m/rna-/){
                    my @sp1=split(/\-/, $splh[1]);
                    $splh[1]=$sp1[1];
                }
        if($splh[1] =~ m/-/){
            $splh[1] =~ s/\-/\_/g;
        }
        print $backsave "$splh[1]\n";
    }
}

close $out_fh;
close $backsave;

print "Lowest scores written to 'Summary.tsv'. \n $distance_count genes were below the threshold of $distance\n\n";

# Find out species we are using:
chomp $species;


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
    
}
close $fh2;


# Sort the data by the third column (numeric comparison)
@data = sort { $a->[1] <=> $b->[1] } @data;

my @last_values = @data[-$distance_count..-1];
my $farthest_output_file = "$species.${distance}\.farthest.txt";
open my $farthest_fh, '>', $farthest_output_file or die "Could not open file '$farthest_output_file': $!";
foreach my $entry (@last_values) {
        my ($gene, $value) = @$entry;
        #print "$gene and $value\n";
        print $farthest_fh "$gene\n";
}


print "Results written to $farthest_output_file and $closest_output_file\n";

#Run GO enrichment analysis on lists

print "ChopGO_VTS2_v12.pl -i $farthest_output_file --GO_file $go_key{$species} -bg $species\.$distance\.bg.txt \n";
print "ChopGO_VTS2_v12.pl -i $closest_output_file --GO_file $go_key{$species} -bg $species\.$distance\.bg.txt \n";
`ChopGO_VTS2_v12.pl -i $farthest_output_file --GO_file $go_key{$species} -bg $species\.$distance\.bg.txt`; 
`ChopGO_VTS2_v12.pl -i $closest_output_file --GO_file $go_key{$species} -bg $species\.$distance\.bg.txt`;


