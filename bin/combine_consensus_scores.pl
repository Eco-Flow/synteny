#!/usr/bin/perl
use warnings;
use strict;

my $sp = "Summary_of_pairwise_comparisons*";
my @files = glob($sp);

if (!@files) {
    die "No files matching the pattern $sp found\n";
}

# Open the first file that matches
open(my $INPUT, "<", $files[0]) or die "Could not open $files[0]: $!\n";

my $outfile = "filec";
open(my $out, ">", $outfile) or die "Could not open $outfile\n";

#Remove header:
my $head=<$INPUT>;
chomp $head;
print $out "$head\tComparison\tTranslocation_junctions\tInversion_junctions\tSame_direction_duplication_junctions\tLoop_direction_duplication_junctions\n";


#Read in the input table
while (my $line=<$INPUT>){
    chomp $line;
    my @sp=split("\t", $line);
    my $comp=$sp[0];
    print $out "$line";

    my $match = "$comp\.lifted.anchors.classified";
    open(my $INPUT_MATCH, "<", $match) or die "Could not open $match\n";
    my $head=<$INPUT_MATCH>;
    my $data=<$INPUT_MATCH>;
    chomp $data;
    print $out "\t$data\n";
}


