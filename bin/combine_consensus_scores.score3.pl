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
print $out "$head\tJunc_inter\tJunc_inver\tJunc_indel\tJunc_indel_lt_5\tJunc_indel_5_20\tJunc_indel_gt_20\n";


#Read in the input table
while (my $line=<$INPUT>){
    chomp $line;
    my @sp=split("\t", $line);
    my $comp=$sp[0];
    print $out "$line";

    my $match = "$comp\_junction_summary.tsv";
    open(my $INPUT_MATCH, "<", $match) or die "Could not open $match\n";
    my $head=<$INPUT_MATCH>;

    my $inter=<$INPUT_MATCH>;
    chomp $inter;
    my @inter_s=split("\t",$inter);

    my $inver=<$INPUT_MATCH>;
    chomp $inver;
    my @inver_s=split("\t",$inver);

    my $indel=<$INPUT_MATCH>;
    chomp $indel;
    my @indel_s=split("\t",$indel);

    my $indel_0_5=<$INPUT_MATCH>;
    chomp $indel_0_5;
    my @indel_0_5_s=split("\t",$indel_0_5);

    my $indel_5_20=<$INPUT_MATCH>;
    chomp $indel_5_20;
    my @indel_5_20_s=split("\t",$indel_5_20);

    my $indel_20_=<$INPUT_MATCH>;
    chomp $indel_20_;
    my @indel_20_s=split("\t",$indel_20_);

    print $out "\t$inter_s[1]\t$inver_s[1]\t$indel_s[1]\t$indel_0_5_s[1]\t$indel_5_20_s[1]\t$indel_20_s[1]\n";
}


