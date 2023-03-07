#!/usr/bin/perl
#Bed_to_histogram.pl
use warnings;
use strict;

die "Needs the input go file with extra title info to be removed\n" if (@ARGV!=1); 

my $gff = $ARGV[0];
my $outfile= "$gff\_simple.txt";

open(my $IN, "<", $gff)   or die "Could not open $gff \n";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

while (my $line=<$IN>){
    chomp $line;
    my @split=split("\t", $line);
    my @spacesplit=split("\ ", $split[0]);

    print $outhandle "$spacesplit[0]\t$split[1]\n";
}
