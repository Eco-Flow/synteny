#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;


print "Please be in folder with all the anchor files\n";

my $outfile="My_scores.tsv";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
my @files=`ls *.anchors`;

print $outhandle "Sp1\tSp2\tSyn_lines\tMax_length\tAverage_length\n";

foreach my $file (@files){
    chomp $file;
	my @split_name=split(/\./, $file);
	my $sp1=$split_name[0];
	my $sp2=$split_name[1];
	#my $lines =`wc -l $file`;

    #my $ORTHOFINDER = $ARGV[0];
    open(my $filein, "<", $file)   or die "Could not open $file\n";
    my $maxcount=0;
    my $maxbiggest=0;
    my @lengths;
    my $lines=0;
    while (my $line=<$filein>){
        chomp $line;
        if ($line =~ m/^#/){
            #new orthoblock
            if ($maxcount >= $maxbiggest){
                $maxbiggest = $maxcount;
            }
            push (@lengths, $maxcount);
            $maxcount=0;
        }
        else{
            $lines++;
            $maxcount++;
        }
    }

    my $sum=0;
    for my $each (@lengths) {
        $sum += $each;
    }

    my $average_length=$sum/scalar(@lengths);
	print $outhandle "$sp1\t$sp2\t$lines\t$maxbiggest\t$average_length\n";
}


sub average {
    my $ArrayAverage = @_;
    for (@_) {
        my $total += $_;
    }

    return $ArrayAverage;
} 