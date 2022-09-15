#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;


print "Please be in folder with all the anchor files\n";

my @files=`ls *.anchors`;

my $outfile="My_scores.tsv";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
my $outfile2="Synteny_matrix.tsv";
open(my $outhandle2, ">", $outfile2)   or die "Could not open $outfile2 \n";


#Print header to outfile
print $outhandle "Sp1\tSp2\tSyn_lines\tMax_length\tAverage_length\n";

#Initiate hash for matrix store
my %sim_hash;

#Loop thru the files to get the synteny info.
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


    #Calculate a matrix of values (total length of syntenic regions)
    $sim_hash{$sp1}{$sp2}=$lines;
}



#Go thru the synteny line hash and print a matrix
foreach my $sp1 (sort keys %sim_hash ){
    print $outhandle2 "\t$sp1";
}
print $outhandle2 "\n";

foreach my $sp1 (sort keys %sim_hash ){
    print $outhandle2 "$sp1";
    foreach my $sp2 (sort keys %sim_hash ){
        if ($sp1 eq $sp2){
            print $outhandle2 "\tNA";
        }
        else{
            print $outhandle2 "\t$sim_hash{$sp1}{$sp2}";
        }
    }

    print $outhandle2 "\n";
}


sub average {
    my $ArrayAverage = @_;
    for (@_) {
        my $total += $_;
    }

    return $ArrayAverage;
} 