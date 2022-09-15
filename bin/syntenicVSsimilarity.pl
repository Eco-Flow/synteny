#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;


print "Please be in folder with all the similarity and synteny files\n";

my $outfile="My_comp_synteny_similarity.tsv";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

my $in1= "Synteny_matrix.tsv";
my $in2= "My_sim_cores.tsv";
open(my $filein1, "<", $in1)   or die "Could not open $in1\n";
open(my $filein2, "<", $in2)   or die "Could not open $in2\n";

my %FILE1;


my $headM=<$filein1>;
chomp $headM;
my @sp=split("\t", $headM);
my $waste=shift(@sp);
my %file1store;

while (my $lineM=<$filein1>){
    chomp $lineM;
    print "line here: $lineM\n";
    my @spline=split("\t", $lineM);
    my $species=shift(@spline);
    my $n=0;
    foreach my $score (@spline){
        my $sp1=$sp[$n];
        my $sp2=$species;
        print "1 $sp1 2 $sp2 3 $score\n";
        $file1store{$sp1}{$sp2}=$score;
        $n++;
    }
}

my $headSy=<$filein2>;
chomp $headSy;
my @sp2=split("\t", $headSy);
my $waste2=shift(@sp2);
my %file2store;

while (my $lineSy=<$filein2>){
    chomp $lineSy;
    my @spline2=split("\t", $lineSy);
    my $species=shift(@spline2);
    my $n=0;
    foreach my $score (@spline2){
        my $sp1=$sp[$n];
        my $sp2=$species;
        $file2store{$sp1}{$sp2}=$score;
        $n++;
    }
}



print $outhandle "$headSy\n";

#Now compare the two numbers:
#print $outhandle "\t";

foreach my $species1 (@sp){
    print $outhandle "$species1";
    foreach my $species2 (@sp){
        print $outhandle "\t$file1store{$species1}{$species2}\/$file2store{$species1}{$species2}";
    }
    print $outhandle "\n";
}
