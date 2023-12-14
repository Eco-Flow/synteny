#!/usr/bin/perl
#Bed_to_histogram.pl
use warnings;
use strict;

die "Provide the name of the Chromopaint.txt file\n" if (@ARGV!=3); 

print "Script is running\n";

my $BED_file1 = $ARGV[0];
open(my $BED1, "<", $BED_file1) or die "Could not open $BED_file1\n";
my $BED_file2 = $ARGV[1];
open(my $BED2, "<", $BED_file2) or die "Could not open $BED_file2\n";
my $ANCHORS = $ARGV[2];
open(my $ANCH, "<", $ANCHORS) or die "Could not open $ANCHORS\n";

my $outfile="seqids_karyotype.txt";
open(my $out, "> $outfile") or die "error opening $outfile. $!";

#read in the first bed file
my %SP1;
while (my $line=<$BED1>) {
    chomp $line;
	my @spl=split ("\t", $line);
    $SP1{$spl[0]}={$spl[3]};
}

#read in the secon bed file
my %SP2;
while (my $line2=<$BED2>) {
    chomp $line2;
	my @spl=split ("\t", $line2);
    $SP2{$spl[0]}={$spl[3]};
}

#Calc the scaffold name in the anchors file
my %chromo_pairs;

while (my $line3=<$ANCH>) {
    chomp $line3;
    my @spl=split ("\t", $line3);
    my $len=scalar(@spl);
    if ( $len == 3 ){
        my $chr1=$SP1{$spl[0]};
        my $chr2=$SP2{$spl[1]};
        my $join="$chr1\_$chr2";
        $chromo_pairs{"$join"}="hit";
    }
    else{
        #probably a # start line between syntenic regions.
    }
}

=cut

foreach my $chrome_sp1 ( keys %SP1) {
    my $best_score=0;
    foreach my $chrome_sp2 ( keys % { $SP1{$chrome_sp1} } ) {
        my $score = $SP1{$chrome_sp1}{$chrome_sp2};
        if ($score >= $best_score){
            $best_chromo{$chrome_sp1}=$chrome_sp2;
            $best_score=$score;
            #print "$chrome_sp1 to $chrome_sp2  score $best_score\n";
        }
    }
}

#Work out new chromo name based on the previous data.
my %done;
my %Final_names;
foreach my $chrome_sp1 ( keys %best_chromo) {
    if ($done{$chrome_sp1}){
        $Final_names{$chrome_sp1}="$best_chromo{$chrome_sp1}\_$done{$chrome_sp1}";
        $done{$chrome_sp1}++;
    }
    else{
        $done{$chrome_sp1}=2;
        $Final_names{$chrome_sp1}=$best_chromo{$chrome_sp1};
    }
}


#finally print of the new chromo name hash.

foreach my $sp1 ( keys %Final_names) {
    print $out "$sp1\t$Final_names{$sp1}\n";
}