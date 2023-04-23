#!/usr/bin/perl
#Bed_to_histogram.pl
use warnings;
use strict;

die "Provide the name of the Chromopaint.txt file\n" if (@ARGV!=1); 

print "Script is running\n";

my $BED_file1 = $ARGV[0];
open(my $BED1, "<", $BED_file1) or die "Could not open $BED_file1\n";

my $outfile="Chromo_equivalent.txt";
open(my $out, "> $outfile") or die "error opening $outfile. $!";

#read in the first bed file
my %SP1;
while (my $line=<$BED1>) {
    chomp $line;
	my @spl=split ("\t", $line);
    my $size = $spl[2]-$spl[1];

    if ( $SP1{$spl[0]}{$spl[3]} ){
        my $old=$SP1{$spl[0]}{$spl[3]};
        $SP1{$spl[0]}{$spl[3]}=$size+$old;
    }
    else{
        $SP1{$spl[0]}{$spl[3]}=$size;
    }

}

#Calculate the best chromo hit for each chromo in species 1.
my %best_chromo;
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