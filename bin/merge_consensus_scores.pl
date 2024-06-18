#!/usr/bin/perl
use warnings;
use strict;

my $sp = "filec";
open(my $INPUT, "<", $sp) or die "Could not open $sp\n";

my $outfile = "filed";
open(my $out, ">", $outfile) or die "Could not open $outfile\n";

#Remove header:
my $head=<$INPUT>;
print $out "$head";

my %pairs;
#Read in the input table, and detect the pairs in the data. Put them into hash entries (%pairs).
while (my $line=<$INPUT>){
    chomp $line;
    my @sp=split("\t", $line);
    my @nsp=split(/\./, $sp[0]);
    #print "$nsp[0] $nsp[1]\n";
    my @sorted = sort @nsp;
    print "$sorted[0]\t$sorted[1]\n";
    my $comb="$sorted[0]\t$sorted[1]";
    if ($pairs{"$comb"}){
    	my $old_l=$pairs{"$comb"};
    	my $new_l="$old_l\n$line";
    	$pairs{"$comb"}=$new_l;
    }
    else{
    	$pairs{"$comb"}=$line;
    }
}

#Go through the pairs of entries, and average all the values.
#Also make sure there are no entries with more or less than 2 entries (would be an error).
#Should work if all species names are unique in input csv.

foreach my $key ( keys %pairs ) {
	my $value=$pairs{$key};
	my @spl=split("\n", $value);
	my @line_run1=split("\t", $spl[0]);
	my @line_run2=split("\t", $spl[1]);
	print $out "$line_run1[0]\t$line_run1[1]\t$line_run1[2]\t$line_run1[3]\t$line_run1[4]\t$line_run1[5]\t$line_run1[6]\t";
	my $avergage_syntenic_blocks=($line_run1[7]+$line_run2[7])/2;
	print $out "$avergage_syntenic_blocks\t";
	my $average_trans_score=($line_run1[8]+$line_run2[8])/2;
	print $out "$average_trans_score\t";
	my $average_inver_score=($line_run1[9]+$line_run2[9])/2;
	print $out "$average_inver_score\t";
	my $average_perc_identity=($line_run1[10]+$line_run2[10])/2;
	print $out "$average_perc_identity\t";
	print $out "$line_run1[11]\t$line_run1[12]\t";
	my $length_transloca_score=($line_run1[13]+$line_run2[13])/2;
	print $out "$length_transloca_score\t";
	my $length_inversion_score=($line_run1[14]+$line_run2[14])/2;
	print $out "$length_inversion_score\t";
	my $TSNL=($line_run1[15]+$line_run2[15])/2;
	print $out "$TSNL\t";
	my $TMN=($line_run1[16]+$line_run2[16])/2;
	print $out "$TMN\t";
	my $IE=($line_run1[17]+$line_run2[17])/2;
	print $out "$IE\t";
	print $out "$line_run1[18]\t";
	my $TJ=($line_run1[19]+$line_run2[19])/2;
	print $out "$TJ\t";
	my $IJ=($line_run1[20]+$line_run2[20])/2;
	print $out "$IJ\t";
	my $SDDJ=($line_run1[21]+$line_run2[21])/2;
	print $out "$SDDJ\t";
	my $LDDJ=($line_run1[22]+$line_run2[22])/2;
	print $out "$LDDJ\t";
	print $out "\n";
}
