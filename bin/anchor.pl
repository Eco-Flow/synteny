#!/usr/bin/perl
#Bed_to_histogram.pl
use warnings;
use strict;

die "need three inputs :sp1 bed file (followed by a space), sp2 bed file, followed by the anchor file\n" if (@ARGV!=3); 

print "Script is running\n";

my $BED_file1 = $ARGV[0];
open(my $BED1, "<", $BED_file1) or die "Could not open $BED_file1\n";
my $BED_file2 = $ARGV[1];
open(my $BED2, "<", $BED_file2) or die "Could not open $BED_file2\n";
my $anchor=$ARGV[2];
open(my $ANC, "<", $anchor) or die "Could not open $anchor\n";

my $outfile="Chromopaint.txt";
open(my $out, "> $outfile") or die "error opening $outfile. $!";


#read in the colour hex list
my $colours="unique_hex2";
open(my $COL, "<", $colours) or die "Could not open $colours\n";
my $col_outfile="colour.idmap";
open(my $col_out, ">", $col_outfile) or die "Could not open $col_outfile\n";

my @colour_board;
while (my $line4=<$COL>) {
    chomp $line4;
	push(@colour_board, "$line4");
}

#read in the first bed file
my %SP1;
while (my $line=<$BED1>) {
    chomp $line;
	my @spl=split ("\t", $line);
	$SP1{$spl[3]}="$spl[0]\t$spl[1]\t$spl[2]";
    #print "$spl[3]\n";
}


#print "We got here\n";
#read in the second bed file.
my %SP2;
while (my $line2=<$BED2>) {
    chomp $line2;
	my @spl=split ("\t", $line2);
	$SP2{$spl[3]}="$spl[0]\t$spl[1]\t$spl[2]";
    #print "$spl[3]\n";
}

my $temp_lo;
my $temp_hi;
my $sp1_bed;
my $sp2_bed;
my $chrome1;
my $chrome2;

my %colour_chrome;

while (my $line3=<$ANC>) {
    chomp $line3;
	my @spl=split ("\t", $line3);
    #print "$spl[0]\n";
	


	if ($line3=~ m/###/g){
        if($sp1_bed){
            print $out "$chrome2\t$temp_lo\t$temp_hi\t$chrome1\t0\t\+\n";
            $colour_chrome{"$chrome1"}="done";
            $sp1_bed=();
            $temp_lo=();
            $temp_hi=();
            $sp2_bed=();
            $chrome1=();
            $chrome2=();
        }
	}
	else{
		$sp1_bed=$SP1{$spl[0]};
        #print "Here $spl[0]  :  $sp1_bed\n";
		my @split_1=split("\t", $sp1_bed);
		$chrome1=$split_1[0];

        #print "$spl[1]\n";
		$sp2_bed=$SP2{$spl[1]};
		my @split_2=split("\t", $sp2_bed);
		$chrome2=$split_2[0];
		my $start=$split_2[1];
		my $end=$split_2[2];
        #print "$temp_lo > $start\n";
		if ($temp_lo){
            if($temp_lo > $start){
			    $temp_lo=$start;
		    }
        }
        else{
            $temp_lo=$start;
        }

        if ($temp_hi){
            if ($temp_hi < $end){
			    $temp_hi=$end;
		    }
        }
		else{
            $temp_hi=$end;
        }
	}	
}

print "Script has completed\n";

#Now make the chrome idmap

my $count=0;
foreach my $chr_id ( keys %colour_chrome) {
    print $col_out "$chr_id\t$chr_id\t$colour_board[$count]\n";
    $count++;
}
