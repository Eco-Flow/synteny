#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Basic qw(:all);

my $sp = "species.csv";
open(my $SPEC, "<", $sp) or die "Could not open $sp\n";

my $outfile = "layout_all";
open(my $out, ">", $outfile) or die "Could not open $outfile\n";
print $out "# y, xstart, xend, rotation, color, label, va,  bed\n";

my $outseq = "seqids_karyotype_all.txt";
open(my $out2, ">", $outseq) or die "Could not open $outseq\n";

my $line=<$SPEC>;

my @species=split(",", $line);

my $len=scalar(@species);

# Run through the species in the order they are presented, and run the set up for the karyotype plots. 
my @store_comparisons;
my $n=0;
my $p=90;
my $intervals=int(80/$len);

my $previous;
for (@species){
	chomp $_;

	#print "here $_ \n";
	if ($previous){
		#print "$previous vs $_\n";
		#print "exists? $previous\.bed.flipped.bed\n";
		#If a flipped version of chromosomes exists use that!
		if (-e "$previous\.bed.flipped.bed"){
			print "syntenous_chromosomes.pl $previous.bed.flipped.bed $_.bed $previous\.$_\.anchors\n";
			`syntenous_chromosomes.pl $previous.bed.flipped.bed $_.bed $previous\.$_\.anchors`;
			`mv seqids_karyotype.txt $previous\.$_\.seqids.txt`;
			#Store comp line information
			my $plus=$n+1;
			my $comp_line="e, $n, $plus, $previous\.$_\.anchors.simple";
			$n++;
			push (@store_comparisons, "$comp_line");

			#Store chromosome order list
			my $fin = "$previous\.$_\.seqids.txt";
			open(my $SPEC, "<", $fin) or die "Could not open $fin\n";
			my $firstline=<$SPEC>;
			my $seconline=<$SPEC>;
			chomp $seconline;
			print $out2 "$seconline\n";
			close $SPEC;
		}
		else{
			print "syntenous_chromosomes.pl $previous.bed $_.bed $previous\.$_\.anchors\n";
			`syntenous_chromosomes.pl $previous.bed $_.bed $previous\.$_\.anchors`;
			`mv seqids_karyotype.txt $previous\.$_\.seqids.txt`;
			#Store comp line information
			my $plus=$n+1;
			my $comp_line="e, $n, $plus, $previous\.$_\.anchors.simple";
			$n++;
			push (@store_comparisons, "$comp_line");

			#Store chromosome order list
			my $fin = "$previous\.$_\.seqids.txt";
			open(my $SPEC, "<", $fin) or die "Could not open $fin\n";
			my $firstline=<$SPEC>;
			my $seconline=<$SPEC>;
			chomp $firstline;
			chomp $seconline;
			print $out2 "$firstline\n$seconline\n";
			close $SPEC;
		}
		
		
		$p=$p-$intervals;
		print $out " .$p,     .1,    .8,      0,      , $_, top, $_.bed.flipped.bed\n";

		#Make sure the previous is now set:
		$previous="$_";
	}
	else{

		print $out " .9,     .1,    .8,      0,      , $_, top, $_.bed\n";
		#No previous to make a consecutive comparison.
		$previous=$_;
	}
}


#Now prit rest of layout file:

print $out "# edges\n";

foreach my $comps (@store_comparisons){
	print $out "$comps\n";
}

#print "$len\t$species[0]\n";

close $out;
close $out2;

