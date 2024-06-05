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
my $actual_previous_2nd;

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
			chomp $firstline;
			chomp $seconline;
			print "FIRST $firstline\n";

			#Match up to previous order. 
			my %prev_hash;
			my $hash = "$previous\.$_\.pairs.txt";
			open(my $PREHASH, "<", $hash) or die "Could not open $hash\n";
			while (my $line= <$PREHASH>){
				chomp $line;
				my @sp=split("\t", $line);
				$prev_hash{$sp[0]}=$sp[1];
			}

			my @value_2nd=split("\,", $actual_previous_2nd);
			my @value_3rd;
			my %done;
			foreach my $order_2nd (@value_2nd){
				print "2nd line values $order_2nd\n";
				if ($prev_hash{$order_2nd}){
					#if we have a match to the previous chromosome number, good,. use it.
					if ($done{"$prev_hash{$order_2nd}"}){
						#already done, don't push
					}
					else{
						push (@value_3rd, $prev_hash{$order_2nd});
						$done{"$prev_hash{$order_2nd}"}="YES";
					}
				}
				else{
					#Else, ignore it, we don't have an equivalent. 
				}
			}

			#Add on any missing chromosomes that are unique to next species:
			my @second_species_all_chromosomes=split("\,", $seconline);

			foreach my $chros (@second_species_all_chromosomes){
				if ($done{"$chros"}){
					#fine its done
				}
				else{
					push (@value_3rd, $chros);
				}
			}


			my $join_new_3rd=join("\,", @value_3rd);

			print $out2 "$join_new_3rd\n";
			$actual_previous_2nd=$join_new_3rd;
			close $SPEC;
			close $PREHASH;
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
			$actual_previous_2nd=$seconline;
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

