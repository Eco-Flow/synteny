#!/usr/bin/perl
use warnings;
use strict;

die "Please specify (1) gff3 file\n" unless(@ARGV==1);
print "Starting script\n\n";

my $gff = $ARGV[0];
my $outfile= "$gff\_single.longest.gene.gff3";

open(my $IN, "<", $gff)   or die "Could not open $gff \n";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

my %gene_hash;
my %best;
my $gene_line;
my $mrna_line;

my %mRNA_hash;

#Set variables that may change during loops:
my $gene;
my $mRNA;
my %mRNA_lines;
my $firstline=0;

while (my $line=<$IN>){
	chomp $line;
	my @split=split("\t", $line);
	my @col9split=split("\;", $split[8]);
	my $type=$split[2];
	print "$line\n";

	if ($firstline == 0){
		#first gene line is entered, but not assessed, as this is the first entry, no gene before to summarise (next step)
		if ($type eq "gene"){
			$gene_line=$line;
			my @ID_split=split("\=", $col9split[0]);
			$gene = $ID_split[1];
			print "gene = $gene\n";
			$firstline = 1;
		}
	}
	else{
		if ($type eq "gene"){
		
			
			my @ID_split=split("\=", $col9split[0]);
			$gene = $ID_split[1];
			#print "gene = $gene\n";
			#If the gene is new:
			if ($gene_hash{$gene}){
				#Fine
				print "Shouldnt happen\n";
			}
			else{
				#This is a new gene, so we can determine the longest isoform of last entry and print it.
				my $longest=0;
				#print "Im here\n";

				my $bestofall;
				foreach my $key (keys %gene_hash ){
					#print "KEY $key\n"; # line should only pull last gene entry
					foreach my $mRNAs_in_gene (keys %{$gene_hash{$key}} ){
						#print "WHATTTTT $key -> $mRNAs_in_gene $gene_hash{$key}{$mRNAs_in_gene}\n";
						my $length=$gene_hash{$key}{$mRNAs_in_gene};
						#print "$length >= $longest\n";
						if ($length >= $longest){
							#print "HEREEEEE: $length >= $longest $mRNAs_in_gene \n ";
							$best{$key}=$mRNAs_in_gene;
							$bestofall=$mRNAs_in_gene;
							print "Best transcript \= $key $mRNAs_in_gene $length\n";
							$longest=$length;
						}
					}
				}



				
					
				


				#For each best mRNA in gene print.

				print $outhandle "$gene_line\n";
				print $outhandle "$mRNA_hash{$bestofall}\n";
				print $outhandle "$mRNA_lines{$bestofall}\n";
				


				foreach my $key (keys %gene_hash ){
					foreach my $mRNAs_in_gene2 (keys %{$gene_hash{$key}} ){
						delete($gene_hash{$key}{$mRNAs_in_gene2});
						print "$key $mRNAs_in_gene2 to be dleeted from hash\n";
					}
				}
			}

			$gene_line=$line;

		}

		if ($type eq "mRNA"){

			$mrna_line=$line;
			my @Parent_split=split("\=", $col9split[1]);
			my @ID_split=split("\=", $col9split[0]);
			$mRNA = $ID_split[1];
			my $parent = $Parent_split[1];
            $mRNA_hash{$mRNA}=$line;
		}

		if ($type eq "exon"){

			my @Parent_split=split("\=", $col9split[1]);
			my $parent = $Parent_split[1];
			my $exon_len=$split[4]-$split[3];
			print "gene $gene parent $parent exon len $exon_len\n";
			if ($gene_hash{"$gene"}{"$parent"}){
				my $old=$gene_hash{$gene}{$parent};
				$gene_hash{$gene}{$parent}=$exon_len+$old;
			}
			else{
				$gene_hash{$gene}{$parent}=$exon_len;
			}
			
		}


		#Add all other mRNA lines into another hash:
		if ($type eq "gene" || $type eq "mRNA"){
			#do nothing.
		}
		else{
			if ($mRNA_lines{$mRNA}){
				my $current=$mRNA_lines{$mRNA};
				my $new="$current\n$line";	
				$mRNA_lines{$mRNA}=$new;
			}
			else{
				$mRNA_lines{$mRNA}=$line;
			}
		}
	}	

	#print $outhandle "$LOC_ID\t$gen_ID\n";
}

