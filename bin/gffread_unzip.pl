#!/usr/bin/perl
use warnings;
use strict;

die "Needs the input sample_id, gff3 file and the fasta file\n" if (@ARGV!=3); 

#print "Script is running\n";

my $sample=$ARGV[0];

if ($ARGV[1] =~ m/.gz$/){
	`zcat $ARGV[1] > genome.fa`;
}
else{
	`cp  $ARGV[1] genome.fa`
}

if ($ARGV[2] =~ m/.gz$/){
	`zcat $ARGV[2] > sample.gff3`;
	`cp sample.gff3 $sample\.gff_for_jvci.gff3`
}
else{
	`cp $ARGV[2] $sample\.gff_for_jvci.gff3`

}

my $gff_head=`head -n 1 $sample\.gff_for_jvci.gff3 | cut -f 2`;
chomp $gff_head;
if ($gff_head=~  m/AUGUSTUS/g){
	`GFF3_to_SingleLongest_mRNA.pl $sample\.gff_for_jvci.gff3`;
	`cp $sample\.gff_for_jvci.gff3\_single.longest.gene.gff3 $sample\.gff_for_jvci.gff3`;
}

`gffread -w $sample\.nucl.fa -g genome.fa $sample\.gff_for_jvci.gff3`
