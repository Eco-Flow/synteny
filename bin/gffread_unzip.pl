use warnings;
use strict;

die "Needs the input sample_id, gff3 file and the fasta file\n" if (@ARGV!=3); 

#print "Script is running\n";

my $sample=$ARGV[0];

if ($ARGV[1] =~ m/.gz$/){
	`zcat $ARGV[1] > genome.fa`;
}
else{
	`cp  $ARGV[1] genome.fa`;
}

if ($ARGV[2] =~ m/.gz$/){
	`zcat $ARGV[2] > sample.gff3`;
	`cp sample.gff3 $sample\.gff_for_jvci.tmp.gff3`;
	`sed 's/?/./' $sample\.gff_for_jvci.tmp.gff3 > $sample\.gff_for_jvci.tmp2.gff3`; #Replace any ? with ., as gffread will not work when a ? is present
}
else{
	`cp $ARGV[2] $sample\.gff_for_jvci.tmp.gff3`;
	`sed 's/?/./' $sample\.gff_for_jvci.tmp.gff3 > $sample\.gff_for_jvci.tmp2.gff3`; #Replace any ? with ., as gffread will not work when a ? is present
}

#Now check the gff3 does not have a start position greater than an end position. Cols 4 and 5. This breaks gffread
my $outfile="$sample\.gff_for_jvci.gff3";
open(my $out, ">", $outfile)   or die "Could not open $outfile \n";
#Read in GFF:
my $inname="$sample\.gff_for_jvci.tmp2.gff3";
open(my $in1, "<", $inname)   or die "Could not open $inname \n";
while ( my $line1 = <$in1> ){
	chomp $line1;
	my @split=split("\t", $line1);
	if ($split[4]){
		if ($split[3] > $split[4]){
			#dont print, this line is not right, and it will fail at gffread
		}
		else{
			print $out "$line1\n";
		}
	}
}


my $gff_head=`head -n 1 $sample\.gff_for_jvci.gff3 | cut -f 2`;
chomp $gff_head;
if ($gff_head=~  m/AUGUSTUS/g){
	`GFF3_to_SingleLongest_mRNA.pl $sample\.gff_for_jvci.gff3`;
	`cp $sample\.gff_for_jvci.gff3\_single.longest.gene.gff3 $sample\.gff_for_jvci.gff3`;
}

`gffread -w $sample\.nucl.fa -g genome.fa $sample\.gff_for_jvci.gff3`;
