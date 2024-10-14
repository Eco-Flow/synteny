use warnings;
use strict;
use List::Util qw(max);
use List::Util qw(min);

my @breaks=`ls *Break_junction_information.txt`;

#Summary of all runs:
my $outname="Trans_Inversion_junction_count.txt";
open(my $OUT, ">", $outname)   or die "Could not open $outname\n";
#Print header of output file- see after 236 where we print this line
print $OUT "Comparison\tTranslocation_junctions\tInversion_junctions\tSame_direction_duplication_junctions\tLoop_direction_duplication_junctions\n";

foreach my $file (@breaks){
	chomp $file;
	#Sort out which name comparison we are making
	my @split_name=split(/\./, $file);
	my $species1=$split_name[0];
	my $species2=$split_name[1];
	my $combName="$species1\.$species2";

	my $outname2="$combName\.Translocation_gene_boundaries.txt";
	open(my $OUT2, ">", $outname2)   or die "Could not open $outname2\n";
	my $outname3="$combName\.Inversion_gene_boundaries.txt";
	open(my $OUT3, ">", $outname3)   or die "Could not open $outname3\n";
	my $outname4="$combName\.Other_gene_boundaries.txt";
	open(my $OUT4, ">", $outname4)   or die "Could not open $outname4\n";

	#Initiate out file handle
	my $outfile="$combName\.Classification_summary.tsv";
	open(my $out, ">", $outfile)   or die "Could not open $outfile \n";

	#Read in syntenic anchor files (produced by MScanX using the program jcvi)
	my $in_break_file="$combName\.Break_junction_information.txt";
	open(my $in_break, "<", $in_break_file)   or die "Could not open $in_break_file\n";

	#Store for each syntenic blocks gene order in sp2. Needed for next script.
	my $in_sp2_order="$combName\.Sp2_synteny_order.txt";
	open(my $in_order, "<", $in_sp2_order)   or die "Could not open $in_sp2_order\n";

	#First read in the synteny order file so we have a hash that tells us exactly where all the genes are location. Based on their order.

	my %gene_order_all;
	my %gene_order_max;
	my %gene_order_min;
	my %position_to_syntenic_block;
	while ( my $line1 = <$in_order> ){
		chomp $line1;
		my @sp=split("\t", $line1);
		my $chr=$sp[0];
		my $syn=$sp[1];
		my $list=$sp[2];
		my @array_list=split("\,", $list);

		if ($line1 =~ m/Synteny_block_number/g){
			#ignore this line, its the header
		}
		else{
			$gene_order_all{$chr}{$syn}=$list;
			$gene_order_max{$chr}{$syn}=max(@array_list);
			$gene_order_min{$chr}{$syn}=min(@array_list);
		}

		#then loop through the genes in each block to assign each position to a syntenic block number
		foreach my $gene_pos (@array_list){
			$position_to_syntenic_block{$chr}{$gene_pos}=$syn;
		}
	     
	}


	my $removeheader=<$in_break>;
	my $odd_junction=0;
	my $same_direction=0;
	my $inversions=0;
	my $translocations=0;

	#For each my apparent syntenic block.
	#Check first if the gene in syntenic block 1, from range min to max, are found within any other sytenic block on the same chromosome.
	while ( my $line2 = <$in_break> ){
		chomp $line2;
		#print "START\n$line2\n";
		my @sp=split("\t", $line2);
		my $chr=$sp[0];
		my $syn1=$sp[1];
		my $syn2=$sp[2];

		my $sta_1=$sp[3];
		my $end_1=$sp[4];
		my $sta_2=$sp[5];
		my $end_2=$sp[6];

		my $type_original=$sp[7];

		#calcualte direction of block
		my $direction_1;
		if ($sta_1 == $end_1){
			$direction_1=0;
		}
		elsif($sta_1 < $end_1){
			$direction_1=1;
		}
		else{
			$direction_1=-1;
		}

		my $direction_2;
		if ($sta_2 == $end_2){
			$direction_2=0;
		}
		elsif($sta_2 < $end_2){
			$direction_2=1;
		}
		else{
			$direction_2=-1;
		}

		
		my $min_block1=$gene_order_min{$chr}{$syn1};
		my $max_block1=$gene_order_max{$chr}{$syn1};
		my $min_block2=$gene_order_min{$chr}{$syn2};
		my $max_block2=$gene_order_max{$chr}{$syn2};

		#try to calculate the smallest gap between the two blocks.
		my $gap;
		my $minmax=min($sta_1,$end_1)-max($sta_2,$end_2);
		my $maxmin=max($sta_1,$end_1)-min($sta_2,$end_2);
		if ($minmax <= $maxmin){
			$gap=$minmax;
		}
		else{
			$gap=$maxmin;
		}
		
		

		if ($type_original eq "INVER"){
			my $junction_type;
			#try to detect if its an inversion:
			if ($direction_1==0 || $direction_2==0){
				$odd_junction++;
				$junction_type="Odd-duplicate";
			}
			elsif($direction_1==1 && $direction_2==1){
				$same_direction++;
				$junction_type="Same_direction";
				print $OUT4 "$sp[-1]\n$sp[-2]\n";
			}
			elsif($direction_1==0 && $direction_2==0){
				$same_direction++;
				$junction_type="Same_direction";
				print $OUT4 "$sp[-1]\n$sp[-2]\n";
			}
			else{
				#rest must be inversions:
				$junction_type="Inversion";
				$inversions++;
				print $OUT3 "$sp[-1]\n$sp[-2]\n";
			}
			print $out "$chr\t$syn1\t$syn2\t$min_block1\t$max_block1\t$sta_1\t$end_1\t$sta_2\t$end_2\t$gap\t$junction_type\n";
		}
		else{
			#Then its a translocation, we cannot calculate a gap
			print $out "$chr\t$syn1\t$syn2\t$min_block1\t$max_block1\t$sta_1\t$end_1\t$sta_2\t$end_2\tNA\tTranslocation\n";
			$translocations++;

			#Save genes that are on the boundary of translocations:
			print $OUT2 "$sp[-1]\n$sp[-2]\n";
		}
	}

	print $OUT "$combName\t$translocations\t$inversions\t$same_direction\t$odd_junction\n";
	print "Translocations: $translocations\nInversions: $inversions\nSame_direction_likely_duplicate: $same_direction\nOdd_direction_likely_duplicates : $odd_junction   \n\#(likely duplicated, one of the syntenic block starts and ends with the same gene)\nSame_direction $same_direction    \n\#Gene order in same direction (so not an inversion [or translocation]), could be caused by duplications\n\n";
	close $out;
	close $in_break;
	close $in_order;
	close $OUT2;
	close $OUT3;
	close $OUT4;
}



