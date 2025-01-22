use warnings;
use strict;



######################################################################################

###A sorter of syntenic junction breaks###

#This script attempts to sort the anchors file by order of that on each chromosome in species 1. 
#It also attempts to resolve situations where there are two or more of the same gene in various blocks.

# Script written by Chris Wyatt (UCL, with MIT license), 18 January 2025.

######################################################################################


# Ensure the correct number of arguments are provided
if (@ARGV != 1) {
    die "Usage: $0 <anchors_file>\n";
}

my $anchors = $ARGV[0];
my @split_anchor=split(/\./, $anchors);

#Get species names:
my $species1=$split_anchor[0];
my $species2=$split_anchor[1];

#get bed files that should be in same directory
my $bed_sp1 = "$species1\.bed";
my $bed_sp2 = "$species2\.bed";

#Read in bed file of species 1 (bed file is expected to be sorted):
open(my $bed1, "<", $bed_sp1)   or die "Could not open $bed_sp1 \n";

#Read in bed file of species 2 (bed file is expected to be sorted):
open(my $bed2, "<", $bed_sp2)   or die "Could not open $bed_sp2 \n";

print "Sort: We are comparing $species1 versus $species2\n";

# Read in the bed files and store for each gene the chromosome, start and stop, plus the position in the genome
my %sp1_coordinate;
my %sp1_positions;
my %sp2_coordinate;
my %sp2_positions;
my %sp2_coord_hits; # a store of positions in sp2 that are in anchors.
my %sp2_coord_to_gene; # check coordiante to gene name
my %sp1_coord_to_gene;

#Save the output to a file:
my $outfile="$anchors\_anchors_sorted.tsv";
open(my $out, ">", $outfile)   or die "Could not open $outfile \n";
#print $out "Break_number\tInter_Intra_chromosome\tSpecies1\tSpecies2\tsp1_gene_b4\tchr\tsta\tend\tsp1_gene_after\tchr\tsta\tend\tsp2_gene_b4\tchr\tsta\tend\tsp2_gene_after\tchr\tsta\tend\tsp1_genepos_before\tsp1_genepos_after\tgenes_in_gap_sp1\ttrend_before_sp2\ttrend_after_sp2\tsp2_genepos_before\tsp2_genepos_after\tgenes_in_gap_sp2\tfinal_classification\tgene_names_in_sp1_gap\n";
my $outfile2="$anchors\_anchors_sorted_full.tsv";
open(my $out2, ">", $outfile2)   or die "Could not open $outfile2 \n";


#Save Bed file of species 1 into hashes:

my $line_number_sp1=1;
my $last_chromo_sp1="Will_list_chromsomes_when_they_are_read_in";
my %sp1_gene_order;

while ( my $line1 = <$bed1> ){
     chomp $line1;
     my @sp=split("\t", $line1);
     my $chr=$sp[0];
     my $sta=$sp[1];
     my $end=$sp[2];
     my $gen=$sp[3];
     #my $cen=($end+$sta)/2;
     if ($sp1_gene_order{$chr}{$sta}){
     	#This should not normally happen. 
     	#But add a small amount so they do not overlap (nad break the pipeline). 
     	$sta=$sta+0.01;
     	$sp1_gene_order{$chr}{$sta}=$gen;
     	#print "ORDER $cen $gen\n";
     }
     else{
     	$sp1_gene_order{$chr}{$sta}=$gen;
     	#print "$chr $cen $gen\n";
     }
     
     #Check if chromosome changes in bed, and reset the line number count, so that each chromsome has a separate count.
     if ($last_chromo_sp1 eq $chr){
          #Do nothing, all good
     }
     else{
          $line_number_sp1=1;
     }

     # Basic score of position based on where the gene starts 
     # This score is not ideal, as genes can overlap. Sometimes, a gene start-stop overlaps entire genes. 
     # This could be a annotation error, or real, but it messes up the calculation of gene order.
     $sp1_positions{$gen}=$line_number_sp1;

     #Add the coordinate, start and stop of each gene:
     $sp1_coordinate{$gen}="$chr\t$sta\t$end";

     #Used to know what genes are associated with each coord
     $sp1_coord_to_gene{$line_number_sp1}=$gen;

     #print "$species1\t$gen\t$chr\t$line_number_sp1\n";

     $last_chromo_sp1=$chr;
     $line_number_sp1++;
}

#Save Bed file of species 2 into hashes:

my $line_number_sp2=1;
my $last_chromo_sp2="Will_list_chromsomes_when_they_are_read_in";

while ( my $line2 = <$bed2> ){
     chomp $line2;
     my @sp=split("\t", $line2);
     my $chr=$sp[0];
     my $sta=$sp[1];
     my $end=$sp[2];
     my $gen=$sp[3];
     
     #Check if chromosome changes in bed, and reset the line number count, so that each chromsome has a separate count.
     if ($last_chromo_sp2 eq $chr){
          #Do nothing, all good
     }
     else{
          $line_number_sp2=1;
     }

     # Basic score of position based on where the gene starts 
     # This score is not ideal, as genes can overlap. Sometimes, a gene start-stop overlaps entire genes. 
     # This could be a annotation eror, or real, but it messes up the calculation of gene order.
     $sp2_positions{$gen}=$line_number_sp2;

     #Add the coordinate, start and stop of each gene:
     $sp2_coordinate{$gen}="$chr\t$sta\t$end";

     #Used later, check how many genes are in an anchored line:
     $sp2_coord_hits{$line_number_sp2}="y";
     #Used to know what genes are associated with each coord
     $sp2_coord_to_gene{$line_number_sp2}=$gen;
     #print "$line_number_sp2 = y\n";

     #print "$species2\t$gen\t$chr\t$line_number_sp2\n";

     $last_chromo_sp2=$chr;
     $line_number_sp2++;
}


#Read in anchor file:
open(my $anch, "<", $anchors)   or die "Could not open $anchors \n";

#Check for any duplicate entries. Which mess up the calculation of where breaks are. 
#The strategy here is just to remove them later. This could be dealt with in a more sophisticated way in future.
my %has_been_used_sp1;
my %has_been_used_sp2;
my %duplicates_sp1;
my %duplicates_sp2;
my %actual_input_gene_match;
while ( my $line4 = <$anch> ){
	chomp $line4;
	my @split=split("\t", $line4);
	if ($split[0] ne '###'){
		if ($has_been_used_sp1{$split[0]}){
			#print "SP1 this has bee used before\n";
			$duplicates_sp1{$split[0]}="yes";
		}
		if ($has_been_used_sp2{$split[1]}){
			#print "SP2 this has bee used before\n";
			$duplicates_sp2{$split[1]}="yes";
		}
		$has_been_used_sp1{$split[0]}=$line4;
		$has_been_used_sp2{$split[1]}=$line4;
		$actual_input_gene_match{$split[0]}=$split[1];
	}
}
close $anch;



#Read in anchor file again and work out which genes are next to each other in the anchor file:
open(my $anch2, "<", $anchors)   or die "Could not open $anchors \n";

my %anch_store_befor;
my %anch_store_after;
my $last;

while ( my $line3 = <$anch2> ){
	chomp $line3;
	my @split=split("\t",$line3);
	if ($line3 eq '###'){     #Means we have found a break. Number the break and save to array.
		#Do nothing
	}
	else{
		if ($last){
			if ($last eq '###'){
				#Do nothing as last one was a break, so no gene before this.
			}
			else{
				$anch_store_befor{$split[0]}=$last;
				$anch_store_after{$last}=$split[0];
			}
		}
		else{
			#Do nothing, as this is the 2nd line. As last doesn;t exxist at that point.
		}
	}
    $last = $split[0];
}
close $anch2;


#Now go through the bed file in order of genes, and check if the synteny matches, or not.

#Read in bed file of species 1 (bed file is expected to be sorted):
close $bed1;
open(my $bed1_again, "<", $bed_sp1)   or die "Could not open $bed_sp1 \n";



# Or go through the hash

# Loop through the nested hash
my $curr_gene_befor;
my $count_nis=0;
my $last_call="notset";
foreach my $chr (keys %sp1_gene_order) {
    print "Chromosome: $chr\n";
    my $count=0;
    foreach my $position (sort { $a <=> $b } keys %{ $sp1_gene_order{$chr} }) {
          my $gene = $sp1_gene_order{$chr}{$position};
          #print "$chr $count\n";
          if ($has_been_used_sp1{$gene}){
	        	#print "$gene\t yes\n";
	        	my $gene_befor=$anch_store_befor{$gene};
	        	if ($gene_befor){
	        		if ($curr_gene_befor){   #if there was a current last gene before, not for first one
		        		if ($gene_befor eq $curr_gene_befor){
		        			#Then yes print next synteny gene in sequence
		        			print $out "$has_been_used_sp1{$gene}\n";
		        			print $out2 "$has_been_used_sp1{$gene}\n";
		        			$curr_gene_befor=$gene;
		        			$last_call="$has_been_used_sp1{$gene}";
		        		}
		        		else{
		        			print $out2 "$gene\t$actual_input_gene_match{$gene}\tNOT_IN_ORDER\n";
		        			#print "###\n";
		        			$curr_gene_befor=$gene;
		        			$count_nis++;
		        			if ($last_call eq "###"){
		        				#do nothing, don't want to block two breaks in a row.
		        			}
		        			else{
		        				print $out "###\n";
		        				print $out2 "###\n";
		        			}
		        			$last_call="###";
		        		}	
				}
				else{
					#If there is no gene in the anchors before a gene in the bed, we consider this the start of a new block
					
					if ($last_call eq "###"){
	        			#do nothing, don't want to block two breaks in a row.
	        		}
	        		else{
						print $out "###\n";
						print $out2 "###\n";
					}
					
					if ($duplicates_sp1{$gene}){
						#Do nothing with a duplicated syntenic gene, as these are hard to place.
					}
					else{
						#print line in synteny file for this gene:
						print $out "$has_been_used_sp1{$gene}\n";
						print $out2 "$has_been_used_sp1{$gene}\n";
						$curr_gene_befor=$gene;
						$last_call="$has_been_used_sp1{$gene}";
					}
				}
			}
			else{
				#No need to reset these, as else it will incorporate "fake" breaks.
				#$curr_gene_befor=$gene;
				#$last_call="$has_been_used_sp1{$gene}";
			}
        	}
        	else{
        		print $out2 "$gene\tNA\tNA\n";
        		#print "NA NA $gene\n";
        	}
        	$count++;
     }
}

#print $out "###\n";     #print my final break line
#print $out2 "###\n";    #print my final break line

print "my not in seq: $count_nis\n";