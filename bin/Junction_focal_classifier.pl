use warnings;
use strict;
use List::Util qw(max);
use List::Util qw(min);


######################################################################################

###A classifier of syntenic junction breaks###

#This script is a second method written to provide basic counts of junctions in the 
#genome that appear to be related to inversions or translocations. 
#This is sometimes difficult to ascertain, due to duplications or unexpected gaps that 
#can break synteny, but in fact are not linked to any structural change such as an 
#inversion (a reverse oriented region without a chromosome) or translocation (movement 
#between chromosomes)

# Input is the anchors file between two species, along with their bed files.

# Script written by Chris Wyatt (UCL, with MIT license), 5 December 2024.

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

my $bed_sp1 = "$species1\.bed";
my $bed_sp2 = "$species2\.bed";

print "We are comparing $species1 versus $species2\n";

# Read in the bed files and store for each gene the chromosome, start and stop, plus the position in the genome
my %sp1_coordinate;
my %sp1_positions;
my %sp2_coordinate;
my %sp2_positions;
my %sp2_coord_hits; # a store of positions in sp2 that are in anchors.
my %sp2_coord_to_gene; # check coordiante to gene name
my %sp1_coord_to_gene;

#Read in bed file of species 1 (bed file is expected to be sorted):
open(my $bed1, "<", $bed_sp1)   or die "Could not open $bed_sp1 \n";

#Read in bed file of species 2 (bed file is expected to be sorted):
open(my $bed2, "<", $bed_sp2)   or die "Could not open $bed_sp2 \n";

#Save the output to a file:
my $outfile="$species1\.$species2\_junction_details.tsv";
open(my $out, ">", $outfile)   or die "Could not open $outfile \n";
print $out "Break_number\tInter_Intra_chromosome\tSpecies1\tSpecies2\tsp1_gene_b4\tchr\tsta\tend\tsp1_gene_after\tchr\tsta\tend\tsp2_gene_b4\tchr\tsta\tend\tsp2_gene_after\tchr\tsta\tend\tsp1_genepos_before\tsp1_genepos_after\tgenes_in_gap_sp1\ttrend_before_sp2\ttrend_after_sp2\tsp2_genepos_before\tsp2_genepos_after\tgenes_in_gap_sp2\tfinal_classification\tgene_names_in_sp1_gap\n";

#Save summary to file too:
my $outfile_sum="$species1\.$species2\_junction_summary.tsv";
open(my $out_sum, ">", $outfile_sum)   or die "Could not open $outfile_sum \n";
print $out_sum "Type\tCount\n";
my $inter_c=0;
my $inver_c=0;
my $indel_c=0;
my $indel_c0_5=0;
my $indel_c5_20=0;
my $indel_c20_more=0;

#Starting Script Proper

#Save Bed file of species 1 into hashes:

my $line_number_sp1=1;
my $last_chromo_sp1="Will_list_chromsomes_when_they_are_read_in";

while ( my $line1 = <$bed1> ){
     chomp $line1;
     my @sp=split("\t", $line1);
     my $chr=$sp[0];
     my $sta=$sp[1];
     my $end=$sp[2];
     my $gen=$sp[3];
     
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
     $sp1_coord_to_gene{$chr}{$line_number_sp1}=$gen;
     #print "$line_number_sp1 equals $gen\n";

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


# Now we read in the anchors file, and start checking where are the gaps, 
# and later,,, what genes are within the break, and check what genes are in the run up to the break and after it.


#Read in bed file of species 1 (bed file is expected to be sorted):
open(my $anch, "<", $anchors)   or die "Could not open $anchors \n";

my $remove_first_hash=<$anch>; #We can safely remove the first ###, as this starts the file. 
my @line_store;   #Here we will store the contents of the previous lines, in order they were found
my $line_number_anch=0;
my $break_number=1;
my %break_locations;

while ( my $line3 = <$anch> ){
     chomp $line3;

     if ($line3 eq '###'){     #Means we have found a break. Number the break and save to array.
          push (@line_store, "$line_number_anch\t$break_number");
          $break_locations{"$break_number"}=$line_number_anch;
          $break_number++;
     }
     else{
          my @sp=split("\t", $line3);
          my $sp1_gene=$sp[0];
          my $sp2_gene=$sp[1];
          push (@line_store, "$line_number_anch\t$sp1_gene\t$sp2_gene");
     }
     $line_number_anch++;
}

#Go through each break and start building up the information about before and after genes.
foreach my $breaks (sort { $a <=> $b } keys %break_locations ){ #For each break we found.
     #print to screen which chromosomes are being processed:
     print "$breaks\t$break_locations{$breaks}\n";
     print "$line_store[$break_locations{$breaks}]\n";
     my $break_loc_number=$break_locations{$breaks};
     my $line_befor_break=$line_store[$break_loc_number-1];
     my $line_after_break=$line_store[$break_loc_number+1];
     my $gene_befor_break_sp1=(split /\t/, $line_befor_break)[1];
     my $gene_after_break_sp1=(split /\t/, $line_after_break)[1];
     print "befor $gene_befor_break_sp1 after $gene_after_break_sp1\n";
     my $gene_befor_break_sp2=(split /\t/, $line_befor_break)[2];
     my $gene_after_break_sp2=(split /\t/, $line_after_break)[2];

     # Check the chromosomes before and after the break in sp1 and sp2
     my $chromosome_befor_sp1=(split /\t/, $sp1_coordinate{$gene_befor_break_sp1})[0];
     my $chromosome_after_sp1=(split /\t/, $sp1_coordinate{$gene_after_break_sp1})[0];
     my $chromosome_befor_sp2=(split /\t/, $sp2_coordinate{$gene_befor_break_sp2})[0];
     my $chromosome_after_sp2=(split /\t/, $sp2_coordinate{$gene_after_break_sp2})[0];

     print "Chromosome: $chromosome_befor_sp1 $chromosome_after_sp1\n";
     print "Chromosome: $chromosome_befor_sp2 $chromosome_after_sp2\n";

     if ($chromosome_befor_sp1 ne $chromosome_after_sp1){
          #This means that the break separates two chromosomes in species 1. So it is not a syntenic break, but a span from one chromosome to another. 
          print "Break $breaks : chromosomal break (not syntenic breaks)\n";
     }
     else{

          #Check if there is a gap between two genes in sp1 across the break (and the size of this gap):
          my $gene_befor_pos=$sp1_positions{$gene_befor_break_sp1};
          my $gene_after_pos=$sp1_positions{$gene_after_break_sp1};
          my $diff=($gene_after_pos-$gene_befor_pos)-1;
          #Check what the genes in gap are:

          print "Here are the beginning and end in sp1 gene names: $gene_befor_break_sp1 and $gene_after_break_sp1\n";
          print "with $gene_befor_pos and $gene_after_pos\n";
          my @sp1_gene_names_in_gap;
          for ( my $i=$gene_befor_pos+1; $i<=$gene_after_pos-1; $i++ ){
               print "GAP $i\:$sp1_coord_to_gene{$chromosome_befor_sp1}{$i}\n";
               if ($sp1_coord_to_gene{$chromosome_befor_sp1}{$i}){
                    push (@sp1_gene_names_in_gap, "$sp1_coord_to_gene{$chromosome_befor_sp1}{$i}");
               }
          }
          my $joint_sp1_gap_genes=join("\,", @sp1_gene_names_in_gap);

          #This means that the break separates species 1's chromosome in species 2.
          if ($chromosome_befor_sp2 ne $chromosome_after_sp2){   # This means that the chromosomes in sp2 are different, despite being the same in species 1. Looks like a translocation break type (not a translocation for sure, as it could be that two chromosomes fused, and then lots of inversions happened).
               # TRANSLOCATIONS:
               print "Break $breaks : represents a translocation type break, sp2 regions before and after are different chromosomes\n";
               print $out "$breaks\tInter\t$species1\t$species2\t$gene_befor_break_sp1\t$sp1_coordinate{$gene_befor_break_sp1}\t$gene_after_break_sp1\t$sp1_coordinate{$gene_after_break_sp1}\t$gene_befor_break_sp2\t$sp2_coordinate{$gene_befor_break_sp2}\t$gene_after_break_sp2\t$sp2_coordinate{$gene_after_break_sp2}\t$gene_befor_pos\t$gene_after_pos\t$diff\tNA\tNA\tNA\tNA\tNA\tInter\t$joint_sp1_gap_genes\n";
               $inter_c++;
          }
          else{
               print "Break $breaks : represents a potential inversion of disruption within the same chromosome\n";
               print $out "$breaks\tIntra\t$species1\t$species2\t$gene_befor_break_sp1\t$sp1_coordinate{$gene_befor_break_sp1}\t$gene_after_break_sp1\t$sp1_coordinate{$gene_after_break_sp1}\t$gene_befor_break_sp2\t$sp2_coordinate{$gene_befor_break_sp2}\t$gene_after_break_sp2\t$sp2_coordinate{$gene_after_break_sp2}\t$gene_befor_pos\t$gene_after_pos\t$diff";

               #Now check if the genes increase/decrease before and after (not likely an inversion), or that genes increase before and decrease after (or vice versa), indicating a likely inversion.
               my @pre_sp2_positions;
               my @post_sp2_positions;
               for (my $i=10; $i>=1; $i-- ){
                    my $line_befor_break_10=$line_store[$break_loc_number-$i];
                    print "B $line_befor_break_10";
                    my @array=split("\t", $line_befor_break_10);
                    my $length = scalar @array; 
                    if ($length == 2){ #if we go back and fall into another break, we have to remove all the entries in the result, and continue again.
                         @pre_sp2_positions = ();
                    }
                    else{
                         my @spl_line=split("\t", $line_befor_break_10);
                         #get the positions in sp2:
                         my $pos_in_sp2=$sp2_positions{$spl_line[2]};
                         #add this to the array:
                         push (@pre_sp2_positions, $pos_in_sp2);
                         print "\t$pos_in_sp2\n";
                    }
               }
               print "here\n";
               my $stop_if_reached_hash=1;   # For next loop, we need to stop printing, if we reach a triple hash, the next break
               for (my $i=1; $i<=10; $i++ ){
                    if ($stop_if_reached_hash){
                         my $line_after_break_10=$line_store[$break_loc_number+$i];
                         my @array=split("\t", $line_after_break_10);
                         my $length = scalar @array; 
                         if ($line_after_break_10){
                              print "A $line_after_break_10";
                              if ($length == 2){ #if we go forward and fall into another break, we have to remove all the entries in the result, and continue again.
                                   #Do nothing , dont print to list, and stop any more lines being added.
                                   $stop_if_reached_hash=0;
                                   print "Yes, we found  [ $line_after_break_10 ] \n";
                              }
                              else{
                                   if ($stop_if_reached_hash){    # if we haven't reached the next break yet, then print
                                        my @spl_line=split("\t", $line_after_break_10);
                                        #get the positions in sp2:
                                        my $pos_in_sp2=$sp2_positions{$spl_line[2]};
                                        #add this to the array:
                                        push (@post_sp2_positions, $pos_in_sp2);
                                        print "\t$pos_in_sp2\n";
                                   }
                              }
                         }
                         else{
                              $stop_if_reached_hash=0;
                         }
                    }
               }
               print "and over here\n";

               #check the trend of the genes before and after the break.

               my $trend_befor = check_trend(@pre_sp2_positions);
               my $trend_after = check_trend(@post_sp2_positions);

               print "Before: $trend_befor\n";
               print "After: $trend_after\n";

               print $out "\t$trend_befor\t$trend_after";

               #Check the gap in SP2, as it may be that the sequence continues after the break (often this is due to duplications, messing up the syntenic block).
               # Calculate the sum of elements in each array
               # Calculate the average for @pre_sp2_positions
               my $sum_before = 0;
               $sum_before += $_ for @pre_sp2_positions;
               my $avg_before = @pre_sp2_positions ? $sum_before / @pre_sp2_positions : 0;

               # Calculate the average for @post_sp2_positions
               my $sum_after = 0;
               $sum_after += $_ for @post_sp2_positions;
               my $avg_after = @post_sp2_positions ? $sum_after / @post_sp2_positions : 0;

               # Compare the averages
               my $max;
               my $min;
               if ($avg_before >= $avg_after) {                       
                    $max=min(@pre_sp2_positions);
                    $min=max(@post_sp2_positions);
                    print $out "\t$max\t$min";
               }
               else{
                    $max=min(@post_sp2_positions);
                    $min=max(@pre_sp2_positions);
                    print $out "\t$min\t$max";
               }
               my $diff_in_sp2=$max-$min;    # Calculate the difference between the closest two numbers
               print "Diff in species 2 : $diff_in_sp2  (of $max and $min)\n";
               # Now check if any of the genes in sp2 are in other syntenic regions (in the anchor file).
               my $in_another_anchor=0;
               my @sp2_gene_names_in_gap;
               for ( my $i=$min+1; $i<=$max-1; $i++ ){
                    #print "H   $i\n";
                    if ($sp2_coord_hits{$i}){
                         $in_another_anchor++;
                         push (@sp2_gene_names_in_gap, "$sp2_coord_to_gene{$i}");
                    }
               }
               print "And $in_another_anchor genes between the two segments in species 2 are in other syntenic blocks\n";
               print $out "\t$in_another_anchor";

               if ($trend_befor eq $trend_after){
                    $indel_c++;
                    if ($in_another_anchor <= 5){
                         print $out "\tindel_tiny\t$joint_sp1_gap_genes\n";
                         $indel_c0_5++;
                    }
                    elsif ($in_another_anchor <= 20){
                         print $out "\tindel_small\t$joint_sp1_gap_genes\n";
                         $indel_c5_20++;
                    }
                    else{
                         print $out "\tindel_large\t$joint_sp1_gap_genes\n";
                         $indel_c20_more++;
                    }
                    
               }
               else{
                    print $out "\tinver\t$joint_sp1_gap_genes\n";
                    $inver_c++;
               }
          }
     }
}


#Print out the 

print $out_sum "Inter\t$inter_c\n";
print $out_sum "Inver\t$inver_c\n";
print $out_sum "Indel\t$indel_c\n";
print $out_sum "Indel_0_5\t$indel_c0_5\n";
print $out_sum "Indel_5_20\t$indel_c5_20\n";
print $out_sum "Indel_20+\t$indel_c20_more\n";

#Subroutines:

sub check_trend {
    my @numbers = @_;
    my ($increasing, $decreasing) = (0, 0);

    # Compare consecutive numbers
    for my $i (1 .. $#numbers) {
        if ($numbers[$i] > $numbers[$i - 1]) {
            $increasing++;
        } elsif ($numbers[$i] < $numbers[$i - 1]) {
            $decreasing++;
        }
    }

    # Determine the general trend
    if ($increasing > $decreasing) {
        return "Increasing";
    } elsif ($decreasing > $increasing) {
        return "Decreasing";
    } else {
        return "No clear trend.";
    }
}
