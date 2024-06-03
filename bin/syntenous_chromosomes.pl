#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Basic qw(:all);

# This script sorts out the best order of chromosomes on the plot. 
# It counts up the number of matches of each pair from the anchor files,
# and takes the first chromosome as the one with the most hits and its
# partner (in sp2). If they have already been used in a pair (it is ignored). 

die "Provide the name of the first bed, second bed, anchor (.new) file file\n" if (@ARGV!=3); 

print "Script is running\n";

my $BED_file1 = $ARGV[0];
open(my $BED1, "<", $BED_file1) or die "Could not open $BED_file1\n";
my $BED_file2 = $ARGV[1];
open(my $BED2, "<", $BED_file2) or die "Could not open $BED_file2\n";
my $ANCHORS = $ARGV[2];
open(my $ANCH, "<", $ANCHORS) or die "Could not open $ANCHORS\n";

my $outfile="seqids_karyotype.txt";
open(my $out, "> $outfile") or die "error opening $outfile. $!";

#Variables to store chromosome lengths:
my %sp1_chromo_lengths;
my %sp2_chromo_lengths;

#read in the first bed file
my %SP1;
my %SP1_start_coord;
while (my $line=<$BED1>) {
    chomp $line;
	my @spl=split ("\t", $line);

    #Save gene to chromosome information:
    $SP1{$spl[3]}=$spl[0];
    $SP1_start_coord{$spl[3]}=$spl[1];

    #Calc max length of chromosome:
    if ($sp1_chromo_lengths{$spl[0]}){
        my $old=$sp1_chromo_lengths{$spl[0]};
        if ($spl[2] >= $old){
            $sp1_chromo_lengths{$spl[0]}=$spl[2];
        }
    }
    else{
        $sp1_chromo_lengths{$spl[0]}=$spl[2];
    }
}

#read in the secon bed file
my %SP2;
my %SP2_start_coord;
while (my $line2=<$BED2>) {
    chomp $line2;
	my @spl=split ("\t", $line2);
    
    #Save gene to chromosome information:
    $SP2{$spl[3]}=$spl[0];
    $SP2_start_coord{$spl[3]}=$spl[1];

    #Calc max length of chromosome:
    if ($sp2_chromo_lengths{$spl[0]}){
        my $old=$sp2_chromo_lengths{$spl[0]};
        if ($spl[2] >= $old){
            $sp2_chromo_lengths{$spl[0]}=$spl[2];
        }
    }
    else{
        $sp2_chromo_lengths{$spl[0]}=$spl[2];
    }
}

#Close Bed files
close $BED1;
close $BED2;

#Calc the scaffold name in the anchors file
my %chromo_pairs;

#Establish the last starting point in bed file.
my $last_sp1="NA";
my $last_sp2="NA";

#Chromo pairs direction score.
my %chromo_sp1_scores;
my %chromo_sp2_scores;

while (my $line3=<$ANCH>) {
    chomp $line3;
    my @spl=split ("\t", $line3);
    my $len=scalar(@spl);
    if ( $len == 3 ){
        my $chr1=$SP1{"$spl[0]"};
        my $chr2=$SP2{"$spl[1]"};
        my $join="$chr1\_\@\_$chr2";
        $chromo_pairs{"$join"}++;
        #my current SP1 bed start is:
        my $curr1=$SP1_start_coord{$spl[0]};
        #my current SP2 bed start is:
        my $curr2=$SP2_start_coord{$spl[1]};

        if ($last_sp1 eq "NA"){
            #Do nothing, there is no previous to compare.
            $last_sp1=$curr1;
        }
        else{
            
            #if the next line in synteny block is gt or lt the previous, it gives us a clue to the orientation of the two species chromosomes.
            if ($curr1 >= $last_sp1){
                #we are on the forward strand direction.
                $chromo_sp1_scores{$chr1}++;
            }
            else{
                #we are on the reverse strand direction.
                $chromo_sp1_scores{$chr1}--;
            }
            $last_sp1=$curr1;
        }

        if ($last_sp2 eq "NA"){
            #Do nothing, there is no previous to compare.
            $last_sp2=$curr2;
        }
        else{

            #if the next line in synteny block is gt or lt the previous, it gives us a clue to the orientation of the two species chromosomes.
            if ($curr2 >= $last_sp2){
                #we are on the forward strand direction.
                $chromo_sp2_scores{$chr2}++;
            }
            else{
                #we are on the reverse strand direction.
                $chromo_sp2_scores{$chr2}--;
            }
            $last_sp2=$curr2;
        }
    }
    else{
        #should be a ###, which represents a syntenic break.
        $last_sp1="NA";
        $last_sp2="NA";
    }
}

close $ANCH;


#Try to check the syntenic block order and see if it makes more sense to flip the chromosome. 
my %Orientation_score;

print "Checking if chromosomes are in the correct orientation\n\n";

open(my $ANCH2, "<", $ANCHORS) or die "Could not open $ANCHORS\n";

my $last_chrom1="NA";
my $last_chrom2="NA";
my $n=0;

while (my $line4=<$ANCH2>) {
    chomp $line4;
    my @spl=split ("\t", $line4);
    my $len=scalar(@spl);
    if ($line4 =~ "###"){
        my $first_line_of_block=<$ANCH2>;
        chomp $first_line_of_block;
        my @spl2=split ("\t", $first_line_of_block);
        #print "$spl2[0]\n";
        my $chr_sp1=$SP1{"$spl2[0]"};
        my $chr_sp2=$SP2{"$spl2[1]"};
        my $join="$chr_sp1\_\@\_$chr_sp2";
        my $coord_sp1=$SP1_start_coord{"$spl2[0]"};
        my $coord_sp2=$SP2_start_coord{"$spl2[1]"};

        if ($last_chrom1 eq "NA"){
            #Do nothing, except register that the last block was present. This is after our firts run:
            $last_chrom1=$chr_sp1;
            $last_chrom2=$chr_sp2;
            #print "Should only happen once\n";
        }
        else{
            #If there was a last chromosome, it means it follows the previous block,,, but may not be on the same chromosome.
            if ($last_chrom1 ne $chr_sp1){
                #Means that the next syntenic block is not on same chromosome, so we can ignore this transition:
                #print "These blocks are not on same chromosome, $last_chrom1 ne $chr_sp1\n";
                $last_chrom1=$chr_sp1;
                $last_chrom2=$chr_sp2;
                #print "$n $line4\n$last_chrom2 NE $chr_sp2\n\n";
            }
            else{
                #Means that the last chromosome in sp1 was the same, we are in the same block.
                if ($last_chrom2 eq $chr_sp2){
                    #Means last sp2 chromosome was also on same one as before. 
                    #print "$n $line4\n$last_chrom2 eq $chr_sp2\n\n";

                    #First check the lengths of the chromosomes:
                    my $sp1_chr_len=$sp1_chromo_lengths{$chr_sp1};
                    my $sp2_chr_len=$sp2_chromo_lengths{$chr_sp2};

                    #print "$sp1_chr_len  len   $sp2_chr_len\n";

                    #Then we can check their relative positions on their chromosomes:

                    #print "$coord_sp1 / $sp1_chr_len\n";
                    #print "$coord_sp2 / $sp2_chr_len\n";
                    my $pos_sp1=$coord_sp1/$sp1_chr_len;
                    my $pos_sp2=$coord_sp2/$sp2_chr_len;


                    if ($Orientation_score{$join}){
                        my $old=$Orientation_score{$join};
                        my $joint="$old\n$pos_sp1\t$pos_sp2";
                        $Orientation_score{$join}="$joint";
                    }
                    else{
                        $Orientation_score{$join}="$pos_sp1\t$pos_sp2";
                    }

                    $last_chrom1=$chr_sp1;
                    $last_chrom2=$chr_sp2;
                    $n++;
                }
                else{
                    #print "$n $line4\n$last_chrom2 NE $chr_sp2\n\n";
                    $last_chrom1=$chr_sp1;
                    $last_chrom2=$chr_sp2;
                    $n++;
                }
            }
        }
    }
}

close $ANCH2;

#Now check orientation data and run correlation analysis, to check if a chromosome is in general positively or negatively correlated, which will tell us the best orientation.
my %Chromosomes_in_sp2_to_be_flipped;

foreach my $key ( keys %Orientation_score ) {
    my @pairs=split("\n", $Orientation_score{$key});
    #print "$pairs[0]\n";
    my @sp1_coords;
    my @sp2_coords;

    #print "Correlation ( $key ):\n";

    foreach my $block (@pairs){
        my @spl=split("\t", $block);
        push (@sp1_coords, $spl[0]);
        push (@sp2_coords, $spl[1]);
        #print "$block\n";
    }
    my $result= correlation(\@sp1_coords,\@sp2_coords);
    if ($result <= -0.5){
        my @sp_key=split("\_\@\_", $key);
        #print "$result\n";
        $Chromosomes_in_sp2_to_be_flipped{$sp_key[1]}="TRUE";
    }
}


#Now we make new bed files with some chromosomes reversed where we get a strong negative correlation. Maybe > -0.5


#open(my $BED1_tbc, "<", $BED_file1) or die "Could not open $BED_file1\n";  <We only need to flip sp2.>
open(my $BED2_tbc, "<", $BED_file2) or die "Could not open $BED_file2\n";
my $outfile_bed="$BED_file2\.flipped.bed";
open(my $out_bed, "> $outfile_bed") or die "error opening $outfile_bed. $!";

#Read through species 2 bed file and flip the coordinates of chromosomes in wrong order. 
while (my $line=<$BED2_tbc>) {
    chomp $line;
    my @spl=split ("\t", $line);
    if ($Chromosomes_in_sp2_to_be_flipped{$spl[0]}){
        #print "$spl[0] $spl[3]\n";
        #If chromosome exists in the listto be flipped. 
        my $chromo_len=$sp2_chromo_lengths{$spl[0]};
        #print "$chromo_len\n";
        my $new_start=$chromo_len-$spl[2];   #The new start coordinate of the gene is the chromosome length minus the 2nd corrdinate
        my $new_end=$chromo_len-$spl[1];     #The new end coordinate of the gene is the chromosome length minus the 1st corrdinate
        #print "$line\n";
        #print "$spl[0]\t$new_start\t$new_end\t$spl[3]\t$spl[4]\t$spl[5]\n";
        print $out_bed "$spl[0]\t$new_start\t$new_end\t$spl[3]\t$spl[4]\t$spl[5]\n";
    }
    else{
        print $out_bed "$line\n";
    }
}

print "\n\nStep 2\n\n";

print "These are the pairs that were evaluated in this script: (and number of gene pairs)\n";

my %done_check_sp1;
my %done_check_sp2;
my @final_list_sp1;
my @final_list_sp2;
foreach my $key ( sort {$chromo_pairs{$b} <=> $chromo_pairs{$a}} keys %chromo_pairs ) {

    print "$key $chromo_pairs{$key}\n";
    my @spl=split("\_\@\_", $key);
    if ($done_check_sp1{$spl[0]}){
        #do nothing
    }
    else{
        push (@final_list_sp1, "$spl[0]");
        $done_check_sp1{$spl[0]}="yes"
    }

    if ($done_check_sp2{$spl[1]}){
        #do nothing
    }
    else{
        push (@final_list_sp2, "$spl[1]");
        $done_check_sp2{$spl[1]}="yes"
    }
}

my $join_sp1=join("\,", @final_list_sp1);
my $join_sp2=join("\,", @final_list_sp2);

print $out "$join_sp1\n$join_sp2\n";

sub correl{
        my ($ssxx, $ssyy, $ssxy) = @_;
        my $sign = $ssxy/abs($ssxy);
        my $correl = $sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
        return $correl;
}

