use warnings;
use strict;
use List::Util qw(max);
use List::Util qw(min);


######################################################################################

###A classifier of syntenic breaks###

#This script has been written to provide basic counts of junctions in the genome that 
#appear to be related to inversions or translocations. 
#This is sometimes difficult to ascertain, due to duplications or unexpected gaps that 
#can break synteny, but in fact are not linked to any structural change such as an 
#inversion (a reverse oriented region without a chromosome) or translocation (movement 
#between chromosomes)

######################################################################################


#Loop through all the anchors files in this directory, and work out the files we need to run the whole script.
my @anchors=`ls *.anchors`;
print "Read in the syntenic anchor files.\n\n";
foreach my $file2 (@anchors){
     chomp $file2;
     #Sort out which name comparison we are making
     my @split_name=split(/\./, $file2);
     my $species1=$split_name[0];
     my $species2=$split_name[1];
     my $combName="$species1\.$species2";
     print "$species1 and $species2\n";

     #Provide a bed file of each species, as well as the unfiltered last (pair wise "blast" hits), plus the anchor files, both with lifted hits and without.

     #Please specify (1)bed sp1 , (2)bed sp2, (3)last sp1 to sp2 (not last filtered, we need to consider all hits), (4)anchor sp1 to sp2 lifted.
     my $bed1 = "$species1\.bed";
     my $bed2 = "$species2\.bed";
     my $last = "$combName\.last.filtered";
     my $anc2 = $file2;


     ### Input files into file handles ###


     #Initiate out file handle
     my $outfile="$combName\.Bed_line_by_line_anchored_equivalent.tsv";
     open(my $out, ">", $outfile)   or die "Could not open $outfile \n";

     #Read in GFF/bed file pairs of species 1:
     open(my $in1, "<", $bed1)   or die "Could not open $bed1 \n";

     #Read in GFF/bed file pairs of species 2:
     open(my $in2, "<", $bed2)   or die "Could not open $bed2 \n";

     #Read in last file with 1 to 1 matches to make sure the gene is orthogous between sp1 and sp2.
     open(my $in_last, "<", $last)   or die "Could not open $last \n";

     #Initate gene duplicate store, where a single gene is located in more than one synteny block.
     my $out_doubles="$combName\.Double_gene_entries.txt";
     open(my $out_d, ">", $out_doubles)   or die "Could not open $out_doubles \n";

     #Read out, break junction info
     my $out_break="$combName\.Break_junction_information.txt";
     open(my $outb, ">", $out_break)   or die "Could not open $out_break\n";

     #Store for each syntenic blocks gene order in sp2. Needed for next script.
     my $out_sp2_order="$combName\.Sp2_synteny_order.txt";
     open(my $out_order, ">", $out_sp2_order)   or die "Could not open $out_sp2_order\n";

     #Read in syntenic anchor files (produced by MScanX using the program jcvi) that allows liftover of events
     open(my $in_anc2, "<", $anc2)   or die "Could not open $anc2\n";

     print "Starting script!\n\n";



     #put all of bed into a sorted hash (%loc_store1), 
     #Where chromosome is the first key, then start position which is all equal gene name.
     #And record the gene position on the chromosome (%gene_order_hash).
     #So as we go through each chromosome we have a count ($gene_order++).
     #which is added to the hash (%gene_order_hash) to store the order information
     my %loc_store1;
     my $last_chromo1="STARTING FROM NOTHING";
     my $gene_order=0;
     my %gene_order_hash;
     while ( my $line1 = <$in1> ){
          chomp $line1;
          my @sp=split("\t", $line1);
          if ($loc_store1{$sp[0]}{$sp[1]}){  #If gene start already exists on this chromosome, we have overlapping gene regions. This is probably a weird annotation thing (maybe two isoforms of the same gene?)
               $loc_store1{$sp[0]}{$sp[1]+1}=$sp[3]; #Then add a single base pair later the new gene start position (this allows us to keep this entry, but in future versions we could remove this step, or merge the "genes")
          }
          else{
               $loc_store1{$sp[0]}{$sp[1]}=$sp[3];
          }
          
          $gene_order++;
          $gene_order_hash{$sp[3]}=$gene_order;
          if ($last_chromo1 ne $sp[0]){
               $gene_order=0;
          }
          $last_chromo1=$sp[0];
          #temporary test, to remove next section
          #if ($sp[3] eq "transcript:ENSBTST00005057984"){
          #     print "HEREEE:  $sp[0]  $sp[1]  $sp[3]\n";
          #}
          
     }

     #Then repeat the previous part with 
     #put all of bed into a sorted hash, by chromosome then start position, to equal gene name (%loc_store2).
     #Find out the sp2 species genes to chromosomes (%SP2_gene_to_chromo).
     #And record the gene position on the chromosome (%gene_order_hash).
     my %loc_store2;
     my $last_chromo2="STARTING FROM NOTHING";
     my %SP2_gene_to_chromo;
     my $gene_order2=0;
     my %gene_order_hash2;
     my %gene_order_hash2_rev;
     while ( my $line2 = <$in2> ){
          chomp $line2;
          my @sp=split("\t", $line2);
          $loc_store2{$sp[0]}{$sp[1]}=$sp[3];
          $SP2_gene_to_chromo{$sp[3]}=$sp[0];

          $gene_order2++;
          $gene_order_hash2{$sp[3]}=$gene_order2;
          $gene_order_hash2_rev{$gene_order2}=$sp[3];
          #print "$sp[3] $gene_order2\n";
          if ($last_chromo2 ne $sp[0]){
               $gene_order2=0;
          }
          $last_chromo2=$sp[0];
     }


     #Reading in the last info, to see which genes are "pairs", homogolous (1to1).

     my %last_hits; #gene name SP1 = gene name SP2 (best score)
     my %last_perc_score; # A hash to store the information of the percentage identity matched to top hit of each pair of genes from last.
     my %scores;
     while ( my $line3 = <$in_last> ){
          chomp $line3;
          my @sp=split("\t", $line3);
          if ($line3 =~ m/^#/){
               #Do nothing, this line is just the comments
          }
          else{
               if ($last_hits{$sp[0]}){
                    my $old=$scores{$sp[0]};
                    if ($sp[11] >= $old){
                         $last_hits{$sp[0]}="$sp[1]";
                         $scores{$sp[0]}=$sp[11];
                    }
                    $last_perc_score{$sp[0]}{$sp[1]}=$sp[2];
               }
               else{
                    $last_hits{$sp[0]}=$sp[1];
                    $scores{$sp[0]}=$sp[11];
                    #print "$sp[0] $sp[1]\n";
                    $last_perc_score{$sp[0]}{$sp[1]}=$sp[2];
               }
          }
     }


     #Read in Synteny block file

     ### Which should look like the following (your gene names could look completely different):
     #transcript:ENSBTST00005074599 rna-XM_006571530.3  1770
     #transcript:ENSBTST00005057529 rna-NM_001278335.1  5980

     my %synteny_genes_to_block_lift;       #Save the syntenic block number for each gene
     my %synteny_genes_to_block_lift_sp2;   #and the same in SP2 genes.
     my %equival_genes_to_sp2_lift;
     my %gene_to_synteny_block;
     my $synt_block_lift=0;
     my $gene_double_count=0;
     my $gene_double_count_withinblock=0;

     my %gene_lifted_status; #Record if the gene pairs were lifted or not, these are often poor, and duplicate existing matches.
     while ( my $line5 = <$in_anc2> ){
          chomp $line5;

          my @sp=split("\t", $line5);
          my $gene=$sp[0];
          my $gene2=$sp[1];
          #if there is a lift in paired genes, record this to filter it out later, gene with thrid column with L are lifted, and should not be chosen above non lifted genes. According to jcvi issue 595

          if ($line5 =~ m/\#/){
              $synt_block_lift++;
          }
          else{

               #Check if the gene pair was lifted.
               my $lifted=0;
               if ($sp[2] =~ m/L/g){
                    $lifted=1;
                    $gene_lifted_status{$gene}="Lifted";
               }
               # Removed as I think it needs to be computed later. $gene_to_synteny_block{$gene}=$synt_block_lift;

               #If the genes have been assigned before to a syntenic block:
               if ($synteny_genes_to_block_lift{$sp[0]}){
                    if ($lifted){  #if the current entry was lifted, then no real need to enter it again:
                         #ignore
                    }
                    else{
                         #if current entry is not a liftover, then overwrite the existing entry
                         $synteny_genes_to_block_lift{$sp[0]}=$synt_block_lift;
                    }

                    #Now calculate how many genes may be duplicated in the anchors file (may not be real duplications, just uncertainty of gene in anchor)
                    my $old=$synteny_genes_to_block_lift{$sp[0]};
                    if ($old eq $synt_block_lift){ #then its values is the contents of this has for this gene
                         #do nothing, all good
                         $gene_double_count_withinblock++; #count how many times this happens to priont to screen later. 
                    }
                    else{
                         #print "WARNING 2: the same gene was found in two different synteny blocks, which can confuse the classification\te.g.$sp[0]\n";
                         print $out_d "$sp[0]\n";
                         $gene_double_count++; #in two synteny blocks
                    }
               }
               else{
                    $synteny_genes_to_block_lift{$sp[0]}=$synt_block_lift;
               }


               #if gene has already been registered with a pair
               #print "Sanity\n";
               if ($equival_genes_to_sp2_lift{$sp[0]}){
                    #print "Sanity\n";
                    #This can happen, it means that a gene ( $sp[0] ) is present twice because of a duplication (probably)\n";
                    my $old=$equival_genes_to_sp2_lift{$sp[0]};

                    if ($lifted){ #if current gene is a lift, ignore it
      
                    }
                    else{
                         $equival_genes_to_sp2_lift{$sp[0]}="$sp[1]";
                    }
               }
               else{  #then we can save it straight away, no conflict.
                    $equival_genes_to_sp2_lift{$sp[0]}=$sp[1];
               }    
          }
     }


     #We report to screen how many times we find the same gene in multiple anchor blocks.
     print "NOTICE: the same gene was found in two different synteny blocks, $gene_double_count times, this can confuse the orthology, and is likely due to duplications that are now in two different synteny blocks\n\n";
     print "NOTICE: the same gene was found in the same synteny block at least twice, $gene_double_count_withinblock times, this is likely due to duplications\n\n";

     #Print off the header for the main output table:
     print $out "Chromosome (SP1)\tGene_name_sp1\tSynteny_block_number-if_lifted\tEquivalent gene in SP2 from jcvi or last (if not in anchor file)\tChromosome (SP2)\tGene Order SP1\tGene Order SP2\tLast_hit (if present)\n";

     print "Now run each chromosome through the classifier\n\nChromosome\tlength\n";

     #Create hash store of for each chromosome in SP1, for the gene position numbers in species 2. 
     #This will help later with determining if a new syntenic block is formed from an inversion or a huge indel (potentially).
     #As we can check if the next set of genes are present in the next syntenic anchor in the same direction, or are missing in the next anchor group.
     my %Chr_to_gene_position_sp2;

     #Store the syntenic block number of all the chromosomes and their respective genes for species 2.
     my %Chr_gene_to_syntenic_block;

     #Store the Chromosome, to sytenic block range numbers.
     my %Chr_syntenic;

     #now loop through the bed files in chromosome order and start testing back to sp2 and last and positions
     foreach my $chromosomes (sort keys %loc_store1 ){ #For each chromosome in SP1
          #print to screen which chromosomes are being processed:
          print "$chromosomes\t";
          my @all_first;
          foreach my $start_pos (keys %{$loc_store1{$chromosomes}}){
               push (@all_first, $start_pos);
          }
          my $max=max(@all_first);
          print "$max\n"; #print the max of the start positions in this chromosome (is a pseudo chromosome length (in bp), based on length from bp0 to last gene start position).
          
          for (my $i=1; $i<$max+1; $i++){                                  #Now we go through base pair by base pair through the chromosome length.

               if ($loc_store1{$chromosomes}{$i}){                         #If we find a start position at a bp locator within that chromosome, we have located a gene (from the bed file, now in $loc_store1 from SP1)
                    
                    my $gene=$loc_store1{$chromosomes}{$i};                #My gene name is called here, in order through genome

                    print $out "$chromosomes\t$gene";                      #We print off the chromosome name followed by the gene name
                    my $synteny_block_lift;
                    if ($synteny_genes_to_block_lift{$gene}){
                         $synteny_block_lift=$synteny_genes_to_block_lift{$gene};
                         print $out "\t$synteny_block_lift";
                         #$Chr_gene_to_syntenic_block{$chromosomes}{$gene}=$synteny_block_lift;
                    }
                    else{
                         print $out "\tNOT";
                         #$Chr_gene_to_syntenic_block{$chromosomes}{$gene}="NOT";
                    }

                    if ($equival_genes_to_sp2_lift{$gene}){              #If there is an equivalent gene in SP2
                         my $gene_in_sp2=$equival_genes_to_sp2_lift{$gene};
                         if ($SP2_gene_to_chromo{$gene_in_sp2}){         #If the gene in SP2 has an associated chromosome.
                              #We have one hit for the gene, proceed.
                              my $chro_in_sp2=$SP2_gene_to_chromo{$gene_in_sp2};   #Then we store chromosome of sp2 gene in variable.
                              if ($equival_genes_to_sp2_lift{$gene}){             #If if also present in the lifted over coordinates.
                                   my $gene_in_sp2_lift=$equival_genes_to_sp2_lift{$gene};
                                   if ($gene_order_hash2{$gene_in_sp2_lift}){     #If the gene is SP2 has a gene order number (lifted)

                                        #Convert the position to a gene name in sp2
                                        my $gen_pos_sp2=$gene_order_hash2{$gene_in_sp2_lift};

                                        #store sytenic break sp 2 gene positions in a hash;
                                        if ($Chr_syntenic{$chromosomes}{$synteny_block_lift}){
                                             my $old_h=$Chr_syntenic{$chromosomes}{$synteny_block_lift};
                                             $Chr_syntenic{$chromosomes}{$synteny_block_lift}="$old_h\t$gen_pos_sp2";
                                        }
                                        else{
                                             $Chr_syntenic{$chromosomes}{$synteny_block_lift}=$gen_pos_sp2;
                                        }
                                        
                                        print $out "\t$gene_in_sp2_lift\t$chro_in_sp2\t$gene_order_hash{$gene}\t$gen_pos_sp2\t";
                                        #We print this to check with single genes (not multiple), which gene positions in sp2 exist, so we can asses later if there are missing genes in a syntenic break or not. 
                                        $Chr_to_gene_position_sp2{$chromosomes}{$gen_pos_sp2}="yes";

                                        

                                        #my $sp2_gen=$gene_order_hash2_rev{};

                                        #no we have the gene name in sp2:    $gene_in_sp2

                                        $Chr_gene_to_syntenic_block{$chromosomes}{$i}= $gene_order_hash2{$gene_in_sp2};
                                   }
                                   else{
                                        #not sure why this statement works.. it is the same print as above, but needed, as both have content. Whether the statement is true or false.
                                        print $out "\t$gene_in_sp2\t$chro_in_sp2\t$gene_order_hash{$gene}\tNo_gene_pos_sp2\t";
                                        $Chr_to_gene_position_sp2{$chromosomes}{$gene_order_hash2{$gene_in_sp2}}="yes";
                                   }
                              }
                              else{
                                   print $out "\t$gene_in_sp2\t$chro_in_sp2\t$gene_order_hash{$gene}\tNo_Sp2_position (1)\t";
                              }
                         }
                         else{
                              #We have multiple hits for this gene
                              print $out "\t";
                              my @split_genes=split("\t", $gene_in_sp2);
                              foreach my $gen (@split_genes){
                                   print $out "$gen\,";

                              }
                              print $out "\t";



                              my @chromo_hits;
                              foreach my $gen2 (@split_genes){
                                   if ($SP2_gene_to_chromo{$gen2}){
                                        push (@chromo_hits, $SP2_gene_to_chromo{$gen2});
                                        print $out "$SP2_gene_to_chromo{$gen2}\,";
                                   }
                                   else{
                                        print $out "No_gene\,";
                                   }
                              }

                              if ($gene_order_hash2{$equival_genes_to_sp2_lift{$gene}}){
                                   print $out "\t$gene_order_hash{$gene}\t$gene_order_hash2{$equival_genes_to_sp2_lift{$gene}}\t";
                              }
                              else{
                                   print $out "\t$gene_order_hash{$gene}\tNo_gene_order_sp2\t";
                              }
                              
                         }
                    }
                    else{
                         #if there is not an equivalent gene in sp2
                         if ($equival_genes_to_sp2_lift{$gene}){

                              my $gene_in_sp2=$equival_genes_to_sp2_lift{$gene};
                              if ($SP2_gene_to_chromo{$gene_in_sp2}){
                                   #We have one hit for the gene, proceed, print gene name and chromosome in sp2.
                                   my $chro_in_sp2=$SP2_gene_to_chromo{$gene_in_sp2};
                                   print $out "\t$gene_in_sp2\t$chro_in_sp2";
                                   if ($equival_genes_to_sp2_lift{$gene}){
                                        if ($gene_order_hash2{$equival_genes_to_sp2_lift{$gene}}){
                                             print $out "\t$gene_order_hash{$gene}\t$gene_order_hash2{$equival_genes_to_sp2_lift{$gene}}\t";
                                        }
                                        else{
                                             #shouldnt happen, means theres no equivalent gene in sp2 base on gene in the order hash.
                                             print $out "\t$gene_order_hash{$gene}\tNo_Sp2_position (2)  $gene_order_hash2{$gene_in_sp2}\t";
                                        } 
                                   }
                                   else{
                                        #shouldnt happen, means theres no equivalent gene in sp2 lift.
                                        print $out "\t$gene_order_hash{$gene}\tNo_Sp2_position (3) $gene_order_hash2{$gene_in_sp2}\t";
                                   }
                              }
                              else{
                                   #If we have multiple hits for this gene
                                   print $out "\t";
                                   my @split_genes=split("\t", $gene_in_sp2);
                                   foreach my $gen (@split_genes){
                                        print $out "$gen\,";

                                   }
                                   print $out "\t";

                                   my @chromo_hits;  #add all the chromosome locations of the hits
                                   foreach my $gen2 (@split_genes){
                                        push (@chromo_hits, $SP2_gene_to_chromo{$gen2});
                                   }

                                   #check if all the genes chromosome hits are the same, if they are print it just once. Else print all the chromosomes in order of the genes in previous column.
                                   if (all_the_same(@chromo_hits)){
                                        print $out "$chromo_hits[0]";
                                   }
                                   else{
                                        foreach my $entries (@chromo_hits){
                                             print $out "$entries\,";
                                        }
                                   }

                                   my $gene_in_sp2_lift=$equival_genes_to_sp2_lift{$gene_in_sp2};

                                   if ($gene_in_sp2_lift){
                                        #print $out "--> $gene_in_sp2_lift\n";
                                        print $out "\t$gene_order_hash{$gene}\tNo_Sp2_position (4)\t";
                                   }
                                   else{
                                        print $out "\t$gene_order_hash{$gene}\tNo_Sp2_position (4)\t";
                                   }
                                   
                              }
                            
                         }
                         else{   #If we don't have a hit on the anchor or lifted anchor file, check the last and print as much info as possible to try to understand what these genes may be.
                              
                              if ($last_hits{$gene}){
                                   my $last_hit_sp2=$last_hits{$gene};
                                   if ($equival_genes_to_sp2_lift{$gene}){
                                        print $out "\t$last_hit_sp2\t$SP2_gene_to_chromo{$last_hit_sp2}\t$gene_order_hash{$gene}\t$gene_order_hash2{$equival_genes_to_sp2_lift{$gene}}\t";
                                   }
                                   else{
                                        if ($gene_order_hash2{$last_hit_sp2}){
                                             print $out "\t$last_hit_sp2\t$SP2_gene_to_chromo{$last_hit_sp2}\t$gene_order_hash{$gene}\tBlast_hit_position: $gene_order_hash2{$last_hit_sp2}\t";
                                        }
                                        else{
                                             #print "$gene  \t  tooo $last_hit_sp2  \t  meee $gene_order_hash{$gene}  \t here $SP2_gene_to_chromo{$last_hit_sp2}\n";
                                             if ($SP2_gene_to_chromo{$last_hit_sp2}){
                                                  print $out "\t$last_hit_sp2\t$SP2_gene_to_chromo{$last_hit_sp2}\t$gene_order_hash{$gene}\tNo_sp2_position (5)\t";
                                             }
                                             else{
                                                  print $out "\t$last_hit_sp2\tNo_sp2_chromo\t$gene_order_hash{$gene}\tNo_sp2_position (5)\t";
                                             }
                                        }
                                        
                                   }
                                   
                              }
                              else{
                                   #These are when no homolog can be found, by jcvi or by last. 
                                   print $out "\tNo_gene_in_sp2\tNo_equiqalent_sp2\t$gene_order_hash{$gene}\tNo_sp2_position (6)\t";
                              }
                         }
                    }

                    #print out the last information, in case helpful to track specific genes.
                    if ($last_hits{$gene}){
                         print $out "$last_hits{$gene}\t";
                    }
                    else{
                         print $out "NO_BLAST\t";
                    }

                    if ($last_hits{$gene}){
                         my $sp2_last_hit=$last_hits{$gene};
                         if ($last_perc_score{$gene}{$sp2_last_hit}){
                              print $out "$last_perc_score{$gene}{$sp2_last_hit}\n";
                         }
                         else{
                              print $out "NO_PERC\n";
                              #print "\tNot this guy: $gene with $sp2_last_hit\n";
                         }
                         
                    }
                    else{
                         print $out "NO_PERC\n";
                    }
               }
          }
     }




     #Close the previous files
     close $out;
     close $in1;
     close $in2;
     close $in_last;
     close $out_d;
     close $in_anc2;



     #Next we can go through the output of the previous step to calculate the numbers of junctions that are likely inversion, those that are translocations and those which are difficult to grade.
     print "\n\nNext: Lets calculate the numbers of junction events\n\n";

     my $outfile2="$combName\.Counts_of_syntenic_junctions.tsv";
     open(my $out2, ">", $outfile2)   or die "Could not open $outfile2 \n";

     my $infile="$combName\.Bed_line_by_line_anchored_equivalent.tsv";
     open(my $in3, "<", $infile)   or die "Could not open $infile \n";

     #My set up the starting chr, block numbers for sp1 and sp2.
     my $previous_chr="start";
     my $previous_block="start";
     my $previous_sp2_chr="start";
     my $previous_sp2_pos="start";

     #My current and previous gene, so we can get coordinates of events
     my $previous_gene;
     my $current_gene;

     #set up counts of scores:
     my $inversion_count=0;
     my $translocation_count=0;
     my $indel_count=0;
     #my $other_count=0;

     #set up a hash to store the order list for each syntenic block
     my %sp2_syntenic_order;

     print $outb "Chromosome_sp1\tSytenic_block_first\tSytenic_block_second\tFirst_start\tFirst_end\tSecond_start\tSecond_end\tTYPE\tFirst_junction_gene\tSecond_juntion_gene\tFirst_junction_gene(sp2)\tSecond_juntion_gene(sp2)\n";

     while ( my $line3 = <$in3> ){
          chomp $line3;
          my @sp=split("\t", $line3);
          my $chromosome=$sp[0];
          my $synteny_block=$sp[2];
          my $sp2_chromosome=$sp[4];
          my $sp2_position=$sp[6];
          $current_gene=$sp[1];
          

          #if the line (gene) is not in a syntenic block (col sp[5]), then skip the line:
          if ($synteny_block eq "NOT"){
               #Skip this line. As it does help understand the junction change from syntenic block A-B, or B-C
               #it is more likely caused by a duplication, or a deletion in sp2 of this particular gene. 
          }
          else{

               #if store for sp2 order exists for this chromosome and block, add to existing, else make a new hash store.
               if ($sp2_syntenic_order{$chromosome}{$synteny_block}){
                    my $old_syn=$sp2_syntenic_order{$chromosome}{$synteny_block};
                    $sp2_syntenic_order{$chromosome}{$synteny_block}="$old_syn\,$sp2_position";
               }
               else{
                    $sp2_syntenic_order{$chromosome}{$synteny_block}=$sp2_position;
               }

               #continue , as the line is processable.:
               #If previous gene is on the same chromosome(eq) 
               #We also remove any previous and current positions that don't have a gene positions (NOT), as this makes it hard to check which genes are missing in the syntenic break. 
               #This is a syntenic break within a chromosome:
               if ($previous_chr eq $chromosome){
                    #if the previous synteny block number is different(ne).
                    if ($previous_block ne $synteny_block  &&  $previous_chr ne "start"){
                         #Add in syntenic range from previous script.
                         my @previous_range=split("\t", $Chr_syntenic{$chromosome}{$previous_block});
                         my @current_range =split("\t", $Chr_syntenic{$chromosome}{$synteny_block});

                         my $previous_start=shift(@previous_range);
                         my $previous_stop=pop(@previous_range);
                         my $current_start=shift(@current_range);
                         my $current_stop=pop(@current_range);
                         
                         if ($sp2_position ne "No_Sp2_position (4)"  && $previous_sp2_pos ne "No_Sp2_position (4)"){
                              #Sp2 position is present.
                              if ($previous_sp2_chr eq $sp2_chromosome){  

                                   #Add genes in sp2 for output table:
                                   my $previous_gene_equivaent_sp2=$equival_genes_to_sp2_lift{$previous_gene};
                                   my $current_gene_equivaent_sp2=$equival_genes_to_sp2_lift{$current_gene};

                                   if (defined $current_stop){
                                        if (defined $previous_stop){
                                             print $outb "$chromosome\t$previous_block\t$synteny_block\t$previous_start\t$previous_stop\t$current_start\t$current_stop\tINVER\t$previous_gene\t$current_gene\t$previous_gene_equivaent_sp2\t$current_gene_equivaent_sp2\n";
                                        }
                                        else{
                                             print $outb "$chromosome\t$previous_block\t$synteny_block\t$previous_start\tNA\t$current_start\t$current_stop\tINVER\t$previous_gene\t$current_gene\t$previous_gene_equivaent_sp2\t$current_gene_equivaent_sp2\n";
                                        }
                                   }
                                   else{
                                        if (defined $previous_stop){
                                             print $outb "$chromosome\t$previous_block\t$synteny_block\t$previous_start\t$previous_stop\t$current_start\tNA\tINVER\t$previous_gene\t$current_gene\t$previous_gene_equivaent_sp2\t$current_gene_equivaent_sp2\n";
                                        }
                                   }

                                   my $max_pos=max($sp2_position, $previous_sp2_pos);
                                   my $min_pos=min($sp2_position, $previous_sp2_pos);
                                   my @positions_missing;
                                   for (my $i=$min_pos; $i<=$max_pos; $i++){
                                        #print "here $i\n";
                                        push (@positions_missing, $i);
                                   }
                                   my $joint_pos=join("\t", @positions_missing);
                                   #print "Joint $joint_pos";
                                   my $does_gap_gene_exist=0;
                                   #For each my missing 
                                   foreach my $missing (@positions_missing){  # For each missing gene between two syntenic blocks, assess what they are.
                                        if ($Chr_to_gene_position_sp2{$chromosome}{$missing}){
                                             #This position was found before, but now we need to check if it existing in another syntenic block, different to that of the current and previous synetic block numbers.
                                             #If it exists in another block, then it seems more like an indel than a inversion, but if there are no missing genes in between assigned to another block, then we assume its an inversion.
                                             #print "So the position exists\n"; #This means it may not be a simple inversion, it could be a deletion or insertion.  
                                             $does_gap_gene_exist=1;

                                             #my genes in sp2 to syntenic block   ===   $equival_genes_to_sp2_lift
                                             #my $gap_gene_syntenic_block_number=$Chr_gene_to_syntenic_block{}{};

                                        }
                                   }
                                   if ($does_gap_gene_exist){
                                        $indel_count++;
                                   }
                                   else{
                                        $inversion_count++;
                                   }

                              }
                              else{
                                   #Then we have a translocation (seemingly), as within the same chromosome it maps to two chromosomes in species 2.
                                   #e.g. ($previous_sp2_chr ne $sp2_chromosome)
                                   $translocation_count++;

                                   #Add genes in sp2 for output table:
                                   my $previous_gene_equivaent_sp2=$equival_genes_to_sp2_lift{$previous_gene};
                                   my $current_gene_equivaent_sp2=$equival_genes_to_sp2_lift{$current_gene};
                                   
                                   if (defined $current_stop){
                                        if (defined $previous_stop){
                                             print $outb "$chromosome\t$previous_block\t$synteny_block\t$previous_start\t$previous_stop\t$current_start\t$current_stop\tTRANS\t$previous_gene\t$current_gene\t$previous_gene_equivaent_sp2\t$current_gene_equivaent_sp2\n";
                                        }
                                        else{
                                             print $outb "$chromosome\t$previous_block\t$synteny_block\t$previous_start\tNA\t$current_start\t$current_stop\tTRANS\t$previous_gene\t$current_gene\t$previous_gene_equivaent_sp2\t$current_gene_equivaent_sp2\n";
                                        }
                                   }
                                   else{
                                        if (defined $previous_stop){
                                             print $outb "$chromosome\t$previous_block\t$synteny_block\t$previous_start\t$previous_stop\t$current_start\tNA\tTRANS\t$previous_gene\t$current_gene\t$previous_gene_equivaent_sp2\t$current_gene_equivaent_sp2\n";
                                        }
                                   }
                              };

                              #Finally set the current values to the previous, for next round.
                              $previous_chr=$chromosome;
                              $previous_block=$synteny_block;
                              $previous_sp2_chr=$sp2_chromosome;
                              $previous_sp2_pos=$sp2_position;
                              $previous_gene=$sp[1];
                         }
                         else{
                              #Means that the species two gene in the pair did not have a location in the anchor file.
                              #These we ignore as we cannot calculate if they are inversions without the gene position information.
                              #print "Does this happen?\n";
                         }
                    }
                    else{
                         #We can ignore this, this happens when the last syntenic block was the same as the last
                         #ie, there is not a syntenic break within a chromosome.
                    }
               }
               else{
                    #Starting on a new chromosome. So we need to set the previous lines information, but we do not count changes between chromosomes as syntenic breaks (they are chromosomal breaks, or breaks in the annotation)
                    $previous_chr=$chromosome;
                    $previous_block=$synteny_block;
                    $previous_sp2_chr=$sp2_chromosome;
                    $previous_sp2_pos=$sp2_position;
                    $previous_gene=$sp[1];
               }   
          } 
     }

     close $out2;
     close $in3;


     #Print to file all the sp2 order numbers into a file
     foreach my $chro (sort keys %sp2_syntenic_order){
          foreach my $syn_blocks (keys %{$sp2_syntenic_order{$chro}}){
               #print "$chro\t$syn_blocks\t$sp2_syntenic_order{$chro}{$syn_blocks}\n";
               print $out_order "$chro\t$syn_blocks\t$sp2_syntenic_order{$chro}{$syn_blocks}\n";
          }
     }


     #Print a help overview message to screen.
     print "We found\n
     $translocation_count Translocation junctions
     $indel_count Indel or inversion junctions
     ";

}



#Sub routine to check all the elements within a hash are the same
sub all_the_same {
    my @test = @_;
    my %seen;
    foreach my $item (@test){
      $seen{$item}++
    }
    my $size = keys %seen;
    if ($size == 1){
        return 1;
    }
    else{
        return 0;
    }
}

