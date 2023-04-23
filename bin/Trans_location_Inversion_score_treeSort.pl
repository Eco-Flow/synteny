#!/usr/bin/perl
use warnings;
use strict;


print "SYNTENY SCORING SCRIPT\n\nPlease be in folder with all the gff files, the anchors and geneScores run\n\n";


my $outname="Trans_location_version.out.txt";
open(my $OUT, ">", $outname)   or die "Could not open $outname\n";
#Print header of output file- see after 236 where we print this line
print $OUT "Comparison\t#sca_total_sp1\t#sca_total_sp2\t#sca_syntenic_sp1\t#sca_syntenic_sp2\t#genes sp1\t#genes sp2\t#syntenic_blocks\ttranslocation_score\tinversion score\tperc_identity\tGenome_length(sp1)\tGenome_length(sp2)\(bps\)\tLengthTransScore\tLengthInversionScore\tTranslocationScore_byNumberAndLength\tTrans_mimimum_count\tInversion_estimate\n";

my @files=`ls *.gff3 `;

my %gene2scaffold;
my %gene2position;
my %genes_per_species;

my %scaffold_min;
my %scaffold_max;

print "1. Running GFF3 files to get lengths of scaffolds\n\n";

#Loop thru the gff3 files to get the synteny info.
foreach my $file (@files){
    chomp $file;
    #get name of species off gff files (sep by dot)
    my @split_name=split(/\./, $file);
    my $species=$split_name[0];
    #print "$species\n";
    open(my $filein, "<", $file)   or die "Could not open $file\n";
    my $current_min=100000000000;
    my $current_max=0;

    while (my $line=<$filein>){
        chomp $line;
        my @split=split("\t", $line);
        my $scaffold=$split[0];
        my $pos1=$split[3];
        my $pos2=$split[4];
        #print "$line\n";
        if ($line =~ /^#/){

        }
        else{
            if (exists $scaffold_min{$species}{$scaffold}){
                if ($pos1 <= $scaffold_min{$species}{$scaffold}){
                    $scaffold_min{$species}{$scaffold}=$pos1;
                    #print "Smallest scaffold corrdinate is smaller than existing one\t $species $scaffold $pos1\n";
                }
            }
            else{
                #print "Didnt exist $species $scaffold, so create a new hash entry with start\n";
                $scaffold_min{$species}{$scaffold}=$pos1;
            }
            
            if (exists $scaffold_max{$species}{$scaffold}){
                if ($pos2 >= $scaffold_max{$species}{$scaffold}){
                    $scaffold_max{$species}{$scaffold}=$pos2;
                    #print "If pos2 is greater than store value, add to hash for this scaffold\t $species $scaffold $pos2\n";
                }
            }
            else{
                $scaffold_max{$species}{$scaffold}=$pos2;
            }
            
        }
        
        if ($split[2]){
            #Take both the gene and mRNA line, both have gene or transcript IDs we can link to a scaffold
            if ($split[2] eq "mRNA" || $split[2] eq "gene" ){
                my @split9=split("\;",$split[8]);
                my @id_spl=split("\=", $split9[0]);
                my $gene=$id_spl[1];
                #Use transcript id from GFF to make a gene to scaffold hash
                $gene2scaffold{$species}{$gene}=$scaffold;
                #We now create a gene position in scaffold hash here.
                $gene2position{$species}{$gene}=$pos1;
            }
        #print "$split[2]\n";       
            if ($split[2] eq "gene" ){
                $genes_per_species{$species}++;
        #print "$split[2] is gene $species\n";
            }
        }
    
    }
    close $filein;
}



#Now make a hash with all the species scaffold length sums, equiv of genome size,, well from first to last gene on a chromosome (added up for all chromosomes).

print "2. Calculating genome lengths from GFF3s, plus chromosome number\n\n";
my %genome_length;
my %species_chromsome_count;
foreach my $species_loop (keys %scaffold_max){
    my $total=0;
    my $chromo_count=0;
    foreach my $scaffold_loop (keys %{$scaffold_max{$species_loop}}){
        #print "$scaffold_max{$species_loop}{$scaffold_loop}\n";
        #print "$scaffold_min{$species_loop}{$scaffold_loop}\n";
        my $len=$scaffold_max{$species_loop}{$scaffold_loop}-$scaffold_min{$species_loop}{$scaffold_loop};
        $total+=$len;
        $chromo_count++;
    }
    $genome_length{$species_loop}=$total;
    $species_chromsome_count{$species_loop}=$chromo_count;
    print "$species_loop = $total on $chromo_count chromosomes\n";
}   

print "\n3. Compiling percent identity scores\n\n";

my %perc_identity_store;
my @identity=`ls *percent.similarity`;

foreach my $file3 (@identity){
    chomp $file3;
    my @split_name=split(/\./, $file3);
    my $species1=$split_name[0];
    my $species2=$split_name[1];
    my $combName="$species1\.$species2";
    my $perc=`cat $file3`;
    chomp $perc;
    $perc_identity_store{$combName}=$perc;
}


print "We have now saved the general information about the annotations.\n\n";

my @anchors=`ls *.anchors`;

print "4. Now we read in the syntenic anchor files.\n\n";

foreach my $file2 (@anchors){
    chomp $file2;
    #Sort out which name comaparison we are making
    my @split_name=split(/\./, $file2);
    my $species1=$split_name[0];
    my $species2=$split_name[1];
    my $combName="$species1\.$species2";
    open(my $filein2, "<", $file2)   or die "Could not open $file2\n";


    #Set up the stores for each file run:
    my %species1_scaffolds_in_synteny;
    my %species2_scaffolds_in_synteny;
    my $number_of_syntenic_blocks=0;
    #Store the sp2 scaffold names and count later
    my %sp1_to_sp2_scaffoldNames;
    #Store the number of syntenic blocks per scaffold.
    my %syntenic_scaffold_count;
    #temp store of last hit.
    my $last;

    #Remove first #, as it will count as an extra sytnenic block:
    my $removed_first_hash=<$filein2>;

    while (my $line=<$filein2>){
        chomp $line;
        
        if ($line =~ m/^#/){
            #Reset block
            $number_of_syntenic_blocks++;
            #if last exists then we have already processed a synteny block. Now summarise this last block.
            #This section is currently not used (next four lines, but could be useful to know).
            if ($last){
                #count the number of syntenic blocks per scaffold. 
                $syntenic_scaffold_count{$combName}{$last}++;
            }
        }
        else{
            
            my @split=split("\t", $line);
            my $gene1=$split[0];
            my $gene2=$split[1];
            #scaffold of sp1 and sp2
            my $scaf_sp1=$gene2scaffold{$species1}{$gene1};
            my $scaf_sp2=$gene2scaffold{$species2}{$gene2};
            #scaffold count for genes in synteny.
            $species1_scaffolds_in_synteny{$scaf_sp1}++;
            $species2_scaffolds_in_synteny{$scaf_sp2}++;

            #The scaffold to scaffold hits (per gene pair). So we can determine if all hit to just one other chromsome, or there has been a translocation.
            if ($sp1_to_sp2_scaffoldNames{$scaf_sp1}{$scaf_sp2}){
                $sp1_to_sp2_scaffoldNames{$scaf_sp1}{$scaf_sp2}++;
            }
            else{
                $sp1_to_sp2_scaffoldNames{$scaf_sp1}{$scaf_sp2}=1;
            }
            
            #Count of last scaffold name of sp 1 
            $last=$scaf_sp1;

        }

    }

    #Section not currently used, but would rpint the number of syntenic blocks per chromosome.
    foreach my $comparison (keys %syntenic_scaffold_count) {
        foreach my $scaffold (keys %{$syntenic_scaffold_count{$comparison}}) {
            #This line could print the number of syntenic blocks per chromsome there are.
            #print "$comparison\t$scaffold\t$syntenic_scaffold_count{$comparison}{$scaffold}\n";
        }
    }

    #This section gets the number of scaffolds in total with some synteny to another in sp2
    my $total_scaffold_species1=keys %species1_scaffolds_in_synteny;
    my $total_scaffold_species2=keys %species2_scaffolds_in_synteny;

    #This gets the genome lengths
    my $total_length_species1=$genome_length{$species1};
    my $total_length_species2=$genome_length{$species2};

    #My inversion score (OLD); simply the number of syntenic blocks in total, divided by the total number of scaffolds that had synteny (in a syntenic block) or by the total length of the syntenous chromosomes (e.g. from first gen to last gene in each chromosome,, in base pairs[bp, ByLength in code]).

    my $inversion_score=$number_of_syntenic_blocks/$total_scaffold_species1;
    my $inversion_score_By_length=$number_of_syntenic_blocks/$total_length_species1;


    #My translocation score (a mean of percentage chromosome identities to chromsomes in sp2. 1=100% identical in sp2s chromosomes, and 0=0% identity to sp2's chromosomes).
    my %all_translocation_scores;
    my $num_scaff_in_sp1=0;
    my $totalNumScaffolds_sp2=0;
    foreach my $sp1scaffold (keys %sp1_to_sp2_scaffoldNames) {
        #print "SCAFF1 : $sp1scaffold in sp1\n";
        $num_scaff_in_sp1++;
        my $total;
        my $best=0;
        
        foreach my $sp2scaffold (keys %{$sp1_to_sp2_scaffoldNames{$sp1scaffold}}) {
            #print "SCAFF2 : $sp2scaffold in sp2\n";
            #print "$sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold}\n";
            $total+=$sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold};
            if ($sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold} > $best){
                $best=$sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold};
            }
            $totalNumScaffolds_sp2++;
        }
        my $scaffold_tranlocation_score=$best/$total;
        $all_translocation_scores{$combName}+=$scaffold_tranlocation_score;
    }

    #Now minus the total num of scaffold hit to in sp2, by the total number of chromosomes in sp1 (as we should expect to see at least the number of chromosomes in sp1).
    my $translocation_score_By_number=$totalNumScaffolds_sp2-$num_scaff_in_sp1;

    my $translocation_score=$all_translocation_scores{$combName}/$total_scaffold_species1;
    my $translocation_score_By_length=$all_translocation_scores{$combName}/$total_length_species1;
    my $translocation_score_By_number_and_length=$translocation_score_By_number/$total_length_species1;

    #My new translocation score. Simpler, just how many is the minimum number of translocations found.
    my $trans_minimum_found=$totalNumScaffolds_sp2-$num_scaff_in_sp1;

    #My new inversion score based on the total number of syntenic blocks in a chromosome, minus any that are caused by translocations
    my $inversion_estimate=$number_of_syntenic_blocks-$trans_minimum_found;

    my $perc=$perc_identity_store{$combName};
    
    print $OUT "$combName\t$species_chromsome_count{$species1}\t$species_chromsome_count{$species2}\t$total_scaffold_species1\t$total_scaffold_species2\t$genes_per_species{$species1}\t$genes_per_species{$species2}\t$number_of_syntenic_blocks\t$translocation_score\t$inversion_score\t$perc\t$total_length_species1\t$total_length_species2\t$translocation_score_By_length\t$inversion_score_By_length\t$translocation_score_By_number_and_length\t$trans_minimum_found\t$inversion_estimate\n";


    close $filein2;
}


`plotting-inversions-treeSort.R > R_output.txt`
