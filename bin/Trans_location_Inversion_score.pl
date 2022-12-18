#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with all the gff files, the anchors  and ... run\n\n";


my $outname="Trans_location_version.out.txt";
open(my $OUT, ">", $outname)   or die "Could not open $outname\n";
print $OUT "Comparison\t#sca_sp1\t#sca_sp2\t#syntenic_blocks\ttranslocation_score\tinversion score\tperc_identity\n";

my @files=`ls *.gff3 `;

my %gene2scaffold;

#Loop thru the files to get the synteny info.
foreach my $file (@files){
    chomp $file;
	my @split_name=split(/\./, $file);
	my $species=$split_name[0];
    open(my $filein, "<", $file)   or die "Could not open $file\n";

    while (my $line=<$filein>){
        chomp $line;
        my @split=split("\t", $line);
        my $scaffold=$split[0];
        if ($split[2]){
            if ($split[2] =~ m/mRNA/g || $split[2] =~ m/gene/g ){
                my @split9=split("\;",$split[8]);
                my @id_spl=split("\=", $split9[0]);
                my $gene=$id_spl[1];
                $gene2scaffold{$species}{$gene}=$scaffold;
            }
        }
    
    }

    close $filein;
}


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


print "We have now save the chromosome information for each gene.\n\n";

my @anchors=`ls *.anchors`;


foreach my $file2 (@anchors){
    chomp $file2;
    my @split_name=split(/\./, $file2);
    my $species1=$split_name[0];
    my $species2=$split_name[1];
    my $combName="$species1\.$species2";
    open(my $filein2, "<", $file2)   or die "Could not open $file2\n";


    #Store for each file run:
    my %species1_scaffolds_in_synteny;
    my %species2_scaffolds_in_synteny;
    my $number_of_syntenic_blocks=0;
    #Store the sp2 scaffold names and count later
    my %sp1_to_sp2_scaffoldNames;
    #Store the number of syntenic blocks per scaffold.
    my %syntenic_scaffold_count;

    #temp store of last hit.
    my $last;

    while (my $line=<$filein2>){
        chomp $line;
        
        if ($line =~ m/^#/){
            #Reset block
            $number_of_syntenic_blocks++;
            #if last exists then we have already processed a synteny block. Now summarise this last one.
            if ($last){
                #count the number of syntenic blocks per scaffold. 
                $syntenic_scaffold_count{$combName}{$last}++;
            }
        }
        else{
            
            my @split=split("\t", $line);
            my $gene1=$split[0];
            my $gene2=$split[1];
            my $scaf_sp1=$gene2scaffold{$species1}{$gene1};
            my $scaf_sp2=$gene2scaffold{$species2}{$gene2};
            #Count of scaffolds in synteny at least once.
            $species1_scaffolds_in_synteny{$scaf_sp1}++;
            $species2_scaffolds_in_synteny{$scaf_sp2}++;

            #The scaffold to scaffold hits . So we can determine if all hit to just one other chromsome, or there has been a translocation.
            if ($sp1_to_sp2_scaffoldNames{$scaf_sp1}{$scaf_sp2}){
                $sp1_to_sp2_scaffoldNames{$scaf_sp1}{$scaf_sp2}++;
            }
            else{
                $sp1_to_sp2_scaffoldNames{$scaf_sp1}{$scaf_sp2}=1;
            }
            

            $last=$species1_scaffolds_in_synteny{$scaf_sp1};

        }

    }


    my $total_scaffold_species1=keys %species1_scaffolds_in_synteny;
    my $total_scaffold_species2=keys %species2_scaffolds_in_synteny;


    #My inversion score; simply the number of syntenic blocks in total, divided by the total number of scaffolds that had synteny (in a syntenic block).

    my $inversion_score=$number_of_syntenic_blocks/$total_scaffold_species1;



    #My translocation score (a mean of percentage chromosome identities to chromsomes in sp2. 1=100% identical in sp2s chromosomes, and 0=0% identity to sp2's chromosomes).
    my %all_translocation_scores;
    foreach my $sp1scaffold (keys %sp1_to_sp2_scaffoldNames) {
        my $total;
        my $best=0;
        foreach my $sp2scaffold (keys %{$sp1_to_sp2_scaffoldNames{$sp1scaffold}}) {
            $total+=$sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold};
            if ($sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold} > $best){
                $best=$sp1_to_sp2_scaffoldNames{$sp1scaffold}{$sp2scaffold};
            }
        }
        my $scaffold_tranlocation_score=$best/$total;
        $all_translocation_scores{$combName}+=$scaffold_tranlocation_score;
    }



    my $perc=$perc_identity_store{$combName};

    my $translocation_score=$all_translocation_scores{$combName}/$total_scaffold_species1;

    #print "$combName\n";
    
    print $OUT "$combName\t$total_scaffold_species1\t$total_scaffold_species2\t$number_of_syntenic_blocks\t$translocation_score\t$inversion_score\t$perc\n";


    close $filein2;
}

