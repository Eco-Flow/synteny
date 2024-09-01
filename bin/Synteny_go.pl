use warnings;
use strict;

print "Please be in folder with all the SpeciesScoreSummary and the Go folder\n";
die "Please specify a percentage cutoff (e.g. 10)\n" unless(@ARGV==1);
my $cutoff = $ARGV[0];

#For each species in your Go folder work out species name, then store in a hash with species name -> path to file
my @gos=`ls Go/*`;
my %go_key;
foreach my $sp (@gos){
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
}

#Read in the species score summary file and work out the species being tested (should just be one per folder).
my $file=`ls *SpeciesScoreSummary.txt`;
chomp $file;
my @split=split(/\./, $file);
my $species=$split[0];

#import original bed used, which we can use to estimate # of zeros, ones that do not show up on the synteny result, even though the gene exists.
# Store in a hash to specify which genes are present.
my $bed=`ls $species\.bed`;
chomp $bed;
my %all_genes;
open(my $bedin, "<", $bed)   or die "Could not open $bed\n";
while (my $lineB = <$bedin>){
    chomp $lineB;
    my @splbed=split("\t",$lineB);
    $all_genes{$splbed[3]}="YES";
}

#print "We run with species $species\n";

#Specify the output file names, and initialise these handles.
my $outname1="$species\.$cutoff\.topSynteny.txt";
open(my $out1, ">", $outname1)   or die "Could not open $outname1\n";
my $outname2="$species\.$cutoff\.botSynteny.txt";
open(my $out2, ">", $outname2)   or die "Could not open $outname2\n";
my $outname3="$species\.$cutoff\.highScore.txt";
open(my $out3, ">", $outname3)   or die "Could not open $outname3\n";
my $outname4="$species\.$cutoff\.lowScore.txt";
open(my $out4, ">", $outname4)   or die "Could not open $outname4\n";
my $outname5="$species\.$cutoff\.averhigh.txt";
open(my $out5, ">", $outname5)   or die "Could not open $outname5\n";
my $outname6="$species\.$cutoff\.averlow.txt";
open(my $out6, ">", $outname6)   or die "Could not open $outname6\n";
my $outname7="$species\.$cutoff\.zeros.txt";
open(my $out7, ">", $outname7)   or die "Could not open $outname7\n";
my $background="$species\.$cutoff\.bg.txt";
open(my $out8, ">", $background)   or die "Could not open $background\n";
my $outname9="$species\.$cutoff\.top_orthologous.txt";
open(my $out9, ">", $outname9)   or die "Could not open $outname9\n";
my $outname10="$species\.$cutoff\.bot_orthologous.txt";
open(my $out10, ">", $outname10)   or die "Could not open $outname10\n";

# Open the SpeciesScoreSummary (first time) and calculate the total number of species total.
# And calculate the total number of genes, so we can work out the number of genes that equate to the top/bottom 10% (or user choice cutoff).
open(my $filein, "<", $file)   or die "Could not open $file\n";
my $header=<$filein>;
my $max_number_species=0;
my $total_genes=0;
while (my $line = <$filein>){
    chomp $line;
    my @splitl=split("\t", $line);
    my $gene=$splitl[1];
    my $score=$splitl[2];
    my $count=$splitl[3];
    #Check for max:
    if ($count >= $max_number_species){
        $max_number_species=$count;
    }
    $total_genes++;
}
close $filein;
#Here is the top number of genes in the $cutoff percentage:
my $top_percentage=$total_genes/$cutoff;

print "Total number of genes: $total_genes\n";
print "What percentage of genes are in the top $cutoff \%: $top_percentage\n";

#Read in the SpeciesSummary again, remove the header, and sort the file by total_score (all scores for distance from break)):
open (my $data , '<', $file)|| die "could not open $file:\n$!";
my $head_del=<$data>;
my @array=(<$data>);
my @sorted=sort {(split(/\t/,$a))[2]<=>(split(/\t/,$b))[2]} @array;

#HIGH SCORE : Print off the x percentage of the distance value, genes most distant from breaks:
for (my $i=0; $i<=$top_percentage; $i++){
    #print "$sorted[$i]";
    my @spl=split("\t", $sorted[$i]);
    my $gene=$spl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    print $out3 "$gene\n";
};

#LOW SCORE : Print off the last percentage of the distance value, genes closest to breaks:
my $top_minus_cut=$total_genes-$top_percentage;
for (my $i=$total_genes-1; $i>=$top_minus_cut; $i--){
    #print "$sorted[$i]";
    my @spl=split("\t", $sorted[$i]);
    my $gene=$spl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    print $out4 "$gene\n";
};


# Sort by total number of species syntenic to target gene. 
# E.g. if there are 5 total species, a score of 4 means it was in a syntenic block to all other species under study.:
my @sorted2=sort {(split(/\t/,$a))[3]<=>(split(/\t/,$b))[3]} @array;

# TOP SYNTENY: Print off the x percentage of the genes with highest numbers of species syntenic:
for (my $i=0; $i<=$top_percentage; $i++){
    #print "$sorted2[$i]";
    my @spl=split("\t", $sorted2[$i]);
    my $gene=$spl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    print $out1 "$gene\n";    
};

# BOT SYNTENY: Print off the x percentage of the genes with the lowest number of species syntenic:
for (my $i=$total_genes-1; $i>=$top_minus_cut; $i--){
    #print "$sorted2[$i]";
    my @spl=split("\t", $sorted2[$i]);
    my $gene=$spl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    print $out2 "$gene\n";
};

#Sort by average distance to break :
my @sorted3=sort {(split(/\t/,$a))[4]<=>(split(/\t/,$b))[4]} @array;

#AVER_HIGH: Print off the x percentage of the distance value:
for (my $i=0; $i<=$top_percentage; $i++){
    #print "$sorted3[$i]";
    my @spl=split("\t", $sorted3[$i]);
    my $gene=$spl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    print $out5 "$gene\n";
};

#AVER_LOW: Print off the last percentage of the distance value:
for (my $i=$total_genes-1; $i>=$top_minus_cut; $i--){
    #print "$sorted3[$i]";
    my @spl=split("\t", $sorted3[$i]);
    my $gene=$spl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    print $out6 "$gene\n";
};


# Sort by number of species the gene is syntenic to.
my @sorted4=sort {(split(/\t/,$a))[3]<=>(split(/\t/,$b))[3]} @array;

# TOP ORTHOLOGOUS SYNTENY: Print off the x percentage of the genes with highest numbers of species syntenic:
for (my $i=0; $i<=$top_percentage; $i++){
    #print "$sorted2[$i]";
    my @spl=split("\t", $sorted2[$i]);
    my $gene=$spl[1];
    my $ortholougous=$spl[5];  #This line says how many species it was orthologous to
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    if ($max_number_species == $ortholougous){
        print $out9 "$gene\n";
    }
    else{
        #Do not print it out, as the number of total comaprable species, was not the same as the ortholgous number of species. 
    }    
};

# BOT ORTHOLOGOUS SYNTENY: Print off the x percentage of the genes with the lowest number of species syntenic:
for (my $i=$total_genes-1; $i>=$top_minus_cut; $i--){
    #print "$sorted2[$i]";
    my @spl=split("\t", $sorted2[$i]);
    my $gene=$spl[1];
    my $ortholougous=$spl[5];  #This line says how many species it was orthologous to
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }
    if ($max_number_species == $ortholougous){
        print $out10 "$gene\n";
    }
    else{
        #Do not print it out, as the number of total comaprable species, was not the same as the ortholgous number of species. 
    }
};

# Go through the Hits from target species to all other species
# Ending SpeciesScoreSummary.txt 
# 1: Species name (target), 2: gene name, 3: count of distance to syntenic break, 4: total species also within a syntenic block, 5: average score (distance to break).

my @hit_genes;

open(my $filein2, "<", $file)   or die "Could not open $file\n";
my $header2=<$filein2>;
while (my $line = <$filein2>){
    chomp $line;
    #print "$line\n";
    my @splitl=split("\t", $line);
    my $gene=$splitl[1];
    
    #Rename weird transcrip ids with :, usually transcipt:ENSGMT0000012, we want just ENSGMT0000012
    if($gene =~ m/\:/){
        my @sp1=split(/\:/, $gene);
        $gene=$sp1[1];
    }
    # Rename weird NCBI id rna- prefix
    if($gene =~ m/rna-/){
        my @sp1=split(/\-/, $gene);
        $gene=$sp1[1];
    }

    #Print out background genes to file:
    print $out8 "$gene\n";
    push (@hit_genes, $gene);
    my $score=$splitl[2];
    my $count=$splitl[3];
    my $average=$splitl[4];
}

close $out8;

#My zero calculate. 
#Lot of genes in each species could be unique, so they would end up with a score of zero for being in or out of synteny.
#This section stores all these genes in a file:
my $score=0;
my $zscore=0;
my $zhit;

foreach my $geneall (keys %all_genes){
	$zhit=0;
	foreach my $hit_gen (@hit_genes){
		if ($hit_gen eq $geneall){
			$zhit=1;
			$score++;
		}
	}
	if ($zhit){
		#DONT PRINT, as this gene matched
		#print $out7 "$geneall\n";
	}
	else{
		$zscore++;
		print $out7 "$geneall\n";
	}
}
print "$score matches and $zscore zeros\n";


#Now run Chopgo
print "Now run ChopGO : e.g. : ChopGO_VTS.pl -i $species\.top.txt --GO_file $go_key{$species}\n";

`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.topSynteny.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.botSynteny.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.highScore.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.lowScore.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.averhigh.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.averlow.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.zeros.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.top_orthologous.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.bot_orthologous.txt --GO_file $go_key{$species} -bg $species\.$cutoff\.bg.txt`;

close $out1;
close $out2;
close $out3;
close $out4;
close $out5;
close $out6;
close $out7;
close $filein2;
