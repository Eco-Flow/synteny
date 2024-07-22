use warnings;
use strict;


print "Please be in folder with all the SpeciesScoreSummary and the Go folder\n";
die "Please specify a percentage cutoff\n" unless(@ARGV==1);
my $cutoff = $ARGV[0];

#For each species in your Go hash folder work out species name
my %go_key;
my @gos=`ls Go/*`;
foreach my $sp (@gos){
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
}

#Read in the species score summary and work out the species being tested.
my $file=`ls *SpeciesScoreSummary.txt`;
chomp $file;
print "FILE: $file\n";
my @split=split(/\./, $file);
my $species=$split[0];


#import original bed used, to estimate # of zeros.
my $bed=`ls $species\.bed`;
chomp $bed;
my %all_genes;

open(my $bedin, "<", $bed)   or die "Could not open $bed\n";
while (my $lineB = <$bedin>){
    chomp $lineB;
    my @splbed=split("\t",$lineB);
    $all_genes{$splbed[3]}="YES";
}


print "We run with species $species\n";

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


# Open the SpeciesScoreSummary and caluclate the total number of species total.
# And calculate the average distance to break scores, so we can find top/bottom 10% (or user choice cutoff).

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
    my $average=$splitl[4];

    #Check for max:
    if ($count >= $max_number_species){
        $max_number_species=$count;
    }
    $total_genes++;
}
close $filein;

my $top_percentage=$total_genes/$cutoff;

print "Total number of genes: $total_genes\n";
print "What percentage of genes are in the top $cutoff \%: $top_percentage\n";

#Sort the file by specific columns:
open (my $data , '<', $file)|| die "could not open $file:\n$!";
my $head_del=<$data>;
my @array=(<$data>);
my @sorted=sort {(split(/\t/,$a))[2]<=>(split(/\t/,$b))[2]} @array;

#TOP : Print off the x percentage of the distance value:
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

#BOT : Print off the last percentage of the distance value:
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

#Sort by total number of species syntenic to:
my @sorted2=sort {(split(/\t/,$a))[3]<=>(split(/\t/,$b))[3]} @array;

#Print off the x percentage of the distance value:
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
#Print off the last percentage of the distance value:
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

#Print off the x percentage of the distance value:
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
#Print off the last percentage of the distance value:
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

`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.topSynteny.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.botSynteny.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.highScore.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.lowScore.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.averhigh.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.averlow.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2_v12.pl -i $species\.$cutoff\.zeros.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;

close $out1;
close $out2;
close $out3;
close $out4;
close $out5;
close $out6;
close $out7;
close $filein2;
