#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with all the SpeciesScoreSummary and the Go folder\n";

#For each species in your Go hash folder work out species name
my %go_key;
my @gos=`ls Go/*`;
foreach my $sp (@gos){
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
    #print "$sp_folder[1] $sp\n";
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
my @hit_genes;
open(my $bedin, "<", $bed)   or die "Could not open $bed\n";
while (my $lineB = <$bedin>){
    chomp $lineB;
    my @splbed=split("\t",$lineB);
    $all_genes{$splbed[3]}="YES";
}




print "We run with species $species\n";


my $outname1="$species\.topSynteny.txt";
open(my $out1, ">", $outname1)   or die "Could not open $outname1\n";
my $outname2="$species\.botSynteny.txt";
open(my $out2, ">", $outname2)   or die "Could not open $outname2\n";
my $outname3="$species\.highScore.txt";
open(my $out3, ">", $outname3)   or die "Could not open $outname3\n";
my $outname4="$species\.lowScore.txt";
open(my $out4, ">", $outname4)   or die "Could not open $outname4\n";
my $outname5="$species\.averhigh.txt";
open(my $out5, ">", $outname5)   or die "Could not open $outname5\n";
my $outname6="$species\.averlow.txt";
open(my $out6, ">", $outname6)   or die "Could not open $outname6\n";
my $outname7="$species\.zeros.txt";
open(my $out7, ">", $outname7)   or die "Could not open $outname7\n";
my $background="$species\.bg.txt";
open(my $out8, ">", $background)   or die "Could not open $background\n";

open(my $filein, "<", $file)   or die "Could not open $file\n";
my $header=<$filein>;
while (my $line = <$filein>){
    chomp $line;
    print "$line\n";
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
    #Print out background geens to file:
    print $out8 "$gene\n";
    push (@hit_genes, $gene);
    my $score=$splitl[2];
    my $count=$splitl[3];
    my $average=$splitl[4];

    my @top;
    if ($count >= 5){
        print "HERE3: $count $score $gene\n";
        push (@top, $gene);
        #print "$gene\n";
    }
    my @bot;
    if ($count <= 2){
        push (@bot, $gene);
    }
    my @high;
    if ($score >= 30){
        push (@high, $gene);
    }
    my @low;
    if ($score <= 5){
        push (@low, $gene);
    }
    my @averhigh;
    if ($average >= 10){
        push (@averhigh, $gene);
    }
    my @averlow;
    if ($average <= 5){
        push (@averlow, $gene);
    }
    
    foreach my $ele1 (@top){
        print $out1 "$ele1\n";
        #print "HERE4 $ele1\n";
    }

    foreach my $ele2 (@bot){
        print $out2 "$ele2\n";
    }

    foreach my $ele3 (@high){
        print $out3 "$ele3\n";
    }

    foreach my $ele4 (@low){
        print $out4 "$ele4\n";
    }

    foreach my $ele5 (@averhigh){
        print $out5 "$ele5\n";
    }

    foreach my $ele6 (@averlow){
        print $out6 "$ele6\n";
    }

}

close $out8;

#My zero calculate:
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

`ChopGO_VTS2.pl -i $species\.topSynteny.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2.pl -i $species\.botSynteny.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2.pl -i $species\.highScore.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2.pl -i $species\.lowScore.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2.pl -i $species\.averhigh.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2.pl -i $species\.averlow.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;
`ChopGO_VTS2.pl -i $species\.zeros.txt --GO_file $go_key{$species} -bg $species\.bg.txt`;

close $out1;
close $out2;
close $out3;
close $out4;
close $out5;
close $out6;
close $out7;
close $filein;
