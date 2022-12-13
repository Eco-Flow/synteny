#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with all the SpeciesScoreSummary and the Go folder\n";


my %go_key;
my @gos=`ls Go/*`;
foreach my $sp (@gos){
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
    #print "$sp_folder[1] $sp\n";
}


my $file=`ls *SpeciesScoreSummary.txt`;
chomp $file;
print "FILE: $file\n";
my @split=split(/\./, $file);
my $species=$split[0];

print "We run with species $species\n";


my $outname1="$species\.top.txt";
open(my $out1, ">", $outname1)   or die "Could not open $outname1\n";
my $outname2="$species\.bot.txt";
open(my $out2, ">", $outname2)   or die "Could not open $outname2\n";
my $outname3="$species\.high.txt";
open(my $out3, ">", $outname3)   or die "Could not open $outname3\n";
my $outname4="$species\.low.txt";
open(my $out4, ">", $outname4)   or die "Could not open $outname4\n";

open(my $filein, "<", $file)   or die "Could not open $file\n";
my $header=<$filein>;
while (my $line = <$filein>){
    chomp $line;
    print "$line\n";
    my @splitl=split("\t", $line);
    my $gene=$splitl[1];
    my $count=$splitl[3];
    my $score=$splitl[2];
    
    #print "HERE: $gene $count $score \n";
    #Rename gene to gene, not transcript with .
    if($gene =~ m/\./){
        my @sp1=split(/\./, $gene);
        $gene=$sp1[0];
        #print "HERE2: $gene $count $score \n";
    }
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
        push (@high, $gene);
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

}

#Now run Chopgo

`ChopGO_VTS.pl -i $species\.top.txt --GO_file $go_key{$species}`;
`ChopGO_VTS.pl -i $species\.bot.txt --GO_file $go_key{$species}`;
`ChopGO_VTS.pl -i $species\.high.txt --GO_file $go_key{$species}`;
`ChopGO_VTS.pl -i $species\.low.txt --GO_file $go_key{$species}`;


close $out1;
close $out2;
close $out3;
close $out4;
close $filein;
