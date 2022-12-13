#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with all the geneScore.tsv files\n";

my @files=`ls *.geneScore.tsv `;

my %score;
my %count;

#Loop thru the files to get the synteny info.
foreach my $file (@files){
    chomp $file;
	my @split_name=split(/\./, $file);
	my $sp1=$split_name[0];
	my $sp2=$split_name[1];

    
	#Now per file:
    open(my $filein, "<", $file)   or die "Could not open $file\n";

    while (my $line=<$filein>){
        chomp $line;

        my @split=split("\t", $line);
        my $sp_name=$split[0];
        my $gene=$split[1];
        my $val=$split[2];

        #Make sure the first line is the right species, sometimes it prints a previous species (to be fixed)
        if ($sp_name eq $sp1){

            #Calculate total scores
            if ($score{$sp1}{$gene}){
                my $old=$score{$sp1}{$gene};
                $score{$sp1}{$gene}=$val+$old;
            }
            else{
                $score{$sp1}{$gene}=$val;
            }
            #Calculate total count
            if ($count{$sp1}{$gene}){
                $count{$sp1}{$gene}++;
            }
            else{
                $count{$sp1}{$gene}=1;
            }
            
        }
    }
    close $filein;
}



foreach my $species (keys %score) {
    my $outname="$species\_geneScoreSummary.txt";
    open(my $outhandle, ">", $outname)   or die "Could not open $outname\n";
    print $outhandle "Species\tgene\ttotal_score\tcount\taverage_score\n";
    foreach my $gene (keys %{$score{$species}}) {
        my $average=$score{$species}{$gene}/$count{$species}{$gene};
        print $outhandle "$species\t$gene\t$score{$species}{$gene}\t$count{$species}{$gene}\t$average\n";
    }
    close $outhandle;
}







