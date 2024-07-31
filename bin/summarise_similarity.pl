use warnings;
use strict;

# This script is to summarise the protein similarity stats 
# of all the pairwise comparisons. This results in a matrix,
# called "My_sim_cores.tsv", which shows the average protein 
# percentage identity across all pairwise orthologous genes.
# This is used as a basic measure of phylogenetic distance

print "Please be in folder with all the similarity files\n";

my $outfile="My_sim_cores.tsv";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
my @files=`ls *.percent.similarity`;

#Initialise hashes
my %done;
my %sim_hash;

#Loop thru files and save the similarity as the value of a hash of the two species IDs.
foreach my $file (@files){
    chomp $file;
	my @split_name=split(/\./, $file);
	my $sp1=$split_name[0];
	my $sp2=$split_name[1];
    my $species="$sp1$sp2";

    open(my $filein, "<", $file)   or die "Could not open $file\n";
    my $similarity=<$filein>;
    chomp $similarity;
    $sim_hash{$sp1}{$sp2}=$similarity;
}

# Now print this into a matrix of protein similarity scores.
foreach my $sp1 (sort keys %sim_hash ){
    print $outhandle "\t$sp1";
}
print $outhandle "\n";

foreach my $sp1 (sort keys %sim_hash ){
    print $outhandle "$sp1";
    foreach my $sp2 (sort keys %sim_hash ){
        if ($sp1 eq $sp2){
            print $outhandle "\tNA";
        }
        else{
            print $outhandle "\t$sim_hash{$sp1}{$sp2}";
        }
    }

    print $outhandle "\n";
}
