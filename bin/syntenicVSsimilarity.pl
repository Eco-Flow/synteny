use warnings;
use strict;


print "Please be in folder with all the similarity and synteny files\n";

my $outfile="My_comp_synteny_similarity.tsv";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
my $outfile2="PairwisePlotData.tsv";
open(my $outhandle2, ">", $outfile2)   or die "Could not open $outfile2 \n";
my $outfile3="R_commands.txt";
open(my $outhandle3, ">", $outfile3)   or die "Could not open $outfile3 \n";

my $in1= "Synteny_matrix.tsv";
my $in2= "My_sim_cores.tsv";
open(my $filein1, "<", $in1)   or die "Could not open $in1\n";
open(my $filein2, "<", $in2)   or die "Could not open $in2\n";

my %FILE1;


my $headM=<$filein1>;
chomp $headM;
my @sp=split("\t", $headM);
my $waste=shift(@sp);
my %file1store;

while (my $lineM=<$filein1>){
    chomp $lineM;
    print "line here: $lineM\n";
    my @spline=split("\t", $lineM);
    my $species=shift(@spline);
    my $n=0;
    foreach my $score (@spline){
        my $sp1=$sp[$n];
        my $sp2=$species;
        print "1 $sp1 2 $sp2 3 $score\n";
        $file1store{$sp1}{$sp2}=$score;
        $n++;
    }
}

my $headSy=<$filein2>;
chomp $headSy;
my @sp2=split("\t", $headSy);
my $waste2=shift(@sp2);
my %file2store;

while (my $lineSy=<$filein2>){
    chomp $lineSy;
    my @spline2=split("\t", $lineSy);
    my $species=shift(@spline2);
    my $n=0;
    foreach my $score (@spline2){
        my $sp1=$sp[$n];
        my $sp2=$species;
        $file2store{$sp1}{$sp2}=$score;
        $n++;
    }
}



print $outhandle "$headSy\n";

#Now compare the two numbers:
print $outhandle2 "Pair\tSyntenous genes\tPercent protein similarity\n";

my %done_recip;

foreach my $species1 (@sp){
    print $outhandle "$species1";
    foreach my $species2 (@sp){
        print $outhandle "\t$file1store{$species1}{$species2}\/$file2store{$species1}{$species2}";


        #Then print a pair wise dataset
        my $joint = "$species1$species2";
        my $joinr = "$species2$species1";
        if ($species1 eq $species2){
            #do nothing
        }
        elsif ($done_recip{$joint} || $done_recip{$joinr}){
            #do nothing
        }
        else{
            print $outhandle2 "$joint\t$file1store{$species1}{$species2}\t$file2store{$species1}{$species2}\n";
            $done_recip{$joint}="DONE";
            $done_recip{$joinr}="DONE";

        }
        
    }
    print $outhandle "\n";
}


print $outhandle3 "
mydata <- read.csv(\"PairwisePlotData.tsv\", h=T, row.names=1, sep=\"\t\")
pdf (\"My_pair_synteny_identity.pdf\")
plot(mydata)
text(mydata\$Syntenous.genes, mydata\$Percent.protein.similarity, labels=rownames(mydata))
dev.off()
";

`R --vanilla < R_commands.txt > output.ofthis.test`;

close $outhandle;
close $outhandle2;
close $outhandle3;



