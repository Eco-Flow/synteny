use warnings;
use strict;

my $in_name="files_in";
open(my $IN_N, "<", $in_name)   or die "Could not open $in_name\n";

my $input=<$IN_N>;
$input =~ s|\[||g;
$input =~ s|\]||g;

print "HERE $input\n";

my @ingos=split("\, ", $input);

foreach my $paths_to_run (@ingos){
    `ln -s $paths_to_run .`;
}

#set the percent valued used in this test;
my $comparison_percentage;

print "Please be in folder with all the Species Go Summarys\n";

my %go_key;
my %go_names;
my %species_list;
my @gos=`ls *results_ALL.tab`;
my $species;
my $subset;
foreach my $sp (@gos){
    chomp $sp;
    print "RUN $sp\n";
    #my @pathsp=split(/\//, $sp);
    my @split=split(/\./, $sp);
    #my @sp_folder=split("\/", $split[0]);
    $species=$split[0];
    $comparison_percentage=$split[1];
    $subset=$split[2];
    #print "$species $subset\n";
    $species_list{$species}="Done";
    open(my $filein, "<", $sp)   or die "Could not open $sp\n";
    my $header=<$filein>;
    while (my $line = <$filein>){
        chomp $line;
        #print "$line\n";
        my @split2=split("\t", $line);
        $go_key{$split2[0]}{$subset}{$species}=$split2[9];
        $go_names{$split2[0]}=$split2[1];
        print "$split2[0] $subset $species $split2[9]\n";
    }
}

#Create a species array:
my @species_array;
foreach my $sp (keys %species_list) {
    push (@species_array, $sp);
}

#Create job type array:
my @job_type_array=("closest","farthest");

#Create output files
my $outname="Go_summary_$comparison_percentage\_all_merged.tsv";
open(my $outhandle, ">", $outname)   or die "Could not open $outname\n";
my $outname2="Go_summary_$comparison_percentage\_closest.tsv";
open(my $outhandle2, ">", $outname2)   or die "Could not open $outname2\n";
my $outname3="Go_summary_$comparison_percentage\_farthest.tsv";
open(my $outhandle3, ">", $outname3)   or die "Could not open $outname3\n";


print $outhandle "GO_ID\tGO_term";
print $outhandle2 "GO_ID\tGO_term";
print $outhandle3 "GO_ID\tGO_term";


foreach my $types ( @job_type_array ) {
    foreach my $species ( @species_array ) {
        print $outhandle "\t$species\_$types";
    }
}

foreach my $species ( @species_array ) {
	print $outhandle2 "\t$species";
    print $outhandle3 "\t$species";
}

print $outhandle "\tCount_closest\tCount_farthest\tCount_significant\n";
print $outhandle2 "\tCount_significant\n";
print $outhandle3 "\tCount_significant\n";

foreach my $goterms (keys %go_key) {
    print $outhandle "$goterms\t$go_names{$goterms}";
    print $outhandle2 "$goterms\t$go_names{$goterms}";
    print $outhandle3 "$goterms\t$go_names{$goterms}";

    #All job_type_array
    #My tests_signif count for all sampled:
    my %count_pval;

    #For each test (e.g. closestSynteny","botSynteny)
    foreach my $types ( @job_type_array ) {

        foreach my $species ( @species_array ) {
            if (exists $go_key{$goterms}{$types}{$species} ){
                print $outhandle "\t$go_key{$goterms}{$types}{$species}";

                #if p value is less than 0.05, add to counts for GO term:
                if ($go_key{$goterms}{$types}{$species} <= 0.05){
                    if ($count_pval{$goterms}{$types}){
                        $count_pval{$goterms}{$types}++;
                    }
                    else{
                        $count_pval{$goterms}{$types}=1;
                    }

                }

                if($types eq "closest"){
                	print $outhandle2 "\t$go_key{$goterms}{$types}{$species}";
                }
                if($types eq "farthest"){
                	print $outhandle3 "\t$go_key{$goterms}{$types}{$species}";
                }                 
            }
            else{
            	#ALL types:
                print $outhandle "\tNA";

                if($types eq "closest"){
                	print $outhandle2 "\tNA";
                }
                if($types eq "farthest"){
                	print $outhandle3 "\tNA";
                }
            }
        }



	    #Now print off the significant counts per row,, for each output table:
	    if($types eq "closest"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle2 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle2 "\t0";
	    	}
	    }
	    if($types eq "farthest"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle3 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle3 "\t0";
	    	}
	    }

    }

    #Now print off the significant counts per row:
    my $comb=0;
    foreach my $types ( @job_type_array ) {
        if ($count_pval{$goterms}{$types}){
            print $outhandle "\t$count_pval{$goterms}{$types}";
            $comb+=$count_pval{$goterms}{$types};
        }
        else{
            print $outhandle "\t0";
        }

    }

    print $outhandle "\t$comb\n";
    print $outhandle2 "\n";
    print $outhandle3 "\n";

}

