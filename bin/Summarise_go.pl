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
my @job_type_array=("topSynteny","botSynteny","averhigh","averlow","highScore","lowScore");

#Create output files
my $outname="Go_summary_$comparison_percentage\_all_merged.tsv";
open(my $outhandle, ">", $outname)   or die "Could not open $outname\n";
my $outname2="Go_summary_$comparison_percentage\_topSynteny.tsv";
open(my $outhandle2, ">", $outname2)   or die "Could not open $outname2\n";
my $outname3="Go_summary_$comparison_percentage\_botSynteny.tsv";
open(my $outhandle3, ">", $outname3)   or die "Could not open $outname3\n";
my $outname4="Go_summary_$comparison_percentage\_averhigh.tsv";
open(my $outhandle4, ">", $outname4)   or die "Could not open $outname4\n";
my $outname5="Go_summary_$comparison_percentage\_averlow.tsv";
open(my $outhandle5, ">", $outname5)   or die "Could not open $outname5\n";
my $outname6="Go_summary_$comparison_percentage\_highScore.tsv";
open(my $outhandle6, ">", $outname6)   or die "Could not open $outname6\n";
my $outname7="Go_summary_$comparison_percentage\_lowScore.tsv";
open(my $outhandle7, ">", $outname7)   or die "Could not open $outname7\n";

print $outhandle "GO_ID\tGO_term";
print $outhandle2 "GO_ID\tGO_term";
print $outhandle3 "GO_ID\tGO_term";
print $outhandle4 "GO_ID\tGO_term";
print $outhandle5 "GO_ID\tGO_term";
print $outhandle6 "GO_ID\tGO_term";
print $outhandle7 "GO_ID\tGO_term";


foreach my $types ( @job_type_array ) {
    foreach my $species ( @species_array ) {

        print $outhandle "\t$species\_$types";

    }
}

foreach my $species ( @species_array ) {
	print $outhandle2 "\t$species";
    print $outhandle3 "\t$species";
    print $outhandle4 "\t$species";
    print $outhandle5 "\t$species";
    print $outhandle6 "\t$species";
    print $outhandle7 "\t$species";
}

print $outhandle "\tCount_topSynteny\tCount_botSynteny\tCount_averhigh\tCount_averlow\tCount_highScore\tCount_lowScore\n";
print $outhandle2 "\tCount_significant\n";
print $outhandle3 "\tCount_significant\n";
print $outhandle4 "\tCount_significant\n";
print $outhandle5 "\tCount_significant\n";
print $outhandle6 "\tCount_significant\n";
print $outhandle7 "\tCount_significant\n";


foreach my $goterms (keys %go_key) {
    print $outhandle "$goterms\t$go_names{$goterms}";
    print $outhandle2 "$goterms\t$go_names{$goterms}";
    print $outhandle3 "$goterms\t$go_names{$goterms}";
    print $outhandle4 "$goterms\t$go_names{$goterms}";
    print $outhandle5 "$goterms\t$go_names{$goterms}";
    print $outhandle6 "$goterms\t$go_names{$goterms}";
    print $outhandle7 "$goterms\t$go_names{$goterms}";
    
    #All job_type_array
    #My tests_signif count for all sampled:
    my %count_pval;

    #For each test (e.g. topSynteny","botSynteny)
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

                if($types eq "topSynteny"){
                	print $outhandle2 "\t$go_key{$goterms}{$types}{$species}";
                }
                if($types eq "botSynteny"){
                	print $outhandle3 "\t$go_key{$goterms}{$types}{$species}";
                }    
                if($types eq "averhigh"){
                	print $outhandle4 "\t$go_key{$goterms}{$types}{$species}";
                }
                if($types eq "averlow"){
                	print $outhandle5 "\t$go_key{$goterms}{$types}{$species}";
                }
                if($types eq "highScore"){
                	print $outhandle6 "\t$go_key{$goterms}{$types}{$species}";
                }
                if($types eq "lowScore"){
                	print $outhandle7 "\t$go_key{$goterms}{$types}{$species}";
                }                     
            }
            else{
            	#ALL types:
                print $outhandle "\tNA";

                if($types eq "topSynteny"){
                	print $outhandle2 "\tNA";
                }
                if($types eq "botSynteny"){
                	print $outhandle3 "\tNA";
                }
                if($types eq "averhigh"){
                	print $outhandle4 "\tNA";
                }
                if($types eq "averlow"){
                	print $outhandle5 "\tNA";
                }                
                if($types eq "highScore"){
                	print $outhandle6 "\tNA";
                }
                if($types eq "lowScore"){
                	print $outhandle7 "\tNA";
                }  
            }
        }

	    #Now print off the significant counts per row:
	    foreach my $types ( @job_type_array ) {
	        if ($count_pval{$goterms}{$types}){
	            print $outhandle "\t$count_pval{$goterms}{$types}";
	        }
	        else{
	            print $outhandle "\t0";
	        }

	    }

	    #Now do for each output table:
	    if($types eq "topSynteny"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle2 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle2 "\t0";
	    	}
	    }
	    if($types eq "botSynteny"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle3 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle3 "\t0";
	    	}
	    }
	    if($types eq "averhigh"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle4 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle4 "\t0";
	    	}
	    }
	    if($types eq "averlow"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle5 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle5 "\t0";
	    	}
	    }                
	    if($types eq "highScore"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle6 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle6 "\t0";
	    	}
	    }
	    if($types eq "lowScore"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle7 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle7 "\t0";
	    	}
	    }


    }

    


    print $outhandle "\n";
    print $outhandle2 "\n";
    print $outhandle3 "\n";
    print $outhandle4 "\n";
    print $outhandle5 "\n";
    print $outhandle6 "\n";
    print $outhandle7 "\n";

}




