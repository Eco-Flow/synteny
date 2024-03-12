use warnings;
use strict;

my @files=`ls *.anchors`;

#Initiate hash for matrix store
my %in_syn_block;
my %last;

my $syn_block=0;

#Loop thru the files to get the synteny info.
foreach my $file (@files){
    chomp $file;
	my @split_name=split(/\./, $file);
	my $sp1=$split_name[0];
	my $sp2=$split_name[1];
    %in_syn_block=();
    %last=();
    my $outfile="$sp1\.$sp2\.geneScore.tsv";
    open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
    print $outhandle "Start\n";

	#Now per file:
    open(my $filein, "<", $file)   or die "Could not open $file\n";
    my $genecount=0;
    
    while (my $line=<$filein>){
        chomp $line;
        print "$sp1 $sp2 $line\n";
        if ($line =~ m/^#/){
            #new orthoblock starts
            $syn_block++;
            my $current_synt_block;
            # now loop thru the previous syntenous block and count backward 
            foreach my $sp_last (keys %last) {
				#print $outhandle "$key";
				$current_synt_block=$sp_last;
				foreach my $gene_last (keys %{$last{$sp_last}}) {
					my $value= $last{$sp_last}{$gene_last};
					my $total_rec=$genecount;
					my $halftot=$total_rec/2;
					if ($value >= $halftot){
						my $new_val=1+($total_rec-$value);
						$in_syn_block{$sp_last}{$gene_last}=$new_val;
					}
					
				}
			}
			#now set the gene count back to zeor, for next syntenous block.
            $genecount=0;
            #Now remove the last hash of genes
            #print "IS : $current_synt_block\n";
            undef %last;

        }
        else{
        	my @sp_gene_line=split("\t", $line);
        	my $sp1_id=$sp_gene_line[0];
        	my $sp2_id=$sp_gene_line[1];
        	$genecount++;
        	$in_syn_block{$sp1}{$sp1_id}=$genecount;
        	$last{$sp1}{$sp1_id}=$genecount;
        }
    }

    foreach my $sp_fin (keys %in_syn_block) {
        foreach my $gene_fin (keys %{$in_syn_block{$sp_fin}}) {
            print $outhandle "$sp_fin\t$gene_fin\t$in_syn_block{$sp_fin}{$gene_fin}\n";
        }
    }

    %in_syn_block=();

}




