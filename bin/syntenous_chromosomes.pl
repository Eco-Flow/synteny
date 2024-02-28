use warnings;
use strict;

die "Provide the name of the first bed, second bed, anchor (.new) file file\n" if (@ARGV!=3); 

print "Script is running\n";

my $BED_file1 = $ARGV[0];
open(my $BED1, "<", $BED_file1) or die "Could not open $BED_file1\n";
my $BED_file2 = $ARGV[1];
open(my $BED2, "<", $BED_file2) or die "Could not open $BED_file2\n";
my $ANCHORS = $ARGV[2];
open(my $ANCH, "<", $ANCHORS) or die "Could not open $ANCHORS\n";

my $outfile="seqids_karyotype.txt";
open(my $out, "> $outfile") or die "error opening $outfile. $!";

#read in the first bed file
my %SP1;
while (my $line=<$BED1>) {
    chomp $line;
	my @spl=split ("\t", $line);
    $SP1{$spl[3]}=$spl[0];
}

#read in the secon bed file
my %SP2;
while (my $line2=<$BED2>) {
    chomp $line2;
	my @spl=split ("\t", $line2);
    #print "$spl[0] and  $spl[3]\n";
    $SP2{$spl[3]}=$spl[0];
}

#Calc the scaffold name in the anchors file
my %chromo_pairs;

while (my $line3=<$ANCH>) {
    chomp $line3;
    my @spl=split ("\t", $line3);
    my $len=scalar(@spl);
    if ( $len == 3 ){
        my $chr1=$SP1{"$spl[0]"};
        my $chr2=$SP2{"$spl[1]"};

        #print "$line3\nch1 $chr1 ch2 $chr2  do they exixsts, and if not, why not\n";
        my $join="$chr1\_\@\_$chr2";
        $chromo_pairs{"$join"}++;
    }
    else{
        #probably a # start line between syntenic regions.
    }
}

my %done_check_sp1;
my %done_check_sp2;
my @final_list_sp1;
my @final_list_sp2;
foreach my $key ( sort {$chromo_pairs{$b} <=> $chromo_pairs{$a}} keys %chromo_pairs ) {

    print "$key $chromo_pairs{$key}\n";
    my @spl=split("\_\@\_", $key);
    if ($done_check_sp1{$spl[0]}){
        #do nothing
    }
    else{
        push (@final_list_sp1, "$spl[0]");
        $done_check_sp1{$spl[0]}="yes"
    }

    if ($done_check_sp2{$spl[1]}){
        #do nothing
    }
    else{
        push (@final_list_sp2, "$spl[1]");
        $done_check_sp2{$spl[1]}="yes"
    }
}

my $join_sp1=join("\,", @final_list_sp1);
my $join_sp2=join("\,", @final_list_sp2);

print $out "$join_sp1\n$join_sp2\n";
