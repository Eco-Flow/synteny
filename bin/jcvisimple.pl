#!/usr/bin/perl
use warnings;
use strict;

my @files=`ls *.anchors`;

my $args=$ARGV[0];  # add in jcvo.compara.synteny options ${params.jcvi_screen_arguments}
open(my $extra_arg, "<", $args)   or die "Could not open $args \n";
my $inarg=<$extra_arg>;
chomp $inarg;

foreach my $file (@files){
	chomp $file;
	print "To convert: $file\n";
	#Make sure we have a simple anchors file:
	`python -m jcvi.compara.synteny screen $inarg --simple $file $file\.new`;
	print "python -m jcvi.compara.synteny screen $inarg --simple $file $file\.new\n";
}
