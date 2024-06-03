#!/usr/bin/perl
use warnings;
use strict;

my @files=`ls *.anchors`;

foreach my $file (@files){
	chomp $file;
	print "To convert: $file\n";
	#Make sure we have a simple anchors file:
	`python -m jcvi.compara.synteny screen --minspan=30 --simple $file $file\.new`;
}
