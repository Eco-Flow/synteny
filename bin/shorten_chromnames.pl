#!/usr/bin/perl
use warnings;
use strict;

# Usage: shorten_chromnames.pl <seqids_file> <bed1> [bed2 ...]
#
# Shortens chromosome names to the minimum number of trailing characters
# needed to remain unique across all chromosomes (minimum 3 chars).
# Example: NC_053021.2 -> 21.2  (because 1.2 was already taken by NC_052531.2)
#
# Outputs:
#   <seqids_file>.short        renamed seqids (same row order)
#   <bed_file>.short           renamed BED files (one per input BED)

die "Usage: shorten_chromnames.pl <seqids_file> <bed1> [bed2 ...]\n" if @ARGV < 2;

my ($seqids_file, @bed_files) = @ARGV;

# Collect all unique chromosome names from the seqids file
open(my $SEQ, "<", $seqids_file) or die "Cannot open $seqids_file: $!\n";
my @all_chroms;
my %seen_orig;
while (my $line = <$SEQ>) {
    chomp $line;
    for my $chr (split(/,/, $line)) {
        push @all_chroms, $chr unless $seen_orig{$chr}++;
    }
}
close $SEQ;

# Build shortest unique names: try trailing substrings of increasing length
my %used_short;
my %name_map;

for my $name (@all_chroms) {
    if (length($name) <= 3) {
        $name_map{$name} = $name;
        $used_short{$name}++;
    } else {
        my $chosen = $name;  # fallback: full name unchanged
        for my $len (3 .. length($name)) {
            my $candidate = substr($name, -$len);
            unless ($used_short{$candidate}) {
                $chosen = $candidate;
                last;
            }
        }
        $name_map{$name} = $chosen;
        $used_short{$chosen}++;
    }
}

# Write renamed seqids file
open($SEQ, "<", $seqids_file) or die;
open(my $OUT_SEQ, ">", "$seqids_file.short") or die "Cannot open $seqids_file.short: $!\n";
while (my $line = <$SEQ>) {
    chomp $line;
    print $OUT_SEQ join(",", map { $name_map{$_} // $_ } split(/,/, $line)) . "\n";
}
close $SEQ;
close $OUT_SEQ;

# Write renamed BED files
for my $bed (@bed_files) {
    open(my $BED_IN, "<", $bed) or do {
        warn "Cannot open $bed: $!, skipping\n";
        next;
    };
    open(my $BED_OUT, ">", "$bed.short") or die "Cannot open $bed.short: $!\n";
    while (my $line = <$BED_IN>) {
        chomp $line;
        my @f = split(/\t/, $line);
        $f[0] = $name_map{$f[0]} // $f[0];
        print $BED_OUT join("\t", @f) . "\n";
    }
    close $BED_IN;
    close $BED_OUT;
}
