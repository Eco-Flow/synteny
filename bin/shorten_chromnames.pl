#!/usr/bin/perl
use warnings;
use strict;

# Usage: shorten_chromnames.pl <seqids_file> <bed1> [bed2 ...]
#
# For each species (one line in seqids_file), tries to derive clean short labels:
#   1. Strip common biological prefixes (chr, scaffold_, LG, etc.)
#   2. If all resulting labels are <= 3 chars and unique within the species -> use them
#   3. Otherwise -> assign sequential integers (1, 2, 3, ...) for that species
#      and write a chromosome_name_legend.tsv mapping table
#
# Outputs:
#   <seqids_file>.short        renamed seqids (same row order)
#   <bed_file>.short           renamed BED files (one per input BED)
#   chromosome_name_legend.tsv (only written if any species needed sequential numbering)

die "Usage: shorten_chromnames.pl <seqids_file> <bed1> [bed2 ...]\n" if @ARGV < 2;

my ($seqids_file, @bed_files) = @ARGV;

my $MAX_LABEL_LEN = 3;

# Read seqids file: each line = one species, comma-separated chromosome names
open(my $SEQ, "<", $seqids_file) or die "Cannot open $seqids_file: $!\n";
my @species_chroms;
while (my $line = <$SEQ>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    push @species_chroms, [split(/,/, $line)];
}
close $SEQ;

# Read chromosome sets from each BED file so we can identify which BED = which seqids line
my %bed_chrom_set;
for my $bed (@bed_files) {
    open(my $B, "<", $bed) or do { warn "Cannot open $bed: $!\n"; next; };
    while (<$B>) {
        chomp;
        my ($chrom) = split(/\t/, $_, 2);
        $bed_chrom_set{$bed}{$chrom} = 1 if defined $chrom;
    }
    close $B;
}

# Derive a display-friendly species label from a BED filename
sub species_label {
    my $bed = shift;
    (my $sp = $bed) =~ s/\.bed(?:\.flipped\.bed)?$//i;
    $sp =~ s{.*/}{};
    return $sp;
}

# Strip common biological prefixes from a chromosome name
sub strip_prefix {
    my $name = shift;
    (my $s = $name) =~ s/^(?:chromosome|scaffold|contig|chr|lg|un)[-_.]?//i;
    return length($s) ? $s : $name;  # never return empty string
}

# Build the global name_map and collect legend entries
my %name_map;
my @legend_entries;

for my $i (0 .. $#species_chroms) {
    my @chroms = @{$species_chroms[$i]};

    # Identify the species name by matching chromosomes to a BED file
    my $species_name = "species_" . ($i + 1);
    for my $bed (@bed_files) {
        my $hits = grep { $bed_chrom_set{$bed}{$_} } @chroms;
        if ($hits > 0) {
            $species_name = species_label($bed);
            last;
        }
    }

    # Try prefix-stripping
    my @stripped = map { strip_prefix($_) } @chroms;

    # Accept stripped names only if ALL are <= MAX_LABEL_LEN and unique within species
    my %seen;
    my $use_stripped = 1;
    for my $s (@stripped) {
        if (length($s) > $MAX_LABEL_LEN || $seen{$s}++) {
            $use_stripped = 0;
            last;
        }
    }

    if ($use_stripped) {
        for my $j (0 .. $#chroms) {
            $name_map{$chroms[$j]} = $stripped[$j];
        }
    } else {
        # Fall back to sequential integers and record in legend
        for my $j (0 .. $#chroms) {
            my $label = $j + 1;
            $name_map{$chroms[$j]} = "$label";
            push @legend_entries, [$species_name, $label, $chroms[$j]];
        }
    }
}

# Write renamed seqids file
open($SEQ, "<", $seqids_file) or die;
open(my $OUT_SEQ, ">", "$seqids_file.short") or die "Cannot open $seqids_file.short: $!\n";
while (my $line = <$SEQ>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    print $OUT_SEQ join(",", map { $name_map{$_} // $_ } split(/,/, $line)) . "\n";
}
close $SEQ;
close $OUT_SEQ;

# Write renamed BED files
for my $bed (@bed_files) {
    open(my $BED_IN, "<", $bed) or do { warn "Cannot open $bed: $!, skipping\n"; next; };
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

# Write legend only if sequential numbering was needed
if (@legend_entries) {
    open(my $LEG, ">", "chromosome_name_legend.tsv")
        or die "Cannot open chromosome_name_legend.tsv: $!\n";
    print $LEG "Species\tDisplay_label\tFull_name\n";
    for my $e (@legend_entries) {
        print $LEG join("\t", @$e) . "\n";
    }
    close $LEG;
}
