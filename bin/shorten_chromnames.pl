#!/usr/bin/perl
use warnings;
use strict;

# Usage: shorten_chromnames.pl <seqids_file> <bed1> [bed2 ...]
#
# For each species (one line in seqids_file), tries to derive clean short labels:
#   1. Strip common biological prefixes (chr, scaffold_, LG, HiC_scaffold_, etc.)
#   2. If all resulting labels are <= 3 chars and unique within the species -> use them
#   3. Otherwise -> assign sequential integers (1, 2, 3, ...) for that species
#      and write a chromosome_name_legend.tsv mapping table
#
# Species order is read from species.csv (written by ribbon.pl in the same directory)
# so that each seqids line is matched to the correct species BED file even when
# multiple species share the same scaffold naming convention.
#
# Outputs:
#   <seqids_file>.short        renamed seqids (same row order)
#   <bed_file>.short           renamed BED files (one per input BED)
#   chromosome_name_legend.tsv (only written if any species needed sequential numbering)

die "Usage: shorten_chromnames.pl <seqids_file> <bed1> [bed2 ...]\n" if @ARGV < 2;

my ($seqids_file, @bed_files) = @ARGV;

my $MAX_LABEL_LEN = 3;

# Read species order from species.csv so we can map seqids lines to species names
my @species_order;
if (-e "species.csv") {
    open(my $SP, "<", "species.csv") or die "Cannot open species.csv: $!\n";
    my $line = <$SP>;
    close $SP;
    if (defined $line) {
        chomp $line;
        @species_order = split(/,/, $line);
        chomp $_ for @species_order;
    }
}

# Read seqids file: each non-blank line = one species, comma-separated chrom names
open(my $SEQ, "<", $seqids_file) or die "Cannot open $seqids_file: $!\n";
my @species_chroms;
while (my $line = <$SEQ>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    push @species_chroms, [split(/,/, $line)];
}
close $SEQ;

# Strip common biological prefixes from a chromosome name
sub strip_prefix {
    my $name = shift;
    (my $s = $name) =~ s/^(?:hic_scaffold|chromosome|scaffold|contig|chr|lg|un)[-_.]?//i;
    return length($s) ? $s : $name;
}

# Derive species name from a BED filename (strips path and .bed / .bed.flipped.bed)
sub species_from_bed {
    my $bed = shift;
    (my $sp = $bed) =~ s/\.bed(?:\.flipped\.bed)?$//i;
    $sp =~ s{.*/}{};
    $sp =~ s/\s+$//;
    return $sp;
}

# Build one name_map per species line
my @per_species_map;   # per_species_map[i] = { orig_chrom => short_label }
my @legend_entries;

for my $i (0 .. $#species_chroms) {
    my @chroms    = @{$species_chroms[$i]};
    my $sp_name   = ($i < @species_order) ? $species_order[$i] : "species_" . ($i + 1);

    # Try stripping common prefixes
    my @stripped = map { strip_prefix($_) } @chroms;

    # Accept only if every stripped label is <= MAX_LABEL_LEN and unique within species
    my %seen;
    my $use_stripped = 1;
    for my $s (@stripped) {
        if (length($s) > $MAX_LABEL_LEN || $seen{$s}++) {
            $use_stripped = 0;
            last;
        }
    }

    my %local_map;
    if ($use_stripped) {
        for my $j (0 .. $#chroms) {
            $local_map{$chroms[$j]} = $stripped[$j];
        }
    } else {
        for my $j (0 .. $#chroms) {
            my $label = $j + 1;
            $local_map{$chroms[$j]} = "$label";
            push @legend_entries, [$sp_name, $label, $chroms[$j]];
        }
    }

    push @per_species_map, \%local_map;
}

# Map species name -> index so BED files can be matched to the right per-species map
my %species_to_idx;
for my $i (0 .. $#species_order) {
    $species_to_idx{$species_order[$i]} = $i;
}

# Write renamed seqids file (one line per species, in original order)
open($SEQ, "<", $seqids_file) or die;
open(my $OUT_SEQ, ">", "$seqids_file.short") or die "Cannot open $seqids_file.short: $!\n";
my $line_idx = 0;
while (my $line = <$SEQ>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my $map = $per_species_map[$line_idx] // {};
    print $OUT_SEQ join(",", map { $map->{$_} // $_ } split(/,/, $line)) . "\n";
    $line_idx++;
}
close $SEQ;
close $OUT_SEQ;

# Write renamed BED files, using each species' own map
for my $bed (@bed_files) {
    my $sp = species_from_bed($bed);
    my $map;
    if (exists $species_to_idx{$sp}) {
        $map = $per_species_map[$species_to_idx{$sp}] // {};
    } else {
        warn "Warning: could not match BED file '$bed' to a species in species.csv; skipping rename\n";
        $map = {};
    }

    open(my $BED_IN, "<", $bed) or do { warn "Cannot open $bed: $!, skipping\n"; next; };
    open(my $BED_OUT, ">", "$bed.short") or die "Cannot open $bed.short: $!\n";
    while (my $line = <$BED_IN>) {
        chomp $line;
        my @f = split(/\t/, $line);
        $f[0] = $map->{$f[0]} // $f[0];
        print $BED_OUT join("\t", @f) . "\n";
    }
    close $BED_IN;
    close $BED_OUT;
}

# Write legend only when sequential numbering was used
if (@legend_entries) {
    open(my $LEG, ">", "chromosome_name_legend.tsv")
        or die "Cannot open chromosome_name_legend.tsv: $!\n";
    print $LEG "Species\tDisplay_label\tFull_name\n";
    for my $e (@legend_entries) {
        print $LEG join("\t", @$e) . "\n";
    }
    close $LEG;
}
