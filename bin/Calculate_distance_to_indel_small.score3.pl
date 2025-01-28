#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Glob ':glob';

# Subroutine to read known genes from a file
sub read_known_genes {
    my ($file) = @_;
    my %known_genes;
    open my $fh, '<', $file or die "Could not open '$file': $!";
    while (my $gene = <$fh>) {
        chomp $gene;
        $known_genes{$gene} = 1;
    }
    close $fh;
    return \%known_genes;
}

# Subroutine to process the BED file and calculate scores
sub process_bed_file {
    my ($bed_file, $known_genes) = @_;
    my %gene_positions;

    # Read the BED file
    open my $bed_fh, '<', $bed_file or die "Could not open '$bed_file': $!";
    while (my $line = <$bed_fh>) {
        chomp $line;
        my ($chr, $start, $end, $gene) = split "\t", $line;
        push @{$gene_positions{$chr}}, $gene;
    }
    close $bed_fh;

    # Calculate scores
    my %gene_scores;
    foreach my $chr (keys %gene_positions) {
        my @genes = @{$gene_positions{$chr}};

        # Traverse each gene and calculate the minimum distance to a known gene
        for my $i (0 .. $#genes) {
            my $gene = $genes[$i];
            my $min_distance = 'NA';

            # If gene is a known gene, its score is 0
            if (exists $known_genes->{$gene}) {
                $gene_scores{"$chr\t$gene"} = 0;
                next;
            }

            # Calculate distance to nearest known gene
            for my $j (0 .. $#genes) {
                next if $i == $j;
                my $other_gene = $genes[$j];
                if (exists $known_genes->{$other_gene}) {
                    my $distance = abs($i - $j);
                    if ($min_distance eq 'NA' or $distance < $min_distance) {
                        $min_distance = $distance;
                    }
                }
            }

            # Store the calculated score
            $gene_scores{"$chr\t$gene"} = $min_distance;
        }
    }

    return \%gene_scores;
}

# Find all Translocation_gene_boundaries.txt files
my @translocation_files = bsd_glob("*indel_small.txt");

# Process each species
my %species_known_genes;
my $sp1;
my $sp2;
my $done=0;
my $type;
foreach my $file (@translocation_files) {
    $done=1;
    my ($species) = split(/\./, $file);
    push @{$species_known_genes{$species}}, read_known_genes($file);
    my @split=split(/\./, $file);
    $sp1=$split[0];
    $sp2=$split[1];
    $type=$split[3];
}

# Find all bed files
my @bed_files = bsd_glob("*.bed");

# Process each BED file and calculate scores
foreach my $bed_file (@bed_files) {
    my ($species) = fileparse($bed_file, qr/\.[^.]*/);
    my $known_genes_combined = {};

    # Combine known genes for the species
    if (exists $species_known_genes{$species}) {
        foreach my $known_genes (@{$species_known_genes{$species}}) {
            @$known_genes_combined{keys %$known_genes} = values %$known_genes;
        }
    } else {
        warn "No translocation files found for species '$species'. Skipping...\n";
        next;
    }

    # Process the BED file
    my $gene_scores = process_bed_file($bed_file, $known_genes_combined);

    # Print the results
    open my $output_fh, '>', "$sp1\.$sp2\.$type\.gene_scores.txt" or die "Could not open output file: $!";
    foreach my $key (sort keys %$gene_scores) {
        print $output_fh "$key\t$gene_scores->{$key}\n";
    }
    close $output_fh;

    print "Processed $bed_file and output to ${species}_gene_scores.txt\n";
}

if ($done){
    `echo $sp1 > species_tested_$sp1`
}
