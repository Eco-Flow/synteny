#!/usr/bin/env python3

"""
Build a concatenated (supermatrix) protein alignment from per-orthogroup
alignments, for species tree inference.

Takes a directory of aligned FASTA files, one per orthogroup, whose sequences are
named after their species (as written by extract_single_copy.py), plus the
OrthoFinder Orthogroups.tsv giving the expected species list.

An orthogroup is used only if it is strictly single-copy and complete: every
species must be represented exactly once. Those orthogroups are concatenated in
a stable order and a RAxML-style partition file is emitted alongside, so
IQ-TREE can fit a separate model per orthogroup.

Pure standard library (no BioPython/ete3) to match the rest of bin/.
"""

import argparse
import os
import sys
from collections import OrderedDict


def read_fasta(path):
    """Read a FASTA file into an OrderedDict of {header_first_token: sequence}."""
    seqs = OrderedDict()
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def read_species(orthogroups_tsv):
    """Return the species list from the column headers of Orthogroups.tsv."""
    with open(orthogroups_tsv) as fh:
        header = fh.readline().rstrip('\n\r').split('\t')
    # First column is the orthogroup ID; the rest are species.
    return header[1:]


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate single-copy orthogroup alignments into a supermatrix'
    )
    parser.add_argument('-m', '--msa-dir', required=True,
                        help='Directory of per-orthogroup aligned FASTA files')
    parser.add_argument('-g', '--orthogroups', required=True,
                        help='OrthoFinder Orthogroups.tsv (used to map genes to species)')
    parser.add_argument('-o', '--out-fasta', required=True,
                        help='Output concatenated FASTA alignment')
    parser.add_argument('-p', '--out-partitions', required=True,
                        help='Output RAxML-style partition file')
    parser.add_argument('--model', default='AA',
                        help="Model/datatype field written for each partition (default: AA, "
                             "letting IQ-TREE's ModelFinder choose per partition)")
    parser.add_argument('--min-orthogroups', type=int, default=1,
                        help='Fail if fewer than this many usable orthogroups are found (default: 1)')

    args = parser.parse_args()

    if not os.path.isdir(args.msa_dir):
        sys.exit("ERROR: alignment directory '%s' not found." % args.msa_dir)

    species = read_species(args.orthogroups)
    if not species:
        sys.exit("ERROR: no species columns found in %s" % args.orthogroups)
    species = sorted(species)
    species_set = set(species)

    alignment_files = sorted(
        os.path.join(args.msa_dir, f)
        for f in os.listdir(args.msa_dir)
        if f.endswith(('.fa', '.faa', '.fasta', '.aln'))
    )
    if not alignment_files:
        sys.exit("ERROR: no alignment files found in %s" % args.msa_dir)

    # Accumulate per-species sequence blocks plus the partition ranges.
    blocks = OrderedDict((sp, []) for sp in species)
    partitions = []
    offset = 0
    skipped_multicopy = 0
    skipped_incomplete = 0
    skipped_ragged = 0

    for path in alignment_files:
        og = os.path.splitext(os.path.basename(path))[0]
        seqs = read_fasta(path)
        if not seqs:
            continue

        # Group this orthogroup's sequences by species.
        by_species = {}
        unmapped = False
        for name, seq in seqs.items():
            if name not in species_set:
                unmapped = True
                break
            by_species.setdefault(name, []).append(seq)
        if unmapped:
            skipped_incomplete += 1
            continue

        if any(len(v) > 1 for v in by_species.values()):
            skipped_multicopy += 1
            continue
        if len(by_species) != len(species):
            skipped_incomplete += 1
            continue

        lengths = {len(v[0]) for v in by_species.values()}
        if len(lengths) != 1:
            # Not a rectangular alignment — should not happen, but never
            # silently corrupt the supermatrix coordinates.
            skipped_ragged += 1
            continue
        width = lengths.pop()

        for sp in species:
            blocks[sp].append(by_species[sp][0])
        partitions.append((og, offset + 1, offset + width))
        offset += width

    if len(partitions) < args.min_orthogroups:
        sys.exit(
            "ERROR: only %d usable single-copy orthogroup(s) found across %d species "
            "(minimum required: %d). A supermatrix cannot be built."
            % (len(partitions), len(species), args.min_orthogroups)
        )

    with open(args.out_fasta, 'w') as fh:
        for sp in species:
            fh.write('>%s\n' % sp)
            seq = ''.join(blocks[sp])
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')

    with open(args.out_partitions, 'w') as fh:
        for og, start, end in partitions:
            fh.write('%s, %s = %d-%d\n' % (args.model, og, start, end))

    print("Supermatrix built successfully")
    print("   Species:                 %d" % len(species))
    print("   Orthogroups used:        %d" % len(partitions))
    print("   Alignment length:        %d" % offset)
    print("   Skipped (multi-copy):    %d" % skipped_multicopy)
    print("   Skipped (incomplete):    %d" % skipped_incomplete)
    if skipped_ragged:
        print("   Skipped (ragged):        %d" % skipped_ragged)
    print("   Alignment written to:    %s" % args.out_fasta)
    print("   Partitions written to:   %s" % args.out_partitions)


if __name__ == '__main__':
    main()
