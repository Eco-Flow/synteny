#!/usr/bin/env python3
"""
Combine AGORA's ancestral gene-order reconstruction (contiguous ancestral
regions, CARs) with each species' extant BUSCO gene order to estimate a
gene-order fragmentation index (FI) per extant chromosome.

NOTE: the manuscript (Maulana et al. 2026, bioRxiv 2026.07.17.739156) gives
the FI formula as a typeset equation that was not recoverable as text from
the supplied methods (it renders as an "Embedded Image" placeholder). This
script computes and reports the three raw components the paper describes --
B (observed gene-order blocks), A (AGORA ancestral units on the chromosome),
M (informative BUSCO markers on the chromosome) -- plus a PLACEHOLDER FI
column using a plausible (B - A) / (M - A) rescaling. Treat the FI column as
provisional until the real equation from the paper/supplement is substituted.
"""
import argparse
import bz2
import glob
import gzip
import lzma
import os
from collections import defaultdict


def open_maybe_compressed(path):
    """AGORA compresses its output by default (bz2, confirmed against a real
    run; HowTo.md also documents gzip/lzma/xz as possible), so ancGenome files
    need transparent decompression regardless of extension."""
    if path.endswith(".bz2"):
        return bz2.open(path, "rt")
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith((".xz", ".lzma")):
        return lzma.open(path, "rt")
    return open(path)


def load_ancgenome_assignments(ancgenome_paths):
    """Map extant gene id -> ancestral block/CAR name from AGORA's ancGenome
    output. Format per row (confirmed against a real AGORA run): block_name,
    start, end, orientation, space-separated gene entries -- the first is the
    ancestral gene id, the rest are "<species>.<species>_<busco_id>" (AGORA
    prefixes each extant copy with its species name and a literal '.'; split
    on the first '.', not '_', since species names can themselves contain
    underscores)."""
    assignment = {}
    for path in ancgenome_paths:
        with open_maybe_compressed(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 5:
                    continue
                block = fields[0]
                entries = fields[4].split()
                for entry in entries[1:]:
                    gene_id = entry.split(".", 1)[1] if "." in entry else entry
                    assignment[gene_id] = block
    return assignment


def load_filtered_table(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            busco_id, sequence, start, end, strand = line.split("\t")
            rows.append((busco_id, sequence, int(start)))
    return rows


def count_blocks(labels):
    """Number of maximal runs of identical consecutive labels."""
    if not labels:
        return 0
    blocks = 1
    for prev, cur in zip(labels, labels[1:]):
        if cur != prev:
            blocks += 1
    return blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tables",
        nargs="+",
        required=True,
        help="species=path pairs to filtered BUSCO tables (busco_id, "
        "sequence, start, end, strand)",
    )
    parser.add_argument(
        "--ancgenome-glob",
        required=True,
        help="Glob pattern matching AGORA ancGenome output file(s). AGORA writes one "
        "file per reconstructed ancestor; the paper's FI is defined relative to a "
        "single target ancestor (\"the last common ancestor\"), so this should "
        "usually be scoped to that one ancestor's file, e.g. "
        "'**/ancGenome.<root_node_name>.list*' -- if the glob matches files from "
        "multiple ancestor levels, a gene present in more than one (e.g. under both "
        "a shallower and a deeper ancestor) gets whichever assignment is loaded "
        "last, in sorted-path order, not necessarily the most recent one",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ancgenome_paths = sorted(
        p for p in glob.glob(args.ancgenome_glob, recursive=True) if os.path.isfile(p)
    )
    assignment = load_ancgenome_assignments(ancgenome_paths)

    with open(args.output, "w") as out:
        out.write("species\tchromosome\tM\tA\tB\tFI_placeholder\n")
        for pair in args.tables:
            species, path = pair.split("=", 1)
            rows = load_filtered_table(path)

            by_chrom = defaultdict(list)
            for busco_id, sequence, start in rows:
                by_chrom[sequence].append((start, busco_id))

            for chromosome, markers in by_chrom.items():
                markers.sort(key=lambda x: x[0])
                labels = []
                for _, busco_id in markers:
                    gene_id = f"{species}_{busco_id}"
                    block = assignment.get(gene_id)
                    if block is not None:
                        labels.append(block)

                m = len(labels)
                a = len(set(labels))
                b = count_blocks(labels)

                if m > a:
                    fi = (b - a) / (m - a)
                else:
                    fi = 0.0

                out.write(
                    f"{species}\t{chromosome}\t{m}\t{a}\t{b}\t{fi:.4f}\n"
                )


if __name__ == "__main__":
    main()
