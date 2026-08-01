#!/usr/bin/env python3
"""
Draw one bootstrap resample of BUSCO markers from a set of per-species
filtered BUSCO tables (busco_filter.py's 5-column output), then rewrite each
species' table restricted to the resampled marker set.

Markers are drawn with replacement from the union of marker (BUSCO) IDs
across all input species, up to the size of that union, then deduplicated to
a set -- a Syngraph marker is a single graph node and can't meaningfully
appear twice, so this is the standard translation of "resample with
replacement" (the same style used for gene-resampling bootstraps over
concatenated phylogenetic loci) to a graph/clustering setting.

This is an original implementation, not derived from any other codebase, but
the underlying idea -- resampling markers to put a confidence value on each
Syngraph ALG call -- follows the `boot10k` bootstrap robustness check in
Maulana et al. 2026 (bioRxiv 2026.07.17.739156), whose companion repository
is https://github.com/Obscuromics/coleoptera-ALGs (`scripts/FigS8.plot.syngraph.boot10k.R`,
not consulted here -- Syngraph itself has no built-in bootstrap flag, so both
that script and this one reimplement the resampling independently). See
"Citation" in README.md.
"""
import argparse
import os
import random
import sys


def read_marker_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            ids.append(line.split("\t", 1)[0])
    return ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("tables", nargs="+", help="Per-species *.filtered.tsv files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    all_markers = sorted({
        marker
        for table in args.tables
        for marker in read_marker_ids(table)
    })
    if not all_markers:
        sys.exit("ERROR: no markers found across input tables")

    rng = random.Random(args.seed)
    sampled = {rng.choice(all_markers) for _ in range(len(all_markers))}

    for table in args.tables:
        out_path = os.path.join(args.outdir, os.path.basename(table))
        with open(table) as fh, open(out_path, "w") as out:
            for line in fh:
                if line.split("\t", 1)[0] in sampled:
                    out.write(line)


if __name__ == "__main__":
    main()
