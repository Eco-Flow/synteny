#!/usr/bin/env python3
"""
Cap the number of single-copy orthogroups carried into the species-tree
supermatrix (--max_orthogroups), keeping the longest ones.

Strictly single-copy orthogroups become far more numerous for closely
related species sets than for distantly related ones (the "every species
has exactly one copy" filter is easier to satisfy the more similar the
species are), so an uncapped run can end up with a wildly different number
of IQ-TREE partitions -- and runtime, since -m MFP model selection runs
per partition -- depending only on how divergent the input species happen
to be, not on any deliberate choice. This caps that at a fixed number,
ranked by orthogroup length (the standard, simple proxy for phylogenetic
informativeness: more sites, more signal) rather than by a conservation
score, since a conservation score would need the alignment/tree this step
precedes in order to compute.
"""
import argparse
import os
import shutil
import sys


def mean_sequence_length(path):
    """Mean sequence length across the orthogroup's single-copy records."""
    lengths = []
    length = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if length:
                    lengths.append(length)
                length = 0
            else:
                length += len(line.strip())
        if length:
            lengths.append(length)
    return sum(lengths) / len(lengths) if lengths else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max", type=int, required=True)
    parser.add_argument("--indir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    orthogroups = sorted(f for f in os.listdir(args.indir) if f.endswith(".fa"))
    if not orthogroups:
        sys.exit(f"ERROR: no orthogroup FASTA files found in {args.indir}")

    ranked = sorted(
        orthogroups,
        key=lambda f: (-mean_sequence_length(os.path.join(args.indir, f)), f),
    )
    kept = ranked[: args.max]

    os.makedirs(args.outdir, exist_ok=True)
    for f in kept:
        shutil.copy(os.path.join(args.indir, f), os.path.join(args.outdir, f))

    sys.stderr.write(
        f"Selected {len(kept)} of {len(orthogroups)} single-copy orthogroups "
        f"(--max_orthogroups {args.max}), ranked by mean sequence length\n"
    )


if __name__ == "__main__":
    main()
