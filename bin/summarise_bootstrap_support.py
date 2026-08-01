#!/usr/bin/env python3
"""
Estimate bootstrap support for Syngraph's ALG (ancestral linkage group)
calls: for each reconstructed ancestral node and each marker Syngraph placed
into an ALG on the full (reference) dataset, what fraction of bootstrap
replicates -- each run on a marker set resampled by resample_markers.py --
place that marker back with the same group of markers.

ALG labels (e.g. "n5_2") are assigned independently in each replicate run and
are not comparable by name across runs, so for each node we first match each
replicate's ALGs to the reference run's ALGs by maximum marker overlap
(greedy best match), then check, per marker, whether the replicate's
re-mapped ALG call agrees with the reference call. A replicate only "votes"
on markers it actually retained after resampling.

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
import sys
from collections import Counter, defaultdict


def parse_table(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    if header[0] != "marker":
        sys.exit(f"ERROR: unexpected first column '{header[0]}' in {path}")

    taxa = []
    col_index = {}
    i = 1
    while i < len(header):
        name = header[i]
        if not name.endswith("_seq"):
            sys.exit(f"ERROR: unexpected column '{name}' in {path}")
        taxon = name[: -len("_seq")]
        taxa.append(taxon)
        col_index[taxon] = i
        i += 3

    return taxa, col_index, rows


def is_ancestral(taxon, col_index, rows):
    start_col = col_index[taxon] + 1
    seen_any = False
    for row in rows:
        if row[col_index[taxon]] == "NA":
            continue
        seen_any = True
        if row[start_col] != "NA":
            return False
    return seen_any


def marker_alg_by_node(path):
    """Return {ancestral_node: {marker: alg}}."""
    taxa, col_index, rows = parse_table(path)
    ancestral_taxa = [t for t in taxa if is_ancestral(t, col_index, rows)]

    result = {}
    for anc in ancestral_taxa:
        anc_col = col_index[anc]
        marker_alg = {}
        for row in rows:
            alg = row[anc_col]
            if alg != "NA":
                marker_alg[row[0]] = alg
        result[anc] = marker_alg
    return result


def match_algs(ref_marker_alg, rep_marker_alg):
    """
    Greedily map each replicate ALG to the reference ALG it shares the most
    markers with, among markers present in both. Returns {rep_alg: ref_alg}.
    """
    overlap = defaultdict(Counter)
    for marker, rep_alg in rep_marker_alg.items():
        ref_alg = ref_marker_alg.get(marker)
        if ref_alg is not None:
            overlap[rep_alg][ref_alg] += 1

    return {
        rep_alg: counts.most_common(1)[0][0]
        for rep_alg, counts in overlap.items()
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, help="Full-dataset algo.table.tsv")
    parser.add_argument("--output", required=True)
    parser.add_argument("replicates", nargs="+", help="Per-replicate algo.table.tsv files")
    args = parser.parse_args()

    ref_by_node = marker_alg_by_node(args.reference)
    if not ref_by_node:
        sys.exit(f"ERROR: no reconstructed ancestral node columns found in {args.reference}")

    present = defaultdict(Counter)    # (node) -> marker -> n replicates the marker appears in
    supported = defaultdict(Counter)  # (node) -> marker -> n of those where the call agreed

    for replicate_path in args.replicates:
        rep_by_node = marker_alg_by_node(replicate_path)
        for node, ref_marker_alg in ref_by_node.items():
            rep_marker_alg = rep_by_node.get(node, {})
            rep_to_ref = match_algs(ref_marker_alg, rep_marker_alg)

            for marker, rep_alg in rep_marker_alg.items():
                ref_alg = ref_marker_alg.get(marker)
                if ref_alg is None:
                    continue
                present[node][marker] += 1
                if rep_to_ref.get(rep_alg) == ref_alg:
                    supported[node][marker] += 1

    with open(args.output, "w") as out:
        out.write("ancestral_node\tmarker\talg\tn_replicates\tn_supported\tsupport\n")
        for node, ref_marker_alg in ref_by_node.items():
            for marker, alg in sorted(ref_marker_alg.items()):
                n = present[node][marker]
                s = supported[node][marker]
                support = f"{s / n:.4f}" if n else "NA"
                out.write(f"{node}\t{marker}\t{alg}\t{n}\t{s}\t{support}\n")


if __name__ == "__main__":
    main()
