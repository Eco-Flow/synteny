#!/usr/bin/env python3
"""
Summarise Syngraph's per-marker table (algo.table.tsv, from `syngraph tabulate`)
into two things Syngraph itself doesn't report directly:

1. The total number of ancestral linkage groups (ALGs) at each reconstructed
   ancestral node.
2. Per species, per ALG: whether that ALG's markers all still sit on one
   chromosome ("intact"), have split across multiple chromosomes ("split"),
   and/or share a chromosome with another ALG's markers ("fused"). An ALG can
   be both split and fused at once (e.g. one piece stayed put, the other
   fused into a different chromosome), so these are independent flags, not a
   single exclusive category.

algo.table.tsv's columns are #marker, then three columns per taxon (extant
species and reconstructed ancestral nodes alike): <taxon>_seq, <taxon>_start,
<taxon>_end. Ancestral nodes are identified by having 'NA' start/end
coordinates throughout (they have no real genomic coordinates, just a
chromosome/linkage-group label), rather than by name, since taxon names are
arbitrary.
"""
import argparse
import sys
from collections import defaultdict


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, help="Syngraph algo.table.tsv")
    parser.add_argument(
        "--alg-summary", required=True,
        help="Output: total ALG count per reconstructed ancestral node",
    )
    parser.add_argument(
        "--status-summary", required=True,
        help="Output: per-species, per-ALG intact/split/fused status",
    )
    args = parser.parse_args()

    taxa, col_index, rows = parse_table(args.table)
    ancestral_taxa = [t for t in taxa if is_ancestral(t, col_index, rows)]
    extant_taxa = [t for t in taxa if t not in ancestral_taxa]

    if not ancestral_taxa:
        sys.exit(f"ERROR: no reconstructed ancestral node columns found in {args.table}")

    with open(args.alg_summary, "w") as alg_out, open(args.status_summary, "w") as status_out:
        alg_out.write("ancestral_node\tn_algs\tn_markers\n")
        status_out.write(
            "ancestral_node\talg\tspecies\tstatus\tchromosomes\tn_markers\tshared_with_algs\n"
        )

        for anc in ancestral_taxa:
            anc_col = col_index[anc]
            # marker -> ALG label at this ancestral node (markers absent at
            # this node, i.e. not shared this far back, are skipped)
            marker_alg = {}
            for row in rows:
                alg = row[anc_col]
                if alg != "NA":
                    marker_alg[row[0]] = alg

            algs = sorted(set(marker_alg.values()))
            alg_out.write(f"{anc}\t{len(algs)}\t{len(marker_alg)}\n")

            for species in extant_taxa:
                sp_col = col_index[species]
                chrom_algs = defaultdict(set)   # chromosome -> ALGs found on it
                alg_chroms = defaultdict(set)   # ALG -> chromosomes its markers are on
                alg_marker_count = defaultdict(int)

                for row in rows:
                    marker = row[0]
                    alg = marker_alg.get(marker)
                    chrom = row[sp_col]
                    if alg is None or chrom == "NA":
                        continue
                    chrom_algs[chrom].add(alg)
                    alg_chroms[alg].add(chrom)
                    alg_marker_count[alg] += 1

                for alg in algs:
                    chroms = alg_chroms.get(alg, set())
                    if not chroms:
                        status_out.write(f"{anc}\t{alg}\t{species}\tabsent\t\t0\t\n")
                        continue

                    shared_with = set()
                    for chrom in chroms:
                        shared_with |= (chrom_algs[chrom] - {alg})

                    statuses = []
                    if len(chroms) > 1:
                        statuses.append("split")
                    if shared_with:
                        statuses.append("fused")
                    if not statuses:
                        statuses.append("intact")

                    status_out.write(
                        f"{anc}\t{alg}\t{species}\t{'+'.join(statuses)}\t"
                        f"{','.join(sorted(chroms))}\t{alg_marker_count[alg]}\t"
                        f"{','.join(sorted(shared_with))}\n"
                    )


if __name__ == "__main__":
    main()
