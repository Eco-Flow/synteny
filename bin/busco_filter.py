#!/usr/bin/env python3
import argparse
import sys


def read_exclusions(path, species):
    excluded = set()
    if path is None:
        return excluded
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            row_species, scaffold = fields[0], fields[1]
            if row_species == species:
                excluded.add(scaffold)
    return excluded


def main():
    parser = argparse.ArgumentParser(
        description="Filter a BUSCO full_table.tsv down to single-copy (Complete) "
        "markers and reformat to a 5-column TSV (busco_id, sequence, start, end, "
        "strand) shared by the Syngraph and AGORA input-prep steps."
    )
    parser.add_argument("species", help="Species/sample identifier")
    parser.add_argument("full_table", help="Path to BUSCO full_table.tsv")
    parser.add_argument("output", help="Path to write the filtered TSV")
    parser.add_argument(
        "--exclude-scaffolds",
        default=None,
        help="Optional TSV of species<TAB>scaffold_name rows to drop "
        "(e.g. unlocalised scaffolds or a Y chromosome)",
    )
    args = parser.parse_args()

    excluded = read_exclusions(args.exclude_scaffolds, args.species)

    kept = 0
    dropped_status = 0
    dropped_excluded = 0

    with open(args.full_table) as fh, open(args.output, "w") as out:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            busco_id, status, sequence, start, end, strand = fields[:6]
            if status != "Complete":
                dropped_status += 1
                continue
            if sequence in excluded:
                dropped_excluded += 1
                continue
            out.write(f"{busco_id}\t{sequence}\t{start}\t{end}\t{strand}\n")
            kept += 1

    sys.stderr.write(
        f"{args.species}: kept {kept} single-copy markers "
        f"(dropped {dropped_status} not Complete, {dropped_excluded} excluded scaffolds)\n"
    )


if __name__ == "__main__":
    main()
