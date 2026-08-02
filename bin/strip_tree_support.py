#!/usr/bin/env python3
"""
Strip internal-node branch-support labels (e.g. IQ-TREE's "100/100" combined
SH-aLRT/UFBoot value) from a Newick tree, leaving topology, branch lengths,
and leaf names untouched.

Syngraph's own tree loader (ete3.Tree(), called by `syngraph infer`) parses
an internal-node label strictly as a single plain support value and rejects
IQ-TREE's "<SH-aLRT>/<UFBoot>" combined format outright with a NewickError --
and Syngraph doesn't use support values for anything (only topology and
branch lengths, per its own docs), so the simplest fix is to remove them
rather than pick one of the two numbers to keep.

A label can only appear immediately after a closing ')' in Newick (leaf names
never do), so anchoring on that is enough to identify it without a full
parser.
"""
import argparse
import re

_SUPPORT_VALUE_RE = re.compile(
    r"\)[0-9]+(\.[0-9]+)?(/[0-9]+(\.[0-9]+)?)?(?=[:,);])"
)


def strip_support_values(newick):
    return _SUPPORT_VALUE_RE.sub(")", newick)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input Newick tree")
    parser.add_argument("-o", "--output", required=True, help="Output Newick tree")
    args = parser.parse_args()

    with open(args.input) as fh:
        newick = fh.read().strip()

    with open(args.output, "w") as out:
        out.write(strip_support_values(newick) + "\n")


if __name__ == "__main__":
    main()
