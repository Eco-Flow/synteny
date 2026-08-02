#!/usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict


def load_filtered_table(path):
    """Read a species' filtered BUSCO TSV: busco_id, sequence, start, end, strand."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            busco_id, sequence, start, end, strand = line.split("\t")
            rows.append((busco_id, sequence, int(start), int(end), strand))
    return rows


def strand_to_agora(strand):
    return "1" if strand == "+" else "-1"


_SUPPORT_VALUE_RE = re.compile(r"^[0-9]+(\.[0-9]+)?(/[0-9]+(\.[0-9]+)?)?$")


def label_internal_nodes(newick):
    """AGORA requires unique names on every internal node. Auto-name any that
    are missing OR that carry only a branch-support value (e.g. IQ-TREE's
    "97.3/100" SH-aLRT/UFBoot label, as re-emitted by root_tree.py) -- a
    support value labels a bipartition, not a node, so AGORA's own tree
    loader treats it the same as no name and invents its own placeholder
    ("NAME_0", "NAME_1", ...) that this script would otherwise never learn
    about, leaving no matching orthologyGroups.NAME_0.list file for it to
    read. Genuinely-named nodes (a real ancestor name already present, e.g.
    from a user-supplied --species_tree, or an AlgoAncN name from a prior
    run of this function) are left untouched, as are branch lengths."""
    counter = [0]

    def next_name():
        counter[0] += 1
        return f"AlgoAnc{counter[0]}"

    out = []
    i = 0
    n = len(newick)
    while i < n:
        c = newick[i]
        if c == ")":
            out.append(c)
            i += 1
            # Look ahead: is there already a label (name/support and/or branch
            # length) immediately after this closing paren?
            j = i
            label_chars = []
            while j < n and newick[j] not in ",);":
                label_chars.append(newick[j])
                j += 1
            label = "".join(label_chars)
            name_part = re.match(r"^[^:]*", label).group(0)
            rest = label[len(name_part):]  # ':branch_length', if present
            if name_part.strip() == "" or _SUPPORT_VALUE_RE.match(name_part.strip()):
                out.append(next_name())
                out.append(rest)
            else:
                out.append(label)
            i = j
            continue
        out.append(c)
        i += 1
    return "".join(out)


class Node:
    __slots__ = ("name", "children", "leaves")

    def __init__(self, name=""):
        self.name = name
        self.children = []
        self.leaves = None  # set of leaf names in this node's subtree


def parse_newick(text):
    """Parse a (already internally-labelled) Newick string into a Node tree,
    ignoring branch lengths -- only topology and names are needed to compute
    MRCAs."""
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    pos = [0]

    def parse_node():
        node = Node()
        if pos[0] < len(text) and text[pos[0]] == "(":
            pos[0] += 1
            while True:
                node.children.append(parse_node())
                if text[pos[0]] == ",":
                    pos[0] += 1
                    continue
                if text[pos[0]] == ")":
                    pos[0] += 1
                    break
        start = pos[0]
        while pos[0] < len(text) and text[pos[0]] not in ",():":
            pos[0] += 1
        label = text[start:pos[0]]
        node.name = label.split(":")[0]
        if pos[0] < len(text) and text[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(text) and text[pos[0]] not in ",()":
                pos[0] += 1
        return node

    return parse_node()


def annotate_leaves(node):
    """Post-order: fill in node.leaves (the set of leaf names under it)."""
    if not node.children:
        node.leaves = {node.name}
        return node.leaves
    leaves = set()
    for child in node.children:
        leaves |= annotate_leaves(child)
    node.leaves = leaves
    return leaves


def internal_nodes(node):
    """All internal (non-leaf) nodes in this subtree, in no particular order."""
    if not node.children:
        return []
    result = [node]
    for child in node.children:
        result.extend(internal_nodes(child))
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Convert filtered per-species BUSCO tables into AGORA's "
        "gene-list / orthology-group / species-tree input formats."
    )
    parser.add_argument(
        "--tables",
        nargs="+",
        required=True,
        help="species=path pairs, e.g. Species_a=Species_a.filtered.tsv",
    )
    parser.add_argument("--tree", required=True, help="Newick species tree")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genes_dir = os.path.join(args.outdir, "genes")
    os.makedirs(genes_dir, exist_ok=True)

    species_tables = {}
    for pair in args.tables:
        species, path = pair.split("=", 1)
        species_tables[species] = load_filtered_table(path)

    # orthology groups: one BUSCO id -> list of (species, "<species>_<busco_id>") pairs.
    # Species are tracked alongside each gene id, rather than recovered later by
    # splitting the id string, since species names can themselves contain
    # underscores (e.g. "Drosophila_yakuba") and would otherwise be truncated.
    groups = defaultdict(list)
    for species, rows in species_tables.items():
        gene_list_path = os.path.join(genes_dir, f"genes.{species}.list")
        with open(gene_list_path, "w") as out:
            for busco_id, sequence, start, end, strand in rows:
                gene_id = f"{species}_{busco_id}"
                out.write(
                    f"{sequence}\t{start}\t{end}\t{strand_to_agora(strand)}\t{gene_id}\n"
                )
                groups[busco_id].append((species, gene_id))

    with open(args.tree) as fh:
        newick = fh.read().strip()
    labelled_tree = label_internal_nodes(newick)
    tree_path = os.path.join(args.outdir, "species_tree.nwk")
    with open(tree_path, "w") as out:
        out.write(labelled_tree + "\n")

    # AGORA's orthologyGroups input is not one flat file but one file per internal
    # (ancestor) node of the species tree, following the same "%s" per-node
    # convention as its gene-list files (doc/HowTo.md: "the list of orthology
    # groups present on each (internal) node of the species tree" -- confirmed
    # against a real AGORA run, whose ALL.reformatGeneFamilies.py step processes
    # every internal node in turn). Since these are single-copy, duplication-free
    # BUSCO families, each orthogroup is written -- restricted to the genes of the
    # species actually descended from that node -- at *every* ancestor whose
    # descendant species overlap it by 2 or more (a node with only 0-1 of a
    # group's species has nothing to reconcile there).
    root = parse_newick(labelled_tree)
    annotate_leaves(root)

    single_copy_groups = {
        busco_id: pairs for busco_id, pairs in groups.items() if len(pairs) >= 2
    }

    orthology_dir = os.path.join(args.outdir, "orthologyGroups")
    os.makedirs(orthology_dir, exist_ok=True)
    for node in internal_nodes(root):
        lines = []
        for busco_id, pairs in single_copy_groups.items():
            restricted = [gene_id for species, gene_id in pairs if species in node.leaves]
            if len(restricted) >= 2:
                lines.append(" ".join(restricted))
        with open(os.path.join(orthology_dir, f"orthologyGroups.{node.name}.list"), "w") as out:
            out.write("\n".join(lines) + ("\n" if lines else ""))


if __name__ == "__main__":
    main()
