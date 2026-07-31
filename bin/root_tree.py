#!/usr/bin/env python3

"""
Root an unrooted Newick tree, either on a named outgroup or at the midpoint.

IQ-TREE writes unrooted trees (a trifurcation at the notional root), but CAFE5
requires a rooted binary tree, so the ML species tree has to be rooted before it
can be used downstream.

Internal node labels produced by IQ-TREE are branch support values, which label
*bipartitions* rather than nodes. Rerooting moves nodes relative to each other,
so supports are carried on edges here and re-emitted from each node's incoming
edge afterwards. That keeps every support value attached to the split it was
actually computed for.

Pure standard library (no BioPython/ete3) to match the rest of bin/.
"""

import argparse
import sys


class Node(object):
    __slots__ = ('name', 'length', 'children', 'parent')

    def __init__(self, name='', length=None):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None


def parse_newick(text):
    """Parse a Newick string into a rooted Node tree."""
    text = text.strip()
    if not text.endswith(';'):
        raise ValueError("Newick string does not end with ';'")
    text = text[:-1]

    pos = [0]

    def parse_node():
        node = Node()
        if pos[0] < len(text) and text[pos[0]] == '(':
            pos[0] += 1  # consume '('
            while True:
                node.children.append(parse_node())
                if pos[0] >= len(text):
                    raise ValueError('Unbalanced parentheses in Newick string')
                if text[pos[0]] == ',':
                    pos[0] += 1
                    continue
                if text[pos[0]] == ')':
                    pos[0] += 1
                    break
                raise ValueError("Unexpected character '%s' at position %d"
                                 % (text[pos[0]], pos[0]))
        # Label (leaf name, or internal support value)
        start = pos[0]
        while pos[0] < len(text) and text[pos[0]] not in ',():':
            pos[0] += 1
        node.name = text[start:pos[0]].strip()
        # Branch length
        if pos[0] < len(text) and text[pos[0]] == ':':
            pos[0] += 1
            start = pos[0]
            while pos[0] < len(text) and text[pos[0]] not in ',()':
                pos[0] += 1
            node.length = float(text[start:pos[0]])
        for child in node.children:
            child.parent = node
        return node

    root = parse_node()
    if pos[0] != len(text):
        raise ValueError('Trailing characters after tree at position %d' % pos[0])
    return root


def build_adjacency(root):
    """Convert a rooted tree into an undirected adjacency map.

    Returns (adj, leaf_names) where adj maps node -> list of (neighbour, length,
    edge_label). Internal node labels become edge labels on the branch above the
    node they were attached to.
    """
    adj = {}
    leaves = []

    def walk(node):
        adj.setdefault(node, [])
        if not node.children:
            leaves.append(node)
        for child in node.children:
            length = child.length if child.length is not None else 0.0
            # A child's own label is a support value only if it is internal.
            label = child.name if child.children else ''
            adj.setdefault(child, [])
            adj[node].append((child, length, label))
            adj[child].append((node, length, label))
            walk(child)

    walk(root)

    # Collapse a degree-2 root (an already-rooted input): joining its two
    # neighbours keeps the topology and total path lengths intact.
    if len(adj[root]) == 2:
        (a, la, _), (b, lb, _) = adj[root]
        adj[a] = [e for e in adj[a] if e[0] is not root]
        adj[b] = [e for e in adj[b] if e[0] is not root]
        # The two half-edges belong to the same bipartition; keep either label.
        adj[a].append((b, la + lb, ''))
        adj[b].append((a, la + lb, ''))
        del adj[root]

    return adj, leaves


def farthest_leaf(adj, start):
    """Return (leaf, distance, predecessor_map) for the farthest leaf from start."""
    stack = [(start, None, 0.0)]
    prev = {start: None}
    best = (start, 0.0)
    while stack:
        node, parent, dist = stack.pop()
        if len(adj[node]) == 1 and dist > best[1]:
            best = (node, dist)
        for neighbour, length, _ in adj[node]:
            if neighbour is parent:
                continue
            prev[neighbour] = node
            stack.append((neighbour, node, dist + length))
    return best[0], best[1], prev


def edge_length(adj, u, v):
    for neighbour, length, _ in adj[u]:
        if neighbour is v:
            return length
    raise ValueError('No edge between the given nodes')


def edge_label(adj, u, v):
    for neighbour, _, label in adj[u]:
        if neighbour is v:
            return label
    return ''


def clade_leaves(adj, node, parent):
    """Leaf names on `node`'s side of the (parent, node) edge."""
    names = set()
    stack = [(node, parent)]
    while stack:
        current, came_from = stack.pop()
        if len(adj[current]) == 1:
            names.add(current.name)
        for neighbour, _, _ in adj[current]:
            if neighbour is came_from:
                continue
            stack.append((neighbour, current))
    return names


def find_outgroup_edge(adj, outgroup, all_leaves):
    """Find the edge whose bipartition separates exactly the outgroup taxa."""
    missing = outgroup - all_leaves
    if missing:
        raise ValueError(
            'Outgroup taxa not present in the tree: %s\nTree tips are: %s'
            % (', '.join(sorted(missing)), ', '.join(sorted(all_leaves)))
        )
    if outgroup == all_leaves:
        raise ValueError('The outgroup cannot contain every taxon in the tree')

    for u in adj:
        for v, _, _ in adj[u]:
            side = clade_leaves(adj, v, u)
            if side == outgroup:
                return u, v
    raise ValueError(
        'The outgroup {%s} is not monophyletic in this tree, so it cannot be '
        'used to root it. Either pick a different outgroup or use midpoint '
        'rooting (omit --outgroup).' % ', '.join(sorted(outgroup))
    )


def find_midpoint_edge(adj, leaves):
    """Find the edge containing the midpoint of the longest leaf-to-leaf path."""
    start = leaves[0]
    end_a, _, _ = farthest_leaf(adj, start)
    end_b, diameter, prev = farthest_leaf(adj, end_a)

    # Walk back from end_b towards end_a until half the diameter is covered.
    path = []
    node = end_b
    while node is not None:
        path.append(node)
        node = prev[node]

    half = diameter / 2.0
    travelled = 0.0
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        length = edge_length(adj, u, v)
        if travelled + length >= half:
            # Root sits on edge (u, v), `half - travelled` along it from u.
            return u, v, (half - travelled)
        travelled += length
    # Degenerate case (e.g. all branch lengths zero): root on the first edge.
    return path[0], path[1], 0.0


def reroot(adj, u, v, dist_from_u, label):
    """Insert a new root into the edge (u, v) and return the rooted tree."""
    root = Node(name='')
    total = edge_length(adj, u, v)
    dist_from_u = max(0.0, min(total, dist_from_u))

    def build(node, came_from, length, incoming_label):
        new = Node(name='', length=length)
        neighbours = [e for e in adj[node] if e[0] is not came_from]
        if not neighbours:
            new.name = node.name          # leaf
        else:
            new.name = incoming_label     # support for the edge above
            for neighbour, edge_len, lbl in neighbours:
                child = build(neighbour, node, edge_len, lbl)
                child.parent = new
                new.children.append(child)
        return new

    # Both halves of the split edge describe the same bipartition, so they
    # carry the same support value.
    left = build(u, v, dist_from_u, label)
    right = build(v, u, total - dist_from_u, label)
    left.parent = root
    right.parent = root
    root.children = [left, right]
    return root


def write_newick(node):
    parts = []

    def render(n):
        if n.children:
            parts.append('(')
            for i, child in enumerate(n.children):
                if i:
                    parts.append(',')
                render(child)
            parts.append(')')
        if n.name:
            parts.append(n.name)
        if n.length is not None:
            parts.append(':%g' % n.length)

    render(node)
    parts.append(';')
    return ''.join(parts)


def main():
    parser = argparse.ArgumentParser(
        description='Root an unrooted Newick tree on an outgroup or at its midpoint'
    )
    parser.add_argument('-i', '--input', required=True, help='Input Newick tree')
    parser.add_argument('-o', '--output', required=True, help='Output rooted Newick tree')
    parser.add_argument('-g', '--outgroup', default=None,
                        help='Comma-separated tip name(s) to root on. '
                             'Omit to use midpoint rooting.')

    args = parser.parse_args()

    with open(args.input) as fh:
        text = fh.read().strip()
    if not text:
        sys.exit("ERROR: input tree file '%s' is empty" % args.input)

    try:
        tree = parse_newick(text)
    except ValueError as exc:
        sys.exit('ERROR: could not parse %s: %s' % (args.input, exc))

    adj, leaves = build_adjacency(tree)
    if len(leaves) < 3:
        sys.exit('ERROR: need at least 3 tips to root a tree, found %d' % len(leaves))

    try:
        if args.outgroup:
            outgroup = {t.strip() for t in args.outgroup.split(',') if t.strip()}
            all_leaves = {leaf.name for leaf in leaves}
            u, v = find_outgroup_edge(adj, outgroup, all_leaves)
            total = edge_length(adj, u, v)
            label = edge_label(adj, u, v)
            rooted = reroot(adj, u, v, total / 2.0, label)
            method = 'outgroup (%s)' % ', '.join(sorted(outgroup))
        else:
            u, v, offset = find_midpoint_edge(adj, leaves)
            label = edge_label(adj, u, v)
            rooted = reroot(adj, u, v, offset, label)
            method = 'midpoint'
    except ValueError as exc:
        sys.exit('ERROR: %s' % exc)

    with open(args.output, 'w') as fh:
        fh.write(write_newick(rooted) + '\n')

    print('Tree rooted successfully')
    print('   Method:   %s' % method)
    print('   Tips:     %d' % len(leaves))
    print('   Output:   %s' % args.output)


if __name__ == '__main__':
    main()
