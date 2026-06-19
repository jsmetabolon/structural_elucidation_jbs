"""Disjoint-set union-find for merging MS1 feature clusters.

Small utility used by Ms1ClusterAnalyzer to combine precursor clusters that
share feature names after RI/mass grouping.
"""

from collections import defaultdict


class UnionFind:
    """Track disjoint sets and merge clusters of related feature identifiers."""

    def __init__(self, elements=[]):
        self.parent = {e: e for e in elements}

    def add_element(self, e):
        if e not in self.parent:
            self.parent[e] = e

    def find(self, f):
        parent = self.parent
        if f not in parent:
            parent[f] = f
            return f
        if f != parent[f]:
            parent[f] = self.find(parent[f])
        return parent[f]

    def union2(self, f1, f2):
        f1 = self.find(f1)
        f2 = self.find(f2)
        if f1 == f2:
            return
        self.parent[f1] = self.parent[f2]

    def in_same_set(self, e1, e2):
        return self.find(e1) == self.find(e2)

    def union(self, elist):
        first, *rest = list(elist)  # might be a set, dict_keys, etc.
        for e in rest:
            self.union2(first, e)

    def get_clusters(self):
        parent = self.parent
        parent = {f: self.find(f) for f in parent.keys()}
        clusters = defaultdict(list)
        for k, v in parent.items():
            clusters[v].append(k)
        return list(clusters.values())
