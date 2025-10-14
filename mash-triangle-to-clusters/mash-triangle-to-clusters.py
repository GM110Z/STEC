#!/usr/bin/env python3
import sys, os, gzip
from collections import defaultdict, deque

USAGE = f"""\
Usage:
  {sys.argv[0]} <triangle.txt[.gz]> <threshold> <out.tsv> [members.txt]

Reads raw `mash triangle` output (lower-triangular matrix), unions pairs with distance <= threshold,
and writes connected components. Optional members.txt (one ID per line) adds singletons.
"""

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

class DSU:
    __slots__ = ("parent","rank")
    def __init__(self):
        self.parent = []
        self.rank = []
    def add(self):
        i = len(self.parent)
        self.parent.append(i)
        self.rank.append(0)
        return i
    def find(self, x):
        p = self.parent
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb: return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1

def main():
    if len(sys.argv) < 4:
        sys.stderr.write(USAGE + "\n"); sys.exit(1)
    tri_path  = sys.argv[1]
    thresh    = float(sys.argv[2])
    out_path  = sys.argv[3]
    manifest  = sys.argv[4] if len(sys.argv) > 4 else None

    dsu = DSU()
    names = []                 # index -> name
    name2id = {}               # name -> index

    # stream the ragged triangle; union edges <= threshold
    with open_maybe_gz(tri_path) as fh:
        first = True
        row = -1
        for line in fh:
            if first:
                first = False   # skip the count line
                continue
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            cur_name = parts[0]
            row += 1  # 0-based index for this row
            if cur_name in name2id:
                cur_id = name2id[cur_name]
            else:
                cur_id = dsu.add()
                name2id[cur_name] = cur_id
                names.append(cur_name)

            # distances to previous rows: fields 2..NF map to ids 0..row-1
            # field j (2-based) -> previous index (j-2)
            for j in range(1, len(parts)):
                dj = parts[j]
                try:
                    d = float(dj)
                except ValueError:
                    continue
                if d <= thresh:
                    prev_id = j - 1  # because j starts at 1 for parts[1]
                    dsu.union(prev_id, cur_id)

            # optional light progress
            if (row+1) % 50000 == 0:
                sys.stderr.write(f"[INFO] processed rows: {row+1}\n")

    # include singletons from manifest (names not seen in triangle)
    if manifest and os.path.exists(manifest):
        with open(manifest) as mf:
            for line in mf:
                s = line.strip()
                if not s: continue
                if s not in name2id:
                    nid = dsu.add()
                    name2id[s] = nid
                    names.append(s)

    # assign cluster ids by DSU roots
    root2cid = {}
    cid = 0
    with open(out_path, "w") as out:
        out.write("genome\tcluster_id\n")
        for idx, nm in enumerate(names):
            r = dsu.find(idx)
            if r not in root2cid:
                cid += 1
                root2cid[r] = cid
            out.write(f"{nm}\t{root2cid[r]}\n")

    sys.stderr.write(f"[INFO] Threshold: {thresh}\n")
    sys.stderr.write(f"[INFO] Genomes written: {len(names):,}\n")
    sys.stderr.write(f"[INFO] Clusters: {len(root2cid):,}\n")

if __name__ == "__main__":
    main()
