#!/usr/bin/env python3

#
# Copyright 2026, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

"""Compare two single-cell count matrices (e.g. HISAT2 Solo.out vs STARsolo).

Equality is not the target and never will be: the two align to different
indexes with different scoring, so some reads land differently by
construction. What matters is whether the disagreements are *explained*.
This reports concordance and then decomposes the difference into causes, so
a residual can be attributed rather than merely observed.

Matrices are keyed by (barcode, gene_id) rather than by row/column index,
because feature ordering differs between tools.
"""

import argparse
import collections
import math
import os
import sys


def load(mtxdir):
    feats = [l.split('\t')[0] for l in open(os.path.join(mtxdir, 'features.tsv'))]
    bcs = [l.strip() for l in open(os.path.join(mtxdir, 'barcodes.tsv'))]
    d = {}
    percell = collections.Counter()
    pergene = collections.Counter()
    with open(os.path.join(mtxdir, 'matrix.mtx')) as f:
        f.readline()
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        for line in f:
            a = line.split()
            g = feats[int(a[0]) - 1]
            c = bcs[int(a[1]) - 1]
            v = float(a[2])
            d[(c, g)] = v
            percell[c] += v
            pergene[g] += v
    return d, percell, pergene


def pearson(a, b):
    n = len(a)
    if n < 2:
        return float('nan')
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    return num / (da * db) if da * db else float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--a', required=True, help='first raw/ directory (e.g. HISAT2)')
    ap.add_argument('--b', required=True, help='second raw/ directory (e.g. STARsolo)')
    ap.add_argument('--a-name', default='A')
    ap.add_argument('--b-name', default='B')
    ap.add_argument('--cells-a', help='optional filtered/barcodes.tsv for A')
    ap.add_argument('--cells-b', help='optional filtered/barcodes.tsv for B')
    args = ap.parse_args()

    A, Ac, Ag = load(args.a)
    B, Bc, Bg = load(args.b)
    na, nb = args.a_name, args.b_name

    print("=" * 68)
    print("  %s vs %s" % (na, nb))
    print("=" * 68)
    print("\n-- matrix totals --")
    print("%-22s %12s %12s" % ("", na, nb))
    print("%-22s %12d %12d" % ("non-zero entries", len(A), len(B)))
    print("%-22s %12.0f %12.0f" % ("total UMIs", sum(A.values()), sum(B.values())))
    print("%-22s %12d %12d" % ("barcodes with UMIs", len(Ac), len(Bc)))
    print("%-22s %12d %12d" % ("genes detected", len(Ag), len(Bg)))
    depth = sum(A.values()) / max(sum(B.values()), 1)
    print("depth ratio %s/%s: %.4f" % (na, nb, depth))

    keys = set(A) | set(B)
    both = set(A) & set(B)
    print("\n-- (cell, gene) pair agreement --")
    print("shared %d   %s-only %d   %s-only %d   Jaccard %.5f"
          % (len(both), na, len(A) - len(both), nb, len(B) - len(both),
             len(both) / len(keys)))

    va = [A.get(k, 0) for k in keys]
    vb = [B.get(k, 0) for k in keys]
    print("Pearson r (raw)   : %.5f" % pearson(va, vb))
    print("Pearson r (log1p) : %.5f" % pearson([math.log1p(x) for x in va],
                                               [math.log1p(x) for x in vb]))

    # The raw matrix is dominated by near-empty droplets holding one or two
    # UMIs, where agreement is noise-limited rather than informative. Stratify.
    print("\n-- agreement stratified by expression (max of the two) --")
    print("%-14s %10s %10s %10s" % ("count range", "pairs", "r(raw)", "r(log1p)"))
    for lo, hi, label in ((1, 1, "= 1"), (2, 3, "2-3"), (4, 9, "4-9"),
                          (10, 10 ** 9, ">= 10")):
        sel = [k for k in keys if lo <= max(A.get(k, 0), B.get(k, 0)) <= hi]
        if len(sel) < 2:
            continue
        xa = [A.get(k, 0) for k in sel]
        xb = [B.get(k, 0) for k in sel]
        print("%-14s %10d %10.5f %10.5f"
              % (label, len(sel), pearson(xa, xb),
                 pearson([math.log1p(x) for x in xa], [math.log1p(x) for x in xb])))

    print("\n-- per-cell and per-gene totals --")
    cells = set(Ac) | set(Bc)
    print("per-cell total UMI  r = %.5f  (%d barcodes)"
          % (pearson([Ac.get(c, 0) for c in cells], [Bc.get(c, 0) for c in cells]), len(cells)))
    genes = set(Ag) | set(Bg)
    print("per-gene total UMI  r = %.5f  (%d genes)"
          % (pearson([Ag.get(g, 0) for g in genes], [Bg.get(g, 0) for g in genes]), len(genes)))

    # Depth-normalised per-gene comparison: if the shortfall were purely a
    # difference in how many reads mapped, every gene would scale by the same
    # factor and the residual would vanish.
    print("\n-- is the difference just sequencing depth? --")
    strong = [g for g in genes if Bg.get(g, 0) >= 100]
    ratios = sorted(Ag.get(g, 0) / Bg[g] for g in strong)
    if ratios:
        med = ratios[len(ratios) // 2]
        q1, q3 = ratios[len(ratios) // 4], ratios[3 * len(ratios) // 4]
        print("per-gene ratio %s/%s over %d genes with >=100 UMIs in %s:"
              % (na, nb, len(strong), nb))
        print("   median %.4f   IQR %.4f-%.4f   (overall depth ratio %.4f)"
              % (med, q1, q3, depth))
        within = sum(1 for r in ratios if abs(r - depth) <= 0.10 * depth)
        print("   %d/%d genes (%.1f%%) within 10%% of the overall depth ratio"
              % (within, len(ratios), 100.0 * within / len(ratios)))
        outliers = sorted(((Ag.get(g, 0) / Bg[g], g) for g in strong))
        print("   most %s-depleted: %s" % (na, ", ".join("%s(%.2f)" % (g, r) for r, g in outliers[:3])))
        print("   most %s-enriched: %s" % (na, ", ".join("%s(%.2f)" % (g, r) for r, g in outliers[-3:])))

    if args.cells_a and args.cells_b:
        ca = set(l.strip() for l in open(args.cells_a))
        cb = set(l.strip() for l in open(args.cells_b))
        print("\n-- called cells --")
        print("%s %d   %s %d   shared %d   Jaccard %.5f"
              % (na, len(ca), nb, len(cb), len(ca & cb), len(ca & cb) / len(ca | cb)))
        # Restricted to real cells, where the biology is.
        sel = [k for k in keys if k[0] in (ca & cb)]
        xa = [A.get(k, 0) for k in sel]
        xb = [B.get(k, 0) for k in sel]
        print("within shared cells: %d pairs, r(raw) %.5f, r(log1p) %.5f"
              % (len(sel), pearson(xa, xb),
                 pearson([math.log1p(x) for x in xa], [math.log1p(x) for x in xb])))
        onlyb = sorted((Bc[c] for c in (cb - ca)), reverse=True)
        onlya = sorted((Ac[c] for c in (ca - cb)), reverse=True)
        med = sorted(Ac[c] for c in (ca & cb))
        if med:
            print("median UMIs in shared cells: %d" % med[len(med) // 2])
        if onlya:
            print("%s-only cells: %d, UMI range %d-%d" % (na, len(onlya), onlya[-1], onlya[0]))
        if onlyb:
            print("%s-only cells: %d, UMI range %d-%d" % (nb, len(onlyb), onlyb[-1], onlyb[0]))


if __name__ == '__main__':
    main()
