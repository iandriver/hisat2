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

"""EmptyDrops_CR cell calling for a Solo.out raw matrix.

This lives in Python rather than in the aligner because it needs an ambient
profile estimate, a Monte-Carlo multinomial simulation and a Benjamini-Hochberg
correction -- a pile of statistics that is far easier to get right, and to
check, with numpy than in the aligner's C++.

It is a good-faith implementation of the CellRanger 3 / EmptyDrops approach,
not a bit-exact reimplementation. Cells called here will not match CellRanger's
list exactly; the knee filter that hisat2 already applied is usually within a
fraction of a percent, and this mainly rescues small cells that sit below it.
"""

import argparse
import os
import sys

try:
    import numpy as np
except ImportError:
    sys.stderr.write(
        "error: EmptyDrops_CR needs numpy. Install it, or use\n"
        "       --solo-cell-filter CellRanger2.2 for the knee filter alone.\n")
    sys.exit(2)


def read_mtx(path):
    """Reads a MatrixMarket coordinate file as (rows, cols, gene, cell, val)."""
    with open(path) as f:
        f.readline()
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        ngene, ncell, nnz = (int(x) for x in line.split())
        g = np.empty(nnz, dtype=np.int64)
        c = np.empty(nnz, dtype=np.int64)
        v = np.empty(nnz, dtype=np.int64)
        for i, line in enumerate(f):
            a = line.split()
            g[i] = int(a[0]) - 1
            c[i] = int(a[1]) - 1
            v[i] = int(a[2])
    return ngene, ncell, g, c, v


def simulate_pvalues(ambient, totals, obs_ll, niters, seed):
    """Monte-Carlo p-values: how often a multinomial draw from the ambient
    profile is at least as likely as what was observed.

    Barcodes are grouped by total count so one batch of simulations serves
    every barcode of that size, which is what makes this tractable.
    """
    rng = np.random.default_rng(seed)
    logp = np.log(ambient)
    pvals = np.ones(len(totals))
    order = np.argsort(totals)
    uniq, starts = np.unique(totals[order], return_index=True)
    starts = list(starts) + [len(order)]

    for k, t in enumerate(uniq):
        idx = order[starts[k]:starts[k + 1]]
        if t <= 0:
            continue
        # Likelihood of niters random molecules-from-ambient draws of this size.
        draws = rng.multinomial(int(t), ambient, size=niters)
        sim_ll = draws @ logp
        sim_ll.sort()
        # p = P(sim <= observed): a barcode is "ambient-like" when its
        # likelihood under ambient is high, so low p means it deviates.
        pos = np.searchsorted(sim_ll, obs_ll[idx], side='right')
        pvals[idx] = (pos + 1) / (niters + 1)
    return pvals


def bh_adjust(p):
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def main():
    ap = argparse.ArgumentParser(description='EmptyDrops_CR cell calling for Solo.out')
    ap.add_argument('--raw', required=True, help='Solo.out/<feature>/raw directory')
    ap.add_argument('--out', help='output directory (default: ../filtered)')
    ap.add_argument('--expect-cells', type=int, default=3000)
    ap.add_argument('--max-percentile', type=float, default=0.99)
    ap.add_argument('--max-min-ratio', type=int, default=10)
    ap.add_argument('--ambient-max', type=int, default=100,
                    help='barcodes with at most this many UMIs define the ambient profile')
    ap.add_argument('--candidate-min', type=int, default=0,
                    help='candidates are tested down to this UMI count '
                         '(default: just above --ambient-max). Must be below the '
                         'knee threshold, or there is nothing left to rescue.')
    ap.add_argument('--fdr', type=float, default=0.01)
    ap.add_argument('--niters', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=0, help='0 gives reproducible output')
    args = ap.parse_args()

    raw = args.raw.rstrip('/')
    out = args.out or os.path.join(os.path.dirname(raw), 'filtered')

    ngene, ncell, g, c, v = read_mtx(os.path.join(raw, 'matrix.mtx'))
    barcodes = [l.rstrip('\n') for l in open(os.path.join(raw, 'barcodes.tsv'))]
    totals = np.zeros(ncell, dtype=np.int64)
    np.add.at(totals, c, v)
    nonzero = np.flatnonzero(totals > 0)
    sys.stderr.write("barcodes with UMIs: %d\n" % len(nonzero))

    # 1. Knee filter: these are cells regardless of what the test says.
    srt = np.sort(totals[nonzero])[::-1]
    idx = min(int(args.expect_cells * (1 - args.max_percentile)), len(srt) - 1)
    umi_max = srt[idx]
    knee = max(1, umi_max // args.max_min_ratio)
    is_cell = totals >= knee
    sys.stderr.write("knee threshold: %d UMIs -> %d cells\n" % (knee, int(is_cell.sum())))

    # 2. Ambient profile from the many tiny barcodes, with a pseudocount so no
    #    gene has zero probability (an observed count there would otherwise give
    #    a likelihood of -inf).
    amb_mask = (totals > 0) & (totals <= args.ambient_max)
    amb = np.zeros(ngene, dtype=np.float64)
    sel = amb_mask[c]
    np.add.at(amb, g[sel], v[sel])
    if amb.sum() <= 0:
        sys.stderr.write("warning: no ambient barcodes found; keeping the knee call\n")
        called = np.flatnonzero(is_cell)
    else:
        amb += 1.0
        amb /= amb.sum()

        # 3. Candidates: below the knee but not obviously empty.
        # Candidates sit between the ambient ceiling and the knee: barcodes
        # too big to be pure ambient, too small to have been called outright.
        cand_min = args.candidate_min if args.candidate_min > 0 else args.ambient_max + 1
        if cand_min >= knee:
            sys.stderr.write(
                "warning: candidate floor (%d) is at or above the knee (%d), so no "
                "barcode can be rescued; lower --candidate-min or --ambient-max\n"
                % (cand_min, knee))
        cand = np.flatnonzero((~is_cell) & (totals >= cand_min))
        sys.stderr.write("ambient barcodes: %d, candidates: %d\n"
                         % (int(amb_mask.sum()), len(cand)))
        if len(cand) == 0:
            called = np.flatnonzero(is_cell)
        else:
            pos = {b: i for i, b in enumerate(cand)}
            obs_ll = np.zeros(len(cand))
            logamb = np.log(amb)
            keep = np.isin(c, cand)
            for gi, ci, vi in zip(g[keep], c[keep], v[keep]):
                obs_ll[pos[ci]] += vi * logamb[gi]

            p = simulate_pvalues(amb, totals[cand], obs_ll, args.niters, args.seed)
            q = bh_adjust(p)
            rescued = cand[q <= args.fdr]
            sys.stderr.write("rescued by EmptyDrops: %d\n" % len(rescued))
            is_cell[rescued] = True
            called = np.flatnonzero(is_cell)

    # 4. Write the filtered matrix.
    os.makedirs(out, exist_ok=True)
    colof = -np.ones(ncell, dtype=np.int64)
    colof[called] = np.arange(len(called))
    keep = colof[c] >= 0
    gk, ck, vk = g[keep], colof[c[keep]], v[keep]
    order = np.lexsort((gk, ck))

    with open(os.path.join(out, 'barcodes.tsv'), 'w') as o:
        for i in called:
            o.write(barcodes[i] + "\n")
    src_feat = os.path.join(raw, 'features.tsv')
    with open(src_feat) as f, open(os.path.join(out, 'features.tsv'), 'w') as o:
        o.write(f.read())
    with open(os.path.join(out, 'matrix.mtx'), 'w') as o:
        o.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        o.write("%d %d %d\n" % (ngene, len(called), len(order)))
        for i in order:
            o.write("%d %d %d\n" % (gk[i] + 1, ck[i] + 1, vk[i]))

    sys.stderr.write("cells called: %d (written to %s)\n" % (len(called), out))


if __name__ == '__main__':
    main()
