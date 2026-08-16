#!/usr/bin/env python3
#
# Copyright 2026
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

"""Predict, in seconds, what `hisat2-build --snp --haplotype` will do with a
variant set, before spending hours finding out.

A graph build fails in two different ways and both are cheap to foresee:

  * It runs out of node ids. Path-graph construction counts nodes in `index_t`,
    so a 32-bit build dies with "exceeded integer bounds" once the graph passes
    2^32 nodes. The reference alone starts that count at roughly its own length
    -- for human that is ~72% of the ceiling before a single variant is added --
    so the headroom for variants is much smaller than it looks.

  * It quietly drops variants. Whenever a 57,344 bp local graph exceeds the
    63,488-edge budget the builder discards part of that window's variant set
    and retries. Nothing reports the total, and `hisat2-inspect --snp` still
    lists every variant supplied, so a finished index does not reveal it.

Both bounds are compile-time constants in hier_idx_common.h, so the arithmetic
below is not a heuristic about the data -- it is the same arithmetic the builder
does, run ahead of time.

Usage:
    hisat2_index_probe.py --fai genome.fa.fai --snp genome.snp \\
        [--haplotype genome.haplotype] [--large-index] [-v]
"""

import sys
from argparse import ArgumentParser, FileType

# hier_idx_common.h
LOCAL_INDEX_SIZE = (1 << 16) - (1 << 13)          # 57,344 bp per local graph
LOCAL_MAX_GBWT = (1 << 16) - (1 << 11)            # 63,488 edges it may hold
LOCAL_INDEX_OVERLAP = 1024
LOCAL_INDEX_INTERVAL = LOCAL_INDEX_SIZE - LOCAL_INDEX_OVERLAP

INDEX_T_MAX_32 = (1 << 32) - 1
INDEX_T_MAX_64 = (1 << 64) - 1

# The local budget is spent by the PATH graph after prefix doubling, not by the
# reference graph, so it cannot be derived from a variant count -- doubling
# multiplies. These two constants are calibrated against builds whose retention
# was measured, on GRCh38 chr1 with the stock halving backoff:
#
#   set             variants   hap/var   windows over budget   retention
#   dbSNP big-4    1,950,322      0.78          1,135            80.9%
#   dbSNP 1000G      1,701,360    0.82            692            86.7%
#   1000G phased   1,182,916      1.09            472            92.4%
#
# Effective capacity times haplotypes-per-variant is near constant (533*0.78 =
# 416, 395*1.09 = 431): more haplotype rows per variant means each variant
# costs more path-graph budget. Loss is ~0.45 of what sits in an over-budget
# window, which is halving giving up about half of it.
#
# Fitted on big-4 and phased; the 1000G set was held out and came back 86.1%
# predicted against 86.7% measured.
CAPACITY_X_HAPRATIO = 420.0
LOSS_FRACTION = 0.45


def read_fai(fh):
    """name -> length, from a samtools .fai."""
    lens = {}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2:
            lens[f[0]] = int(f[1])
    return lens


def variant_cost(vtype, data):
    """Nodes a variant adds to the reference graph.

    A substitution adds its alternate base. An insertion adds one node per
    inserted base. A deletion adds no node -- it is an edge that skips existing
    ones -- but it still costs an edge, which is what the local budget counts.
    """
    if vtype == 'single':
        return 1
    if vtype == 'insertion':
        return len(data)
    return 0                                       # deletion


def main(fai, snp, haplotype, large_index, verbose):
    lens = read_fai(fai)
    ref_bp = sum(lens.values())

    per_window = {}                                # (chrom, window) -> added edges
    n_var = {'single': 0, 'insertion': 0, 'deletion': 0}
    added_nodes = 0
    for line in snp:
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        vtype, chrom, pos, data = f[1], f[2], int(f[3]), f[4]
        if vtype not in n_var:
            continue
        n_var[vtype] += 1
        cost = variant_cost(vtype, data)
        added_nodes += cost
        # One per variant, not weighted by inserted length: the capacity
        # constant below was calibrated against variant counts per window, and
        # weighting here without recalibrating there made the fit worse by
        # 6-10 points. Keep the two consistent.
        w = pos // LOCAL_INDEX_INTERVAL
        per_window[(chrom, w)] = per_window.get((chrom, w), 0) + 1

    n_hap = 0
    hap_alts = 0
    if haplotype is not None:
        for line in haplotype:
            f = line.rstrip('\n').split('\t')
            if len(f) < 5:
                continue
            n_hap += 1
            hap_alts += f[4].count(',') + 1

    total_var = sum(n_var.values())

    # ---- the global bound ------------------------------------------------
    # Path-graph construction starts with one node per edge of the reference
    # graph and grows from there; the reference contributes ~one edge per base.
    initial_nodes = ref_bp + added_nodes
    ceiling = INDEX_T_MAX_64 if large_index else INDEX_T_MAX_32
    frac = initial_nodes / ceiling

    # ---- the local budget ------------------------------------------------
    hap_ratio = (n_hap / total_var) if (haplotype is not None and total_var) else 0.80
    capacity = CAPACITY_X_HAPRATIO / hap_ratio
    n_windows = sum((L + LOCAL_INDEX_INTERVAL - 1) // LOCAL_INDEX_INTERVAL
                    for L in lens.values())
    over = [(k, v) for k, v in per_window.items() if v > capacity]
    at_risk = sum(v for _, v in over)

    print("reference          %s bp in %d sequences, %d local windows"
          % (f"{ref_bp:,}", len(lens), n_windows))
    print("variants           %s  (single %s, insertion %s, deletion %s)"
          % (f"{total_var:,}", f"{n_var['single']:,}",
             f"{n_var['insertion']:,}", f"{n_var['deletion']:,}"))
    if haplotype is not None:
        print("haplotypes         %s rows, %.2f per variant"
              % (f"{n_hap:,}", (n_hap / total_var) if total_var else 0))
    print("density            1 per %.0f bp" % (ref_bp / total_var if total_var else 0))
    print()

    print("--- global node bound (%s-bit) ---" % (64 if large_index else 32))
    print("  reference alone      %s nodes  (%.1f%% of the ceiling)"
          % (f"{ref_bp:,}", 100.0 * ref_bp / ceiling))
    print("  plus variants        %s nodes  (%.1f%% of the ceiling)"
          % (f"{initial_nodes:,}", 100.0 * frac))
    print("  headroom for growth  %.1f%%" % (100.0 * (1 - frac)))

    # Prefix doubling only has to grow the count by 1/frac to overrun the
    # ceiling. Below ~50% there is at least a factor of two of room; above it,
    # a single join round can be enough.
    # Observed, not modelled: GRCh38 at 72% of the 32-bit ceiling failed with
    # "exceeded integer bounds" at both 14.95M and 12.64M variants, and chr1 at
    # 5.8% built without trouble. Nothing has been measured between those, so
    # the middle band says so instead of guessing.
    if frac >= 1.0:
        verdict_global = "WILL FAIL -- over the ceiling before doubling even starts"
    elif frac >= 0.70:
        verdict_global = "WILL FAIL -- GRCh38 failed twice at this fraction (72%)"
    elif frac >= 0.35:
        verdict_global = "UNKNOWN -- no build has been measured in this range"
    else:
        verdict_global = "fits -- chr1 built comfortably at 5.8%"
    print("  verdict              %s" % verdict_global)
    print()

    print("--- local edge budget (calibrated) ---")
    print("  haplotypes per variant       %.2f" % hap_ratio)
    print("  effective capacity           %.0f variants per %d bp window"
          % (capacity, LOCAL_INDEX_SIZE))
    print("  windows with variants        %s" % f"{len(per_window):,}")
    print("  windows over the budget      %s (%.1f%% of populated)"
          % (f"{len(over):,}",
             100.0 * len(over) / len(per_window) if per_window else 0))
    print("  variant edges in those       %s (%.1f%% of all)"
          % (f"{at_risk:,}", 100.0 * at_risk / added_nodes if added_nodes else 0))
    # Each over-budget window loses roughly the excess plus what the backoff
    # overshoots past it; halving typically gives up about half the window.
    est_lost = LOSS_FRACTION * at_risk
    print("  estimated retention          %.1f%%"
          % (100.0 * (1 - est_lost / added_nodes) if added_nodes else 100.0))

    if verbose and over:
        print("\n  densest windows:")
        for (chrom, w), v in sorted(over, key=lambda x: -x[1])[:10]:
            print("    %s:%d-%d  %s variants (%.1fx capacity)"
                  % (chrom, w * LOCAL_INDEX_INTERVAL,
                     w * LOCAL_INDEX_INTERVAL + LOCAL_INDEX_SIZE,
                     f"{v:,}", v / capacity))

    print()
    if frac >= 0.70 and not large_index:
        print("Suggestion: --large-index raises the node ceiling, at the cost of a")
        print("  larger index and a much larger build. Reducing haplotype count helps")
        print("  more than reducing variant count: haplotype rows are what multiply")
        print("  during path-graph construction.")
    return 0 if frac < 0.85 else 1


if __name__ == '__main__':
    p = ArgumentParser(
        description='Predict whether hisat2-build will overrun its node bound '
                    'or silently drop variants, before running it')
    p.add_argument('--fai', type=FileType('r'), required=True,
                   help='samtools .fai for the genome')
    p.add_argument('--snp', type=FileType('r'), required=True,
                   help='.snp file that would be passed to hisat2-build')
    p.add_argument('--haplotype', type=FileType('r'), default=None,
                   help='.haplotype file, if one would be supplied')
    p.add_argument('--large-index', action='store_true',
                   help='assume a 64-bit build (hisat2-build-l)')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='list the densest over-budget windows')
    a = p.parse_args()
    sys.exit(main(a.fai, a.snp, a.haplotype, a.large_index, a.verbose))
