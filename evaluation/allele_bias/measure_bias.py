#!/usr/bin/env python3
"""Compare reference bias across aligner configurations at the same het sites.

  measure_bias.py --vcf het.vcf --bam linear=linear.bam --bam graph=graph.bam \
                  [--counts "linear + WASP"=wasp_counts.tsv]

Arms given as --bam are counted from the alignments; arms given as --counts read
a per-site table written by wasp_filter.py. Every arm is held to the site list of
the first arm, so a difference in depth between arms cannot be mistaken for a
difference in denominator.
"""
import argparse
import collections
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import allele_bias as ab


def read_counts_tsv(path):
    counts = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        c, p, r, a = line.split()
        counts[(c.replace('chr', ''), int(p))] = [int(r), int(a)]
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--bam', action='append', default=[], metavar='NAME=PATH')
    ap.add_argument('--counts', action='append', default=[], metavar='NAME=PATH')
    ap.add_argument('--min-depth', type=int, default=ab.MIN_DEPTH)
    ap.add_argument('--no-groups', action='store_true', help='skip the MHC split')
    a = ap.parse_args()

    sites = ab.load_sites(a.vcf)
    by_chrom = ab.index_by_chrom(sites)
    print(f"het sites: {len(sites):,}")

    tmp = tempfile.mkdtemp(prefix='allele_bias.')
    beds = {}
    for prefix in (False, True):
        beds[prefix] = os.path.join(tmp, f"sites{'_chr' if prefix else ''}.bed")
        ab.write_bed(sites, beds[prefix], prefix)

    arms = collections.OrderedDict()
    for spec in a.bam:
        name, path = spec.split('=', 1)
        # Detected, not declared: a region file that disagrees with the BAM's
        # contig naming selects nothing and looks like an empty arm.
        counts, _ = ab.scan_bam(path, sites, by_chrom, beds[ab.bam_uses_chr_prefix(path)])
        arms[name] = counts
    for spec in a.counts:
        name, path = spec.split('=', 1)
        arms[name] = read_counts_tsv(path)
    if not arms:
        sys.exit("no arms given")

    groups = [("all sites", list(sites))]
    if not a.no_groups:
        groups += [("MHC", [k for k in sites if ab.in_mhc(k)]),
                   ("non-MHC", [k for k in sites if not ab.in_mhc(k)])]

    first = next(iter(arms))
    for gname, gkeys in groups:
        gset = set(gkeys)
        base_keys, _, base_depth = None, None, None
        print(f"\n[{gname}]  (unbiased = 0.5000)")
        print(ab.HEADER)
        for name, counts in arms.items():
            restrict = base_keys if base_keys is not None else [k for k in gkeys]
            keys, mean, depth = ab.summarize(name, counts, restrict, a.min_depth, base_depth)
            if base_keys is None and name == first:
                base_keys, base_depth = keys, depth


if __name__ == '__main__':
    main()
