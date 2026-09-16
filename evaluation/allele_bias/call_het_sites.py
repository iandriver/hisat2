#!/usr/bin/env python3
"""Call candidate heterozygous sites from a pileup, when genotypes are unknown.

  call_het_sites.py --pileup linear.pileup --snp-sites sites.tsv --out het.vcf

`--snp-sites` is a four-column table (contig, 1-based position, id, alt allele),
built from the graph index's own variant list:

  hisat2-inspect --snp INDEX \\
    | awk -F'\\t' '$2=="single"{split($3,a," "); print a[1]"\\t"($4+1)"\\t"$1"\\t"$5}'

Call the sites from the LINEAR arm. Calling them from the graph arm would select
sites where the graph had already found alternate reads, which is the result the
experiment is supposed to test. Calling them from linear is conservative in the
other direction: sites where the linear aligner recovered no alternate read at
all are never considered, so the graph's advantage is understated.

Without genotypes some called sites will really be homozygous with mapping noise.
That error is shared by every arm, since all arms are scored at the same sites.
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import allele_bias as ab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pileup', required=True)
    ap.add_argument('--snp-sites', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-depth', type=int, default=ab.MIN_DEPTH)
    ap.add_argument('--lo', type=float, default=0.10, help='min alt fraction')
    ap.add_argument('--hi', type=float, default=0.90, help='max alt fraction')
    ap.add_argument('--chr-prefix', action='store_true',
                    help="write contigs as chr1 rather than 1")
    ap.add_argument('--contig-order', metavar='FILE',
                    help='order output by this contig list (STAR wants genome order)')
    ap.add_argument('--sample', default='SAMPLE')
    a = ap.parse_args()

    alt_of = {}
    for line in open(a.snp_sites):
        f = line.rstrip('\n').split('\t')
        alt_of[(f[0].replace('chr', ''), int(f[1]))] = f[3].upper()

    rows = []
    for line in open(a.pileup):
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        key = (f[0].replace('chr', ''), int(f[1]))
        alt = alt_of.get(key)
        ref_base = f[2].upper()
        if alt is None or ref_base not in 'ACGT':
            continue
        c = collections.Counter(ch.upper() for ch in ab.parse_pileup_bases(f[4]))
        ref = c.get('.', 0) + c.get(',', 0)
        n_alt = c.get(alt, 0)
        if ref + n_alt < a.min_depth:
            continue
        frac = n_alt / (ref + n_alt)
        if not (a.lo <= frac <= a.hi):
            continue
        rows.append((key[0], key[1], ref_base, alt))

    if a.contig_order:
        order = {l.strip(): i for i, l in enumerate(open(a.contig_order))}
        rows.sort(key=lambda r: (order.get(('chr' if a.chr_prefix else '') + r[0], 1 << 30), r[1]))
    else:
        rows.sort()

    with open(a.out, 'w') as o:
        o.write("##fileformat=VCFv4.2\n")
        o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        o.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a.sample}\n")
        for c, p, r, alt in rows:
            name = f"chr{c}" if a.chr_prefix else c
            o.write(f"{name}\t{p}\t.\t{r}\t{alt}\t.\tPASS\t.\tGT\t0|1\n")
    print(f"het sites called: {len(rows):,} -> {a.out}")


if __name__ == '__main__':
    main()
