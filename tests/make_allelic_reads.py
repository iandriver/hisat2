#!/usr/bin/env python3
"""Generate barcoded reads with known REF/ALT alleles for validating
--solo-allelic against a hand-computed expectation.

Reads are exact substrings of the reference (optionally with the alternate
base substituted at one SNP), so the true allele of every read is known by
construction rather than inferred.
"""

import argparse
import os
import sys


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].strip().split()[0]
        else:
            seq.append(line.strip())
    return name, ''.join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--snp', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--nsnps', type=int, default=20)
    ap.add_argument('--readlen', type=int, default=100)
    args = ap.parse_args()

    chrom, seq = read_fasta(args.fasta)
    snps = []
    for line in open(args.snp):
        f = line.rstrip('\n').split('\t')
        if len(f) < 5 or f[1] != 'single':
            continue
        snps.append((f[0], int(f[3]), f[4]))   # rsID, 0-based pos, alt base

    half = args.readlen // 2
    # Keep SNPs far enough apart that no read spans two of them; otherwise the
    # expected counts would depend on multi-SNP reads and stop being obvious.
    chosen, lastpos = [], -10**9
    for rs, pos, alt in snps:
        if pos - lastpos < 3 * args.readlen:
            continue
        if pos < half or pos + half >= len(seq):
            continue
        ref = seq[pos].upper()
        if ref == alt.upper() or ref not in 'ACGT':
            continue
        chosen.append((rs, pos, ref, alt.upper()))
        lastpos = pos
        if len(chosen) >= args.nsnps:
            break

    if len(chosen) < args.nsnps:
        print("only found %d usable SNPs" % len(chosen), file=sys.stderr)

    barcodes = ['ACGTACGTACGTACGT', 'TTTTGGGGCCCCAAAA',
                'GGGGCCCCTTTTAAAA', 'CCCCAAAAGGGGTTTT']

    os.makedirs(args.outdir, exist_ok=True)
    fq = open(os.path.join(args.outdir, 'reads.fq'), 'w')
    expected = open(os.path.join(args.outdir, 'expected.tsv'), 'w')
    expected.write("barcode\trsid\tallele\tumi_count\n")

    umi_n = 0
    rid = 0
    for si, (rs, pos, ref, alt) in enumerate(chosen):
        start = pos - half
        base = seq[start:start + args.readlen].upper()
        if len(base) != args.readlen or 'N' in base:
            continue
        off = pos - start
        altseq = base[:off] + alt + base[off + 1:]
        for bi, bc in enumerate(barcodes):
            # A deterministic but varied number of molecules per (cell, variant).
            nref = 1 + (si + bi) % 3
            nalt = 1 + (si + 2 * bi) % 2
            for which, n, s in (('ref', nref, base), ('alt', nalt, altseq)):
                for k in range(n):
                    umi_n += 1
                    # Base-4 encode the counter so every UMI is distinct; a
                    # decimal-digit mapping would collide and silently deflate
                    # the expected counts.
                    v, umi = umi_n, ''
                    for _ in range(12):
                        umi = 'ACGT'[v & 3] + umi
                        v >>= 2
                    rid += 1
                    fq.write("@r%d.%s_%s_%s\n%s\n+\n%s\n" % (rid, which, bc, umi, s, 'I' * len(s)))
                expected.write("%s\t%s\t%s\t%d\n" % (bc, rs, which, n))
    fq.close()
    expected.close()

    # A single gene spanning the region, so every read is assigned and the
    # allelic path is exercised independently of gene-model subtleties.
    with open(os.path.join(args.outdir, 'model.ht2gm'), 'w') as o:
        o.write("#hisat2-gene-model\tv1\n")
        o.write("#ref\t%s\t%d\n" % (chrom, len(seq)))
        o.write("G\t0\tGENE_ALL\tAll\t%s\t+\t0\t%d\n" % (chrom, len(seq)))
        o.write("E\t0\t%s\t0\t%d\n" % (chrom, len(seq)))

    with open(os.path.join(args.outdir, 'whitelist.txt'), 'w') as o:
        for bc in barcodes:
            o.write(bc + "\n")

    print("chrom=%s len=%d snps_used=%d reads=%d" % (chrom, len(seq), len(chosen), rid))


if __name__ == '__main__':
    main()
