#!/usr/bin/env python3
"""WASP-style bias correction, run over whichever aligner you point it at.

WASP (van de Geijn et al. 2015, doi:10.1038/nmeth.3582) removes reference bias
by swapping the allele in every read that covers a heterozygous site, realigning,
and discarding reads that fail to return to the same locus. It corrects by
throwing reads away, so its effect has to be read alongside the depth it costs.

Running it over the same aligner being evaluated is what isolates the correction
method from the aligner:

  wasp_filter.py --vcf het.vcf --bam linear.bam --aligner hisat2 \
                 --aligner-path ./hisat2 --index /path/grch38/genome \
                 --out wasp_linear

writes <out>.counts.tsv for measure_bias.py --counts.
"""
import argparse
import collections
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import allele_bias as ab


def remap_cmd(aligner, path, index, fq, out_prefix, threads):
    if aligner == 'hisat2':
        return [path, '-x', index, '-U', fq, '-p', str(threads),
                '--no-unal', '-S', f'{out_prefix}.remap.sam'], f'{out_prefix}.remap.sam'
    if aligner == 'rustar':
        return [path, '--runMode', 'alignReads', '--genomeDir', index,
                '--readFilesIn', fq, '--runThreadN', str(threads),
                '--outSAMtype', 'SAM', '--outFileNamePrefix', f'{out_prefix}.'], \
               f'{out_prefix}.Aligned.out.sam'
    raise SystemExit(f"unknown aligner {aligner}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--aligner', choices=['hisat2', 'rustar'], required=True)
    ap.add_argument('--aligner-path', required=True)
    ap.add_argument('--index', required=True)
    ap.add_argument('--out', required=True, help='output prefix')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--reuse-remap', metavar='SAM',
                    help='skip realignment and read this SAM instead')
    a = ap.parse_args()

    sites = ab.load_sites(a.vcf)
    by_chrom = ab.index_by_chrom(sites)
    tmp = tempfile.mkdtemp(prefix='wasp_filter.')
    bed = os.path.join(tmp, 'sites.bed')
    ab.write_bed(sites, bed, ab.bam_uses_chr_prefix(a.bam))
    counts, recs = ab.scan_bam(a.bam, sites, by_chrom, bed)
    print(f"reads covering a het site: {len(recs):,}")

    fq = f'{a.out}.swapped.fq'
    with open(fq, 'w') as o:
        for i, r in enumerate(recs):
            ref, alt = sites[(r['chrom'], r['site'])]
            swapped = alt if r['allele'] == 0 else ref
            seq = r['seq'][:r['offset']] + swapped + r['seq'][r['offset'] + 1:]
            o.write(f"@{i}\n{seq}\n+\n{r['qual']}\n")   # name is the record index

    sam = a.reuse_remap
    if not sam:
        cmd, sam = remap_cmd(a.aligner, a.aligner_path, a.index, fq, a.out, a.threads)
        subprocess.run(cmd, check=True,
                       stdout=open(f'{a.out}.remap.out', 'w'),
                       stderr=open(f'{a.out}.remap.log', 'w'))

    best = {}
    for line in open(sam):
        if line.startswith('@'):
            continue
        f = line.split('\t')
        if int(f[1]) & 0x900:
            continue
        best[int(f[0])] = (f[2].replace('chr', ''), int(f[3]), int(f[4]))

    kept = collections.defaultdict(lambda: [0, 0])
    npass = 0
    for i, r in enumerate(recs):
        b = best.get(i)
        # WASP's rule: the allele-swapped read must come back to the same locus.
        if b and b[0] == r['chrom'] and b[1] == r['pos'] and b[2] >= ab.MIN_MAPQ:
            kept[(r['chrom'], r['site'])][r['allele']] += 1
            npass += 1
    print(f"passed remapping: {npass:,}/{len(recs):,} ({npass/max(len(recs),1):.1%})")

    with open(f'{a.out}.counts.tsv', 'w') as o:
        o.write("#chrom\tpos\tref_count\talt_count\n")
        for (c, p), (r, alt) in sorted(kept.items()):
            o.write(f"{c}\t{p}\t{r}\t{alt}\n")
    print(f"wrote {a.out}.counts.tsv")
    print(ab.HEADER)
    ab.summarize('before WASP', counts)
    ab.summarize('after WASP', kept, list(counts))


if __name__ == '__main__':
    main()
