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

"""Remove variants from a .snp/.haplotype pair before building a graph index.

A SNP-aware index turns each variant into an alternate path, so a read matching
that path costs no mismatch.  That is the point of the index -- but where a
region is a near-copy of another expressed region, the alternate paths can make
the copy an *exact* match for reads that belong to the original.  The reads do
not move because they fit better; they stop being uniquely placeable at all.

Processed pseudogenes are the case that matters for RNA quantification.  They
are retrotransposed copies of expressed mRNAs, so a read from the parent gene
already nearly matches them, and variants catalogued inside them are partly
paralogous sequence variants -- artefacts of exactly this mapping ambiguity
during variant calling.  Measured on 10x PBMC data against the prebuilt
grch38_snp index, 83% of the reads that a linear index places uniquely at RPS27
lose their unique placement under the graph, with mismatch counts unchanged,
and ribosomal protein genes shed ~103,000 UMIs to their pseudogenes.

This script drops variants that fall inside gene bodies of the chosen biotypes,
so the index keeps its alternate paths everywhere a real allele needs one (the
HLA locus, above all) and stops manufacturing them inside pseudogene copies.

Typical use, between the extract step and the build step:

    hisat2_extract_snps_haplotypes_UCSC.py genome.fa snp151Common.txt genome
    hisat2_filter_snps.py --gtf annotation.gtf \\
        --haplotype genome.haplotype --haplotype-out genome.filtered.haplotype \\
        genome.snp genome.filtered.snp
    hisat2-build genome.fa --snp genome.filtered.snp \\
        --haplotype genome.filtered.haplotype genome_snp

The .snp format is  <id> <type> <chrom> <pos> <allele>, 0-based, where type is
single/deletion/insertion and, for a deletion, <allele> is the deleted length.
The .haplotype format is  <id> <chrom> <left> <right> <comma-separated snp ids>.

hisat2-build already ignores haplotype entries naming unknown variants, so
filtering the .snp file alone is safe; rewriting the haplotype file as well
keeps its spans honest, which is what bounds local-graph construction.
"""

from sys import stderr, exit
from argparse import ArgumentParser, FileType

# Retrotransposed copies of expressed mRNAs. These are the ones that compete
# with their parent gene for reads; unprocessed pseudogenes are duplications
# that retain intron structure and behave far less like the parent's cDNA.
DEFAULT_BIOTYPES = 'processed_pseudogene,transcribed_processed_pseudogene'

# Pseudogenes are not tidily disjoint from real genes: DHFRP2, a transcribed
# processed pseudogene of DHFR, lies across the 3' end of HLA-B and covers its
# last exon -- the very window 3' single-cell chemistry sequences. Masking it
# naively would strip 30% of HLA-B's variants to fix a ribosomal-protein
# problem. Exons of these biotypes are therefore never masked.
DEFAULT_PROTECT = 'protein_coding'


def attr(field, key):
    """Value of a GTF attribute, or None. Handles quoted and bare values."""
    i = field.find(key + ' "')
    if i >= 0:
        i += len(key) + 2
        j = field.find('"', i)
        return field[i:j] if j > 0 else None
    i = field.find(key + ' ')
    if i < 0:
        return None
    i += len(key) + 1
    j = field.find(';', i)
    return field[i:j if j > 0 else len(field)].strip()


def _merge(by_chrom):
    """{chrom: [(s,e)...]} -> {chrom: (starts, ends)} merged, plus total bp."""
    merged = {}
    total = 0
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        out = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                if e > out[-1][1]:
                    out[-1][1] = e
            else:
                out.append([s, e])
        merged[chrom] = ([iv[0] for iv in out], [iv[1] for iv in out])
        total += sum(e - s for s, e in out)
    return merged, total


def scan_gtf(gtf, biotypes, protect, flank, verbose=False):
    """One pass: gene bodies to mask, and exons to protect from masking."""
    want, keep = set(biotypes), set(protect)
    seen = {}
    mask_by, prot_by = {}, {}
    n_mask = n_prot = 0
    for line in gtf:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            bt = attr(f[8], 'gene_biotype') or attr(f[8], 'gene_type')
            if bt is None:
                continue
            seen[bt] = seen.get(bt, 0) + 1
            if bt in want:
                # GTF is 1-based inclusive; store 0-based half-open, flanked.
                mask_by.setdefault(f[0], []).append(
                    (max(0, int(f[3]) - 1 - flank), int(f[4]) + flank))
                n_mask += 1
        elif f[2] == 'exon' and keep:
            bt = attr(f[8], 'gene_biotype') or attr(f[8], 'gene_type')
            if bt in keep:
                prot_by.setdefault(f[0], []).append(
                    (int(f[3]) - 1, int(f[4])))
                n_prot += 1

    if verbose:
        unknown = want - set(seen)
        if unknown:
            print('warning: no genes with biotype %s in the GTF' %
                  ', '.join(sorted(unknown)), file=stderr)
        print('masking %d genes across %d references' % (n_mask, len(mask_by)),
              file=stderr)

    masked, mask_bp = _merge(mask_by)
    protected, prot_bp = _merge(prot_by)
    return masked, mask_bp, n_mask, protected, prot_bp, n_prot


def overlaps(masked, chrom, start, end):
    """Does [start, end) hit a masked interval on chrom? Binary search."""
    iv = masked.get(chrom)
    if iv is None:
        return False
    starts, ends = iv
    # rightmost interval whose start is <= start
    lo, hi = 0, len(starts)
    while lo < hi:
        mid = (lo + hi) // 2
        if starts[mid] <= start:
            lo = mid + 1
        else:
            hi = mid
    i = lo - 1
    if i >= 0 and start < ends[i]:
        return True
    # an interval starting inside [start, end)
    return lo < len(starts) and starts[lo] < end


def variant_span(vtype, pos, allele):
    """0-based half-open span a variant occupies on the reference."""
    if vtype == 'deletion':
        try:
            return pos, pos + int(allele)
        except ValueError:
            return pos, pos + 1
    if vtype == 'insertion':
        # Occupies no reference base; treat as the junction it sits in.
        return pos, pos + 1
    return pos, pos + len(allele) if allele else pos + 1


def filter_snps(snp_in, snp_out, masked, protected, report=None):
    kept_ids = set()
    n_in = n_out = n_rescued = 0
    dropped_by_type = {}
    for line in snp_in:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        n_in += 1
        # hisat2-inspect --snp writes the full FASTA description as the
        # reference name; hisat2-build matches on the first token.
        chrom = f[2].split()[0]
        try:
            pos = int(f[3])
        except ValueError:
            continue
        start, end = variant_span(f[1], pos, f[4])
        if overlaps(masked, chrom, start, end):
            if overlaps(protected, chrom, start, end):
                n_rescued += 1
            else:
                dropped_by_type[f[1]] = dropped_by_type.get(f[1], 0) + 1
                if report is not None:
                    report.write('%s\t%s\t%s\t%d\n' % (f[0], f[1], chrom, pos))
                continue
        kept_ids.add(f[0])
        snp_out.write('\t'.join([f[0], f[1], chrom] + f[3:5]) + '\n')
        n_out += 1
    return n_in, n_out, kept_ids, dropped_by_type, n_rescued


def filter_haplotypes(ht_in, ht_out, kept_ids, snp_pos):
    """Drop removed variants from each haplotype and retighten its span.

    A haplotype naming no surviving variant is dropped entirely, which is what
    hisat2-build would do anyway; recomputing left/right matters because those
    bounds drive local graph construction.
    """
    n_in = n_out = n_trimmed = 0
    for line in ht_in:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        n_in += 1
        ids = [i for i in f[4].split(',') if i in kept_ids]
        if not ids:
            continue
        if len(ids) != len(f[4].split(',')):
            n_trimmed += 1
        chrom = f[1].split()[0]
        spans = [snp_pos[i] for i in ids if i in snp_pos]
        if spans:
            left = min(s for s, _ in spans)
            right = max(e for _, e in spans) - 1
        else:
            left, right = int(f[2]), int(f[3])
        ht_out.write('%s\t%s\t%d\t%d\t%s\n' %
                     ('ht%d' % n_out, chrom, left, max(left, right), ','.join(ids)))
        n_out += 1
    return n_in, n_out, n_trimmed


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Remove variants inside pseudogenes (or other chosen gene '
                    'biotypes) from a .snp/.haplotype pair before building a '
                    'HISAT2 graph index')
    parser.add_argument('snp_in', type=FileType('r'),
        help='input .snp file (use "-" for stdin)')
    parser.add_argument('snp_out', type=FileType('w'),
        help='output .snp file')
    parser.add_argument('--gtf', dest='gtf', type=FileType('r'), required=True,
        help='GTF annotation supplying the gene bodies to mask')
    parser.add_argument('--biotypes', dest='biotypes', default=DEFAULT_BIOTYPES,
        help='comma-separated gene biotypes to mask (default: %s)' %
             DEFAULT_BIOTYPES)
    parser.add_argument('--protect', dest='protect', default=DEFAULT_PROTECT,
        help='never mask exons of genes with these biotypes, even where a '
             'pseudogene overlaps them; empty string to disable '
             '(default: %s)' % DEFAULT_PROTECT)
    parser.add_argument('--flank', dest='flank', type=int, default=0,
        help='also mask this many bp either side of each gene (default: 0)')
    parser.add_argument('--haplotype', dest='ht_in', type=FileType('r'),
        help='input .haplotype file to filter alongside the .snp file')
    parser.add_argument('--haplotype-out', dest='ht_out', type=FileType('w'),
        help='output .haplotype file; required with --haplotype')
    parser.add_argument('--report', dest='report', type=FileType('w'),
        help='write the dropped variants here, one per line')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print statistics to stderr')

    args = parser.parse_args()
    if args.ht_in and not args.ht_out:
        parser.error('--haplotype requires --haplotype-out')

    biotypes = [b.strip() for b in args.biotypes.split(',') if b.strip()]
    if not biotypes:
        parser.error('--biotypes is empty')

    protect = [b.strip() for b in args.protect.split(',') if b.strip()]

    masked, masked_bp, n_genes, protected, prot_bp, n_ex = scan_gtf(
        args.gtf, biotypes, protect, args.flank, args.verbose)
    if not masked:
        print('error: no genes matched biotypes %s' % ', '.join(biotypes),
              file=stderr)
        exit(1)

    # Retightening haplotype spans needs every variant's coordinates, so with
    # --haplotype the .snp file is read once into memory first.
    snp_pos = {}
    if args.ht_in:
        lines = args.snp_in.readlines()
        for line in lines:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 5:
                try:
                    snp_pos[f[0]] = variant_span(f[1], int(f[3]), f[4])
                except ValueError:
                    pass
        args.snp_in = lines

    n_in, n_out, kept, by_type, n_rescued = filter_snps(
        args.snp_in, args.snp_out, masked, protected, args.report)
    args.snp_out.close()

    if args.verbose:
        print('masked %d genes, %.1f Mb' % (n_genes, masked_bp / 1e6),
              file=stderr)
        if protect:
            print('protected %d %s exons, %.1f Mb; %d variants kept that the '
                  'mask would otherwise have dropped' %
                  (n_ex, '/'.join(protect), prot_bp / 1e6, n_rescued),
                  file=stderr)
        print('variants: %d in, %d out, %d dropped (%.2f%%)' %
              (n_in, n_out, n_in - n_out,
               100.0 * (n_in - n_out) / n_in if n_in else 0.0), file=stderr)
        for t in sorted(by_type):
            print('  dropped %-10s %d' % (t, by_type[t]), file=stderr)

    if args.ht_in:
        h_in, h_out, h_trim = filter_haplotypes(
            args.ht_in, args.ht_out, kept, snp_pos)
        args.ht_out.close()
        if args.verbose:
            print('haplotypes: %d in, %d out, %d retained but trimmed' %
                  (h_in, h_out, h_trim), file=stderr)
