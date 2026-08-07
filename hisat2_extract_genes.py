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

"""Extract a gene model (.ht2gm) from a GTF file for single-cell quantification.

Unlike hisat2_extract_exons.py, which merges exons ACROSS genes and discards
gene identity, this retains gene_id/gene_name and merges exons only WITHIN a
gene.  The result is a gene model suitable for assigning alignments to genes.

Output format (tab-separated, 0-based half-open coordinates):

    #hisat2-gene-model  v1
    #ref   <name>  <length|->
    G  <gene_idx>  <gene_id>  <gene_name>  <chrom>  <strand>  <start>  <end>
    E  <gene_idx>  <chrom>  <start>  <end>
    J  <gene_idx>  <chrom>  <donor>  <acceptor>

G gives the gene body (min exon start to max exon end) used for GeneFull
counting.  E gives the per-gene exon union used for Gene counting.  J gives
annotated junctions, used to classify spliced vs unspliced reads.

The #ref lines let the loader verify the model was built against the same
assembly as the index.  Lengths are only emitted when --fai is supplied.
"""

from sys import stderr, stdout, exit
from collections import defaultdict as dd
from argparse import ArgumentParser, FileType


def parse_gtf(gtf_file):
    """Parse exon records, keyed by transcript, retaining gene identity.

    The line-level parsing here deliberately mirrors
    hisat2_extract_splice_sites.py so that the junctions this emits agree with
    the ones fed to --known-splicesite-infile.
    """
    trans = {}                 # transcript_id -> [chrom, strand, gene_key, exons]
    gene_names = {}            # gene_key -> gene_name
    skipped = dd(int)

    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            skipped['malformed'] += 1
            continue
        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or 'transcript_id' not in values_dict:
            skipped['no_gene_or_transcript_id'] += 1
            continue

        gene_id = values_dict['gene_id']
        transcript_id = values_dict['transcript_id']
        # A gene_id can legitimately appear on more than one sequence (GENCODE
        # PAR_Y).  Key on the triple so those stay distinct rather than being
        # silently fused into one impossible interval spanning two chromosomes.
        gene_key = (gene_id, chrom, strand)
        gene_names.setdefault(gene_key, values_dict.get('gene_name', gene_id))

        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, gene_key, [[left, right]]]
        else:
            trans[transcript_id][3].append([left, right])

    return trans, gene_names, skipped


def merge_intervals(ivs):
    """Union of [start, end) intervals; merges overlapping and abutting ones."""
    ivs.sort()
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def build_gene_model(trans, gene_names):
    """Collapse transcripts into per-gene exon unions, bodies, and junctions."""
    gene_exons = dd(list)      # gene_key -> [[start, end), ...]
    gene_juncs = dd(set)       # gene_key -> {(donor, acceptor)}

    for chrom, strand, gene_key, exons in trans.values():
        exons.sort()
        # GTF is 1-based inclusive; convert to 0-based half-open.
        for left, right in exons:
            gene_exons[gene_key].append((left - 1, right))
        # Junction convention matches hisat2_extract_splice_sites.py: donor is
        # the 0-based last base of the upstream exon, acceptor the 0-based
        # first base of the downstream exon.
        for i in range(1, len(exons)):
            gene_juncs[gene_key].add((exons[i - 1][1] - 1, exons[i][0] - 1))

    genes = []
    for gene_key in gene_exons:
        gene_id, chrom, strand = gene_key
        merged = merge_intervals(gene_exons[gene_key])
        genes.append({
            'gene_id': gene_id,
            'gene_name': gene_names[gene_key],
            'chrom': chrom,
            'strand': strand,
            'start': merged[0][0],
            'end': max(e for _, e in merged),
            'exons': merged,
            'juncs': sorted(gene_juncs[gene_key]),
        })

    # Gene order is GTF order -- gene_exons is a dict, so iterating it above
    # yields genes in the order they first appear in the file. That is what
    # STARsolo and rustar use, so features.tsv lines up row-for-row with
    # theirs and matrices can be compared positionally as well as by name.
    # It is no less deterministic than sorting: the same GTF gives the same
    # order every time. Sorting by coordinate here instead, which this used
    # to do, silently produced a different row order from every other tool.
    return genes


def read_fai(fai_file):
    lengths = {}
    for line in fai_file:
        fields = line.rstrip('\n').split('\t')
        if len(fields) >= 2:
            lengths[fields[0]] = fields[1]
    return lengths


def write_model(genes, lengths, out):
    out.write('#hisat2-gene-model\tv1\n')
    for chrom in sorted({g['chrom'] for g in genes}):
        out.write('#ref\t{}\t{}\n'.format(chrom, lengths.get(chrom, '-')))
    for idx, g in enumerate(genes):
        out.write('G\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            idx, g['gene_id'], g['gene_name'], g['chrom'], g['strand'],
            g['start'], g['end']))
        for s, e in g['exons']:
            out.write('E\t{}\t{}\t{}\t{}\n'.format(idx, g['chrom'], s, e))
        for d, a in g['juncs']:
            out.write('J\t{}\t{}\t{}\t{}\n'.format(idx, g['chrom'], d, a))


def report(genes, skipped, verbose):
    n_exons = sum(len(g['exons']) for g in genes)
    n_juncs = sum(len(g['juncs']) for g in genes)
    print('genes: {}, exon intervals: {}, junctions: {}'.format(
        len(genes), n_exons, n_juncs), file=stderr)

    # A gene_id yielding more than one record means it appears on more than one
    # sequence or strand (GENCODE PAR_Y).  Surfacing it matters: downstream
    # features.tsv will contain the id twice.
    by_id = dd(int)
    for g in genes:
        by_id[g['gene_id']] += 1
    dupes = {k: v for k, v in by_id.items() if v > 1}
    if dupes:
        print('warning: {} gene_id(s) span multiple sequences/strands and were '
              'kept as separate records, e.g. {}'.format(
                  len(dupes), ', '.join(sorted(dupes)[:3])), file=stderr)
    for reason, count in sorted(skipped.items()):
        print('warning: skipped {} line(s): {}'.format(count, reason),
              file=stderr)

    if verbose:
        lens = sorted(g['end'] - g['start'] for g in genes)
        exon_bp = sum(e - s for g in genes for s, e in g['exons'])
        print('gene body length: median {}, max {}'.format(
            lens[len(lens) // 2], lens[-1]), file=stderr)
        print('total exonic bp: {}'.format(exon_bp), file=stderr)
        print('single-exon genes: {}'.format(
            sum(1 for g in genes if len(g['exons']) == 1)), file=stderr)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract a gene model (.ht2gm) from a GTF file')
    parser.add_argument('gtf_file',
        nargs='?',
        type=FileType('r'),
        help='input GTF file (use "-" for stdin)')
    parser.add_argument('-o', '--output',
        dest='output',
        type=FileType('w'),
        default=stdout,
        help='output .ht2gm file (default: stdout)')
    parser.add_argument('--fai',
        dest='fai',
        type=FileType('r'),
        help='samtools faidx .fai for the genome, so reference lengths are '
             'recorded and the loader can verify the model matches the index')
    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)

    trans, gene_names, skipped = parse_gtf(args.gtf_file)
    if not trans:
        print('error: no usable exon records found in GTF', file=stderr)
        exit(1)
    genes = build_gene_model(trans, gene_names)
    lengths = read_fai(args.fai) if args.fai else {}
    write_model(genes, lengths, args.output)
    report(genes, skipped, args.verbose)
