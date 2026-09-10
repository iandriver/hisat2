#!/usr/bin/env python3
"""Build a reference containing a spliced gene and a processed retrogene copy.

A processed retrogene is a reverse-transcribed copy of the MATURE mRNA, so a
read spanning the parent's exon-exon junction matches the retrogene
contiguously at the same raw alignment score. HISAT2 resolves this correctly
and reports only the spliced parent alignment, but a counter that unions genes
over every pre-selection candidate sees both and discards the read as
multi-gene. That defect cost Rps27 97% of its UMIs on real mouse data.

The reference is written here rather than reused from example/, because the
example index has no retrogene to trip over.

Layout (0-based):
    [0, pad)            filler
    parent exon1        E1
    intron              I
    parent exon2        E2
    filler
    retrogene           E1 + E2 concatenated, verbatim
    filler
"""
import argparse, os


def read_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    break
                name = line[1:].split()[0]
            elif name is not None:
                chunks.append(line.strip())
    return name, ''.join(chunks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--nreads', type=int, default=60)
    ap.add_argument('--readlen', type=int, default=90)
    ap.add_argument('--with-snp', action='store_true',
                    help='add one SNP in parent exon 1 and emit REF and ALT reads')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    _, src = read_fasta(a.fasta)
    # Real sequence, so composition and repeat content are realistic; skip the
    # leading N run that reference FASTAs usually start with.
    src = src.upper()
    i = 0
    while i < len(src) and src[i] == 'N':
        i += 1
    src = src[i:]
    src = ''.join(c for c in src if c in 'ACGT')
    if len(src) < 20000:
        raise SystemExit("source sequence too short")

    PAD, EX1, INTRON, EX2, GAP = 600, 400, 500, 400, 600
    e1 = src[1000:1000 + EX1]
    intron = src[6000:6000 + INTRON]
    e2 = src[11000:11000 + EX2]
    pad0 = src[16000:16000 + PAD]
    gap = src[17000:17000 + GAP]
    tail = src[18000:18000 + PAD]

    # The retrogene is the spliced mRNA, copied verbatim -- that identity is
    # the whole point of the fixture.
    mrna = e1 + e2
    seq = pad0 + e1 + intron + e2 + gap + mrna + tail
    chrom = 'testchr'

    p_e1_s = PAD
    p_e1_e = p_e1_s + EX1
    p_e2_s = p_e1_e + INTRON
    p_e2_e = p_e2_s + EX2
    retro_s = p_e2_e + GAP
    retro_e = retro_s + len(mrna)

    with open(os.path.join(a.outdir, 'ref.fa'), 'w') as o:
        o.write('>%s\n' % chrom)
        for k in range(0, len(seq), 60):
            o.write(seq[k:k + 60] + '\n')

    with open(os.path.join(a.outdir, 'model.ht2gm'), 'w') as o:
        o.write("#hisat2-gene-model\tv1\n")
        o.write("#ref\t%s\t%d\n" % (chrom, len(seq)))
        o.write("G\t0\tGENE_PARENT\tParent\t%s\t+\t%d\t%d\n" % (chrom, p_e1_s, p_e2_e))
        o.write("E\t0\t%s\t%d\t%d\n" % (chrom, p_e1_s, p_e1_e))
        o.write("E\t0\t%s\t%d\t%d\n" % (chrom, p_e2_s, p_e2_e))
        o.write("J\t0\t%s\t%d\t%d\n" % (chrom, p_e1_e - 1, p_e2_s))
        o.write("G\t1\tGENE_RETRO\tRetro\t%s\t+\t%d\t%d\n" % (chrom, retro_s, retro_e))
        o.write("E\t1\t%s\t%d\t%d\n" % (chrom, retro_s, retro_e))

    # Splice sites and exons in hisat2_extract_splice_sites.py /
    # hisat2_extract_exons.py format (0-based, inclusive), usable both as
    # --known-splicesite-infile and as hisat2-build --ss/--exon. A graph (--snp)
    # index needs the junction built in: without it the aligner seeds on the
    # retrogene's ungapped copy and never explores the spliced parent, so no
    # read would cover the SNP at all.
    with open(os.path.join(a.outdir, 'splicesites.txt'), 'w') as o:
        o.write("%s\t%d\t%d\t+\n" % (chrom, p_e1_e - 1, p_e2_s))
    with open(os.path.join(a.outdir, 'exons.txt'), 'w') as o:
        o.write("%s\t%d\t%d\t+\n" % (chrom, p_e1_s, p_e1_e - 1))
        o.write("%s\t%d\t%d\t+\n" % (chrom, p_e2_s, p_e2_e - 1))
        o.write("%s\t%d\t%d\t+\n" % (chrom, retro_s, retro_e - 1))

    # One SNP in parent exon 1, placed so every read covers it: reads take
    # 20..69 bases from the end of exon 1, so offset EX1-10 is always inside.
    # The retrogene keeps the reference base, so a REF read still aligns
    # perfectly to both loci -- exactly the case the allelic path used to drop,
    # because it demanded a single raw candidate rather than a single
    # top-scoring one.
    SNP_OFF = EX1 - 10
    snp_pos = p_e1_s + SNP_OFF
    ref_base = seq[snp_pos]
    alt_base = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[ref_base]
    if a.with_snp:
        with open(os.path.join(a.outdir, 'ref.snp'), 'w') as o:
            o.write("rsRETRO\tsingle\t%s\t%d\t%s\n" % (chrom, snp_pos, alt_base))

    # Every read straddles the junction, so each one aligns spliced to the
    # parent and contiguously to the retrogene.
    bcs = ['AAACCCAAGAAACACT', 'AAACCCAAGAAACCAT', 'AAACCCAAGAAACCCA']
    expect = {}
    fq = open(os.path.join(a.outdir, 'reads.fq'), 'w')
    n = 0
    for k in range(a.nreads):
        # Vary the split so the junction sits at different offsets.
        left = 20 + (k * 3) % (a.readlen - 40)
        start = EX1 - left
        r = mrna[start:start + a.readlen]
        if len(r) < a.readlen:
            continue
        bc = bcs[k % len(bcs)]
        # The 6-base code is written twice, so any two UMIs differ in at least
        # two positions and 1MM collapsing cannot merge distinct molecules.
        code = ''.join('ACGT'[(k >> (2 * j)) & 3] for j in range(6))
        umi = code + code
        allele = 'ref'
        if a.with_snp and k % 2 == 1:
            i = SNP_OFF - start          # SNP position within the read
            r = r[:i] + alt_base + r[i + 1:]
            allele = 'alt'
        expect[(bc, allele)] = expect.get((bc, allele), 0) + 1
        fq.write('@r%d_%s_%s\n%s\n+\n%s\n' % (k, bc, umi, r, 'I' * len(r)))
        n += 1
    fq.close()

    with open(os.path.join(a.outdir, 'whitelist.txt'), 'w') as o:
        for b in bcs:
            o.write(b + '\n')

    if a.with_snp:
        with open(os.path.join(a.outdir, 'expected_allelic.tsv'), 'w') as o:
            o.write("barcode\trsid\tallele\tumi_count\n")
            for (bc, allele), c in sorted(expect.items()):
                o.write("%s\trsRETRO\t%s\t%d\n" % (bc, allele, c))

    print("len=%d parent=%d-%d/%d-%d retro=%d-%d reads=%d"
          % (len(seq), p_e1_s, p_e1_e, p_e2_s, p_e2_e, retro_s, retro_e, n))


if __name__ == '__main__':
    main()
