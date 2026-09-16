"""Shared pieces for measuring reference bias at heterozygous sites.

Reference bias is the tendency of an aligner to recover the reference allele
more readily than the alternate one, because an alternate-allele read mismatches
the reference and so aligns worse, or not at all. At a heterozygous site an
unbiased aligner recovers both alleles equally, so the reference fraction sits
at 0.5 and anything above that is bias.

Every arm is counted through `scan_bam` so that only the intended thing differs
between them.
"""
import bisect
import os
import collections
import subprocess

MIN_MAPQ = 10      # uniquely-mapped reads, which is what allele-specific pipelines use
MIN_BASEQ = 13
MIN_DEPTH = 20


def load_sites(vcf):
    """{(chrom, pos): (ref, alt)} from a VCF. Contig names are normalised without
    the 'chr' prefix, because the indexes involved disagree about it."""
    sites = {}
    for line in open(vcf):
        if line.startswith('#'):
            continue
        f = line.split('\t')
        sites[(f[0].replace('chr', ''), int(f[1]))] = (f[3].upper(), f[4].upper())
    return sites


def index_by_chrom(sites):
    by = collections.defaultdict(list)
    for c, p in sites:
        by[c].append(p)
    for c in by:
        by[c].sort()
    return by


def bam_uses_chr_prefix(bam):
    """Whether a BAM's contigs are named chr1 or 1. Passing a region file whose
    names disagree with the header selects nothing, silently, and looks exactly
    like an arm with no coverage."""
    if not os.path.exists(bam):
        raise SystemExit(f"no such BAM: {bam}")
    r = subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        raise SystemExit(f"samtools idxstats failed on {bam}: {r.stderr.strip()}")
    return r.stdout.split('\t', 1)[0].startswith('chr')


def write_bed(sites, path, chr_prefix):
    with open(path, 'w') as o:
        for (c, p) in sorted(sites):
            o.write(f"{'chr' if chr_prefix else ''}{c}\t{p-1}\t{p}\n")


def ref_offset(cigar, start, target):
    """Index within the read of reference position `target`, or None when the
    read does not have a base there (a deletion or an intron)."""
    ref = start
    read = 0
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ''
        if ch in 'M=X':
            if ref <= target < ref + n:
                return read + (target - ref)
            ref += n
            read += n
        elif ch in 'IS':
            read += n
        elif ch in 'DN':
            if ref <= target < ref + n:
                return None
            ref += n
    return None


def parse_pileup_bases(s):
    """Strip mpileup's read-start, read-end and indel run-length markers."""
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == '^':
            i += 2
            continue
        if c == '$':
            i += 1
            continue
        if c in '+-':
            j = i + 1
            n = ''
            while j < len(s) and s[j].isdigit():
                n += s[j]
                j += 1
            i = j + int(n)
            continue
        out.append(c)
        i += 1
    return out


def scan_bam(bam, sites, by_chrom, bed, min_mapq=MIN_MAPQ, min_baseq=MIN_BASEQ):
    """Per-site [ref, alt] counts, plus one record per counted read so a caller
    can rewrite those reads (WASP needs them). Secondary, supplementary and
    duplicate alignments are excluded."""
    counts = collections.defaultdict(lambda: [0, 0])
    recs = []
    if not os.path.exists(bam):
        raise SystemExit(f"no such BAM: {bam}")
    p = subprocess.Popen(
        ["samtools", "view", "-q", str(min_mapq), "-F", "0x904", "-L", bed, bam],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split('\t')
        chrom, pos, cigar, seq, qual = f[2].replace('chr', ''), int(f[3]), f[5], f[9], f[10]
        positions = by_chrom.get(chrom)
        if not positions:
            continue
        end = pos + len(seq) + 200          # generous: an intron can push the far end out
        for i in range(bisect.bisect_left(positions, pos), len(positions)):
            site = positions[i]
            if site > end:
                break
            off = ref_offset(cigar, pos, site)
            if off is None or off >= len(seq):
                continue
            if ord(qual[off]) - 33 < min_baseq:
                continue
            ref, alt = sites[(chrom, site)]
            base = seq[off].upper()
            which = 0 if base == ref else (1 if base == alt else None)
            if which is None:                # a third base: sequencing error, not informative
                continue
            counts[(chrom, site)][which] += 1
            recs.append(dict(name=f[0], chrom=chrom, pos=pos, site=site,
                             allele=which, offset=off, seq=seq, qual=qual))
    p.wait()
    if p.returncode != 0:
        raise SystemExit(f"samtools view failed on {bam}: {p.stderr.read().strip()}")
    if not recs:
        # Almost always a contig-naming mismatch between the BED and the BAM,
        # which otherwise reports as an arm with no coverage.
        raise SystemExit(f"no reads counted in {bam}; check contig naming and the site list")
    return counts, recs


def in_mhc(key):
    return key[0] == '6' and 28_500_000 <= key[1] <= 33_400_000


HEADER = (f"{'arm':<26}{'sites':>7}{'meanREF':>10}{'bias':>9}"
          f"{'REF':>9}{'ALT':>9}{'depth':>9}{'vs 1st':>9}")


def summarize(tag, counts, restrict_to=None, min_depth=MIN_DEPTH, baseline_depth=None):
    """Mean reference fraction over sites with enough depth. `restrict_to` holds
    every arm to the same site list, so depth differences between arms cannot be
    confused with a different denominator."""
    keys = [k for k in (restrict_to if restrict_to is not None else counts)
            if k in counts and sum(counts[k]) >= min_depth]
    if not keys:
        print(f"  {tag:<24}{0:>7}{'n/a':>10}")
        return [], float('nan'), 0
    fracs = [counts[k][0] / sum(counts[k]) for k in keys]
    ref = sum(counts[k][0] for k in keys)
    alt = sum(counts[k][1] for k in keys)
    mean = sum(fracs) / len(fracs)
    depth = ref + alt
    delta = f"{depth/baseline_depth - 1:>+9.1%}" if baseline_depth else ""
    print(f"  {tag:<24}{len(keys):>7}{mean:>10.4f}{mean-0.5:>+9.4f}"
          f"{ref:>9,}{alt:>9,}{depth:>9,}{delta}")
    return keys, mean, depth
