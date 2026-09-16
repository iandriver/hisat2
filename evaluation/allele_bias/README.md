# Reference bias at heterozygous sites

Measures how much each aligner configuration over-recovers the reference allele,
and what it costs to correct. At a heterozygous site an unbiased aligner recovers
both alleles equally, so the reference fraction sits at 0.5; anything above is
bias. These are the scripts behind the table in `README_SOLO.md`.

## Result

10M reads from a 10x human PBMC 5k 3' v3 library against GRCh38, aligned to
JHU's `grch38` (linear) and `grch38_snp` (graph, 15,311,803 SNPs), which differ
only in whether variants are in the index, plus rustar against a
CellRanger-derived index. 1,504 heterozygous sites, counted through one code
path (MAPQ >= 10, base quality >= 13, primary alignments only).

| arm | sites | mean REF | bias | REF | ALT | depth |
|---|---|---|---|---|---|---|
| HISAT2 linear | 1,480 | 0.5426 | +0.0426 | 67,634 | 61,393 | 129,027 |
| HISAT2 linear + WASP | 1,363 | 0.5172 | +0.0172 | 60,135 | 60,036 | 120,171 |
| rustar linear | 1,472 | 0.5248 | +0.0248 | 71,304 | 70,519 | 141,823 |
| rustar + WASP | 1,446 | 0.5108 | +0.0108 | 66,624 | 70,131 | 136,755 |
| HISAT2 graph | 1,471 | **0.5042** | **+0.0042** | 67,835 | 73,204 | 141,039 |

The graph index is the least biased, and the only arm that reduces bias without
losing reads. Both WASP arms correct by discarding reads, costing 6.9% of depth
on HISAT2 and 3.6% on rustar. The graph instead recovers alternate-allele reads:
against HISAT2 linear, REF is flat (67,634 to 67,835) while ALT rises 19%
(61,393 to 73,204).

Mind the baseline. HISAT2's own linear mode is more reference-biased than
rustar's (+0.0426 against +0.0248) and finds 13% fewer alternate reads, so the
improvement reads as ~90% against HISAT2 linear but ~83% against rustar, and the
depth gain disappears against rustar. What holds against every baseline is the
comparison with WASP: less residual bias at more depth.

## Reproducing

Site list from the graph index's own variants:

    hisat2-inspect --snp $GRAPH_INDEX \
      | awk -F'\t' '$2=="single"{split($3,a," "); print a[1]"\t"($4+1)"\t"$1"\t"$5}' \
      > sites.tsv
    cut -f1,2 sites.tsv > sites.pos

Align each arm the same way, then pile up at those sites:

    hisat2 -x $INDEX -U reads.fq -p 10 --no-unal -S - | samtools sort -o ARM.bam -
    samtools index ARM.bam
    samtools mpileup -l sites.pos -f genome.fa -d 100000 -Q 13 -q 10 --no-BAQ ARM.bam > ARM.pileup

`--no-BAQ` matters: BAQ re-scores alignments around mismatches, which is exactly
the signal being compared.

Call heterozygous sites from the **linear** arm, then score every arm at those
same sites:

    call_het_sites.py --pileup linear.pileup --snp-sites sites.tsv --out het.vcf
    wasp_filter.py --vcf het.vcf --bam linear.bam --aligner hisat2 \
        --aligner-path ./hisat2 --index $LINEAR_INDEX --out wasp_linear
    measure_bias.py --vcf het.vcf \
        --bam linear=linear.bam --bam graph=graph.bam \
        --counts "linear + WASP"=wasp_linear.counts.tsv

`measure_bias.py` detects each BAM's contig naming, because a region file that
says `6` against a BAM that says `chr6` selects nothing and looks identical to
an arm with no coverage.

## What is weak about this

- **No genotypes.** Heterozygous sites are inferred from the data, so some are
  really homozygous with mapping noise. Every arm is scored at the same sites, so
  the error is shared, but a genotyped donor would be better.
- **Ascertainment.** Sites come from the linear arm, so sites where the linear
  aligner found no alternate read at all never enter. That understates the graph.
- **WASP is reimplemented here** (swap the allele, realign, require the same
  locus) rather than run as the reference implementation, so that the correction
  method is isolated from the aligner. Pass rates were 94.7% on HISAT2 and 96.7%
  on rustar.
- **The graph's cost is included but not itemised.** It turns 205,099 uniquely
  aligned reads into multimappers (2.39% of linear-unique); the MAPQ filter drops
  all of them, so the table already charges the graph for it.

## STAR

STAR has this built in (`--waspOutputMode SAMtag --varVCFfile het.vcf
--outSAMattributes NH HI AS nM vW`, then keep `vW:i:1`), which would replace
`wasp_filter.py` with the reference implementation. It was not used here: the
arm64 macOS build of STAR 2.7.11b reads zero reads from any FASTQ, including its
own test data, while still exiting 0 and logging "finished successfully". Check
`Number of input reads` in `Log.final.out` before trusting a STAR run on Apple
silicon.
