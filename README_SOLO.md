# HISAT2 single-cell (Solo)

In-aligner single-cell quantification: barcode and UMI handling, gene
assignment, UMI deduplication and cell × gene matrices, in one pass.

The reason to use it over STARsolo is memory and variant awareness. A full
10M-read run takes **52 s at 4.85 GB peak** against STARsolo's 28.3 GB on the
same data, and because HISAT2 aligns to a graph index containing known SNPs,
it can emit per-cell allele-specific counts that a linear-reference aligner
cannot produce without bias.

## Quick start

```bash
# 1. Build the gene model once per annotation (seconds, no index rebuild)
hisat2_extract_genes.py --fai genome.fa.fai -o genes.ht2gm annotation.gtf

# 2. Align and count. For 10x, pass cDNA first and the barcode read second.
hisat2 -x genome_index \
    -1 cDNA_R2.fq.gz -2 barcode_R1.fq.gz --solo-barcode-mate 2 \
    --solo-cb-whitelist 3M-february-2018.txt \
    --gene-annotation genes.ht2gm --gene-strand Reverse \
    --solo-out-dir Solo.out \
    -p 16 -S /dev/null
```

The gene model is a sidecar, not part of the index, so annotation can be
updated without the ~200 GB graph-index rebuild.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--gene-annotation <f>` | — | `.ht2gm` gene model; enables `GX:Z`/`GN:Z` tags |
| `--gene-feature <s>` | `Gene` | `Gene` (exonic) or `GeneFull` (gene body, for nuclei) |
| `--gene-strand <s>` | `Unstranded` | `Unstranded`, `Forward`, `Reverse` |
| `--solo-barcode-mate <1\|2>` | `2` | Which mate carries CB+UMI |
| `--solo-cb-in-readname` | off | Read CB/UMI from a `name_CB_UMI` suffix instead |
| `--solo-cb-whitelist <f>` | — | Barcode whitelist; required for matrix output |
| `--solo-cb-start` / `--solo-cb-len` | `1` / `16` | Barcode offset and length within the barcode read |
| `--solo-umi-start` / `--solo-umi-len` | `17` / `12` | UMI offset and length |
| `--solo-cb-match <s>` | `1MM` | `Exact` or `1MM` |
| `--solo-emit-raw` | off | Also emit `CR:Z` (barcode as sequenced, before correction) |
| `--solo-out-dir <d>` | — | Write matrices here; enables counting |
| `--solo-umi-dedup <s>` | `1MM_CR` | `Exact`, `1MM_CR`, `1MM_All`, `NoDedup` |
| `--solo-cell-filter <s>` | `CellRanger2.2` | `CellRanger2.2`, `TopCells`, `EmptyDrops_CR`, `None` |
| `--solo-expected-cells <n>` | `3000` | Expected cell count, for the knee filter |
| `--solo-top-cells <n>` | `3000` | How many barcodes `TopCells` keeps |
| `--solo-multi-mappers <s>` | `Unique` | `Unique`, `Uniform`, `EM` |
| `--solo-velocyto` | off | Also emit spliced/unspliced/ambiguous matrices |
| `--solo-allelic` | off | Per-cell REF/ALT counts (needs a SNP index) |

## Output

```
Solo.out/
  Gene/  (or GeneFull/)
    Summary.csv
    raw/       barcodes.tsv  features.tsv  matrix.mtx
               UniqueAndMult-{Uniform,EM}.mtx     # with --solo-multi-mappers
    filtered/  barcodes.tsv  features.tsv  matrix.mtx
  Velocyto/raw/  spliced.mtx  unspliced.mtx  ambiguous.mtx
  Allelic/raw/   ref.mtx  alt.mtx            # features are rsIDs
```

Layout matches STARsolo, so `Read10X` and `scanpy.read_10x_mtx` work unchanged.
Multimapper matrices are MatrixMarket `real`, since distributing a molecule
across genes gives fractional counts.

## Read this before your first run

**Get `--gene-strand` right.** A wrong setting is not an error; it produces a
well-formed matrix containing a small fraction of your reads. On the 5' GEM-X
data used below, `Reverse` gives 46.0% of reads assigned to a gene and
`Forward` gives 3.8%. HISAT2 warns when the rate falls below 10%, but check it
anyway. The correct value depends on chemistry and on which mate you pass
first, so determine it empirically on a subset.

**Use `--solo-multi-mappers EM` if you care about gene families.** The default
discards reads compatible with more than one gene, which for near-identical
paralogues means discarding most of the signal. Measured against STARsolo:

| Family | Default | With EM | STARsolo |
|---|---|---|---|
| Rpl10 + Rpl10-ps3 | 1,469 (0.15×) | 12,877 (1.29×) | 9,985 |
| Rps27 + Rps27rt | 1,044 (0.10×) | 12,928 (1.23×) | 10,532 |
| Tmsb10 + Tmsb10b | 11,553 (0.49×) | 25,625 (1.08×) | 23,829 |

Ribosomal proteins, haemoglobins and other duplicated families are affected.
EM recovers about 10% more molecules overall (395,354 on a 10M-read run).

## Concordance with STARsolo

Both tools were run on the same 10M reads from a 5' GEM-X mouse PBMC sample
against GRCm39. Reproduce with `tests/solo_concordance.py`.

| Measure | Result |
|---|---|
| **Per-cell total UMI** | **r = 0.99993** |
| **Called-cell overlap** | **Jaccard 0.99620** (3,673 of 3,684 / 3,676) |
| Per-gene total UMI | r = 0.97749 |
| Per (cell, gene), counts ≥ 10 | r = 0.94969 |
| Valid barcodes | 0.8791 vs 0.8893 |
| Sequencing saturation | 0.1476 vs 0.1467 |
| Fraction of UMIs in cells | 0.8925 vs 0.8935 |
| Total UMIs | 3,924,469 vs 4,067,946 (0.965×) |

Equality is not the target: the two align to different indexes with different
scoring, so ~5% of reads that STAR maps, HISAT2 does not. What matters is that
the disagreements are accounted for.

**On the low whole-matrix correlation.** Pearson r over all pairs is 0.855 in
log space, which looks poor until stratified: 2,024,446 of the 2,585,000 pairs
hold a single UMI, where correlation on a near-constant value is meaningless.
Those singletons are present in *both* matrices 94.1% of the time. Agreement
rises monotonically with expression — r = 0.50 at 2–3 UMIs, 0.70 at 4–9, 0.95
at ≥ 10 — which is the expected shape for sparse count data, not a defect.

**On per-gene differences.** The median per-gene ratio is 1.0036 and 95.8% of
genes with ≥ 100 UMIs fall within 10% of the overall depth ratio. The outliers
are all paralogue pairs, and are the multimapper effect described above rather
than misassignment.

**On cell calling.** The 14 disagreeing barcodes all carry 373–396 UMIs against
a median of 760 for called cells: they sit exactly at the knee, where a small
difference in totals tips a barcode either way.

## Variant-aware single-cell

With a SNP-aware index (`hisat2-build --snp`), `--solo-allelic` writes per-cell
reference and alternate counts keyed by rsID. Reads carrying alternate alleles
align to the graph without penalty, which a linear reference cannot do.
Measured on reads constructed to carry known alleles:

| | ALT reads, zero-mismatch | Total alignment score |
|---|---|---|
| Graph index | 120/120 | 0 |
| Linear index | 0/120 | −720 |

Under a stricter score threshold — the situation for a read carrying several
variants, since HISAT2's default minimum for a 100 bp read admits at most
three — the linear reference loses **every** alternate-allele read while
retaining all reference reads. That is reference bias, and it is what this
avoids.

## Limitations

- Cell calling: the knee filter has perfect precision against CellRanger's list
  (3,684/3,684) with recall 0.955. `EmptyDrops_CR` raises recall to 0.962 but
  drops precision to 0.991, so it is a trade rather than an improvement; it also
  needs numpy and runs for several minutes.
- `--solo-allelic` and Velocyto require an unambiguous assignment and skip
  reads that do not have one.
- Velocyto totals fall ~1% short of `GeneFull` because ambiguous-barcode reads
  are resolved for the gene matrix but not for Velocyto.
- Barcodes with several equally-good whitelist neighbours are resolved against
  observed counts at the end of the run; this recovers most but not all of the
  gap to CellRanger's valid-barcode rate (0.8791 vs 0.8850).
- Only `CB_UMI_Simple`-style geometry is supported. There is no equivalent of
  `CB_UMI_Complex` or `SmartSeq`, and no third-input-stream mode.
