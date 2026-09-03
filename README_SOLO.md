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

### Counting both features in one pass

`--gene-feature Gene,GeneFull` writes both matrices from a single alignment
pass, which is what STARsolo's `--soloFeatures Gene GeneFull` does. Everything
expensive about a read -- alignment, CIGAR, reference blocks -- is independent
of the feature, so the second one costs only an extra interval query per read.
Measured on 10M mouse reads at `-p 16`: **52.9 s and 5.42 GB for both**, against
51.2 s + 51.6 s for two separate runs. Output is byte-identical to running each
feature on its own.

`GX:Z`/`GN:Z` describe one assignment, so they follow the feature listed first.
Order does not otherwise matter: the matrices are identical either way.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--gene-annotation <f>` | — | `.ht2gm` gene model; enables `GX:Z`/`GN:Z` tags |
| `--gene-feature <s>` | `Gene` | `Gene` (exonic), `GeneFull` (gene body, for nuclei), or both as a comma-separated list |
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
| `--solo-internal-priming` | off | Flag molecules primed at A-rich genomic sites (costs ~60% wall time) |
| `--solo-allelic` | off | Per-cell REF/ALT counts (needs a SNP index) |

## Output

```
Solo.out/
  Gene/  and/or  GeneFull/
    Summary.csv
    UMIcloneSize.tsv
    UMIcloneSizeByExpr.tsv
    InternalPriming.tsv                    # with --solo-internal-priming
    raw/       barcodes.tsv  features.tsv  matrix.mtx
               UniqueAndMult-{Uniform,EM}.mtx     # with --solo-multi-mappers
    filtered/  barcodes.tsv  features.tsv  matrix.mtx
  Velocyto/raw/  spliced.mtx  unspliced.mtx  ambiguous.mtx
  Allelic/raw/   ref.mtx  alt.mtx            # features are rsIDs
```

Layout matches STARsolo, so `Read10X` and `scanpy.read_10x_mtx` work unchanged.
Multimapper matrices are MatrixMarket `real`, since distributing a molecule
across genes gives fractional counts.

### Reads per UMI

`UMIcloneSize.tsv` is the distribution of reads per surviving molecule: one
row per clone size, counted over every barcode and over called cells alone.
Reads absorbed by error correction are credited to the molecule that absorbed
them, so the column totals reconcile exactly with `Total UMIs` and with
`Sequencing Saturation`.

Saturation is the mean of this distribution and nothing more. The reason to
keep the shape is the tail. A UMI is supposed to be attached during reverse
transcription, but residual UMI-bearing oligos can reprime during
preamplification PCR and put a second, valid UMI on a molecule that already
carries one. Because that UMI is intact rather than corrupted, no
deduplication mode here can see it: `1MM_CR` and `1MM_All` correct *errors* in
a UMI, and one of these is not an error. What it does leave is a heavier tail
than depth alone produces, which is what the histogram exposes and what
`UMI Clone Size P99` and `UMI Clone Size Max` in `Summary.csv` summarise.

The two columns are kept apart because ambient barcodes are overwhelmingly
single-read molecules and would bury the tail if pooled with real cells. Use
`umis_in_cells` unless you are specifically looking at the background.

`UMIcloneSizeByExpr.tsv` cuts the called-cell distribution ten ways by gene
expression, where `d1` holds the genes accounting for the first tenth of all
molecules. The deciles are cut by molecule mass rather than by gene count,
since a tenth of the genes would hold almost none of the molecules. This
answers the first objection anyone will raise about a heavy tail: a tail
confined to `d1` is a handful of very highly expressed genes behaving
differently, while one present across every decile is a property of the
chemistry.

This is a diagnostic, not a correction. Sugino & Lee report that model-based
correction removes the average count inflation but does not recover distorted
fold-changes, so a heavy tail here is a reason to treat differential
expression from that library with suspicion, not a thing to subtract off. The
branching model that turns the shape into a contamination estimate is in
*The Phantom of the PCR* (bioRxiv 2026, doi:10.64898/2026.08.01.742199).

## Read this before your first run

**Get `--gene-strand` right.** A wrong setting is not an error; it produces a
well-formed matrix containing a small fraction of your reads. On the 5' GEM-X
data used below, `Reverse` gives 46.0% of reads assigned to a gene and
`Forward` gives 3.8%. HISAT2 warns when the rate falls below 10%, but check it
anyway. The correct value depends on chemistry and on which mate you pass
first, so determine it empirically on a subset.

**Filter your annotation, or lose ~10% of your reads.** A gene model built
straight from a full GENCODE GTF makes reads ambiguous wherever annotated genes
overlap — readthrough transcripts, nested lncRNAs, and the like — and the
default `--solo-multi-mappers Unique` then discards them. Measured on 10x PBMC
1k v3 (66.6M read pairs, GRCh38):

| gene model | reads multi-gene, discarded | reads in a unique gene | total UMIs |
|---|---|---|---|
| full GENCODE v50 (78,941 genes) | **10.04%** | 44.10% | 9,130,654 |
| CellRanger gene set (32,364 genes) | **1.58%** | 48.24% | 9,558,062 |

This is why CellRanger and similar pipelines ship a *filtered* reference rather
than the raw GENCODE release. Restricting the GTF to a curated gene set before
running `hisat2_extract_genes.py` costs nothing and recovers those reads; the
alternative is `--solo-multi-mappers EM`, which redistributes them instead.

Like a wrong `--gene-strand`, this failure is silent: the matrix is well formed,
the run reports no error, and only the `Reads Multi-Gene` line in `Summary.csv`
shows what happened. Check it.

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

### Internal priming

An oligo(dT) primer is meant to find the poly(A) tail. It will also anneal to any
A-rich stretch inside a transcript, and the molecule that results is
indistinguishable from a real one downstream: it carries a barcode and a UMI, it
maps, and it is counted. The only evidence is the genome immediately past where
the read ends, which is A-rich for a mispriming event and ordinary for a real 3'
end, whose A's are on the transcript rather than in the reference.

`--solo-internal-priming` reads 20 nt of genome outward from each alignment's 3'
end -- rightward counting A for a forward alignment, leftward counting T for a
reverse one -- and calls it primed at 12 or more, or a run of 6. These are the
thresholds CellRanger, scAPA and Sierra converged on. The flag adds
`Reads Internally Primed` (over reads that reach a count, not over all reads) and
molecule-level counts to `Summary.csv`, and writes `InternalPriming.tsv`:
molecules, primed molecules and the fraction per gene, ordered by molecule count.

Nothing is filtered. The rate is reported so it can be looked at, because a
threshold that silently removed molecules would be the same mistake as counting
them without knowing.

The check is deferred until a read is actually counted, since reading the genome
at an arbitrary position is a random access into three gigabytes and better than
half of all reads never reach a count.

| library | reads primed | interpretation |
|---|---|---|
| human brain snRNA-seq (GSE163577, SRR13278449) | **17.42%** | 4.8x background |
| human PBMC 5k 3' v3, whole cell | **3.34%** | indistinguishable from background |
| random genomic windows, same rule, strand-averaged | 3.64% | the floor |

The background row is what the rule scores on 236,153 random 20 nt windows drawn
from the same `genome.fa` the index was built on. A whole-cell 3' library sits on
that floor; the single-nucleus library sits at nearly five times it. The rule is
finding structure, not firing on ordinary sequence -- and mispriming here is a
property of nuclear RNA, which is mostly intronic and where A-rich stretches
live.

Within the snRNA-seq sample the same split shows up between features:

| | molecules | primed |
|---|---|---|
| `Gene` (exonic) | 405,782 | 13.51% |
| `GeneFull` (gene body) | 1,870,831 | **17.90%** |

Roughly 19% of the 1.47 M molecules `GeneFull` adds over `Gene` are primed,
against 13.5% of the exonic ones. If you count gene bodies for nuclei, part of
what you gain is mispriming, and this says how much. Per gene, `MALAT1` sits at
22.2% and `NEAT1` at 18.2% -- both nuclear, both heavily counted in snRNA-seq.

**Cost.** On the PBMC library the flag is free: 116.9 s with it against 123.9 s
without, on 10 M reads at `-p 10`, with system time unchanged. An earlier
measurement on the snRNA-seq library, before the check was deferred and on a run
holding 10.4 GB resident, cost 64% more wall time with system time rising from
116 s to 688 s -- the signature of page faults rather than computation. Those two
numbers differ in both the code and the library, and the snRNA-seq FASTQs are no
longer on disk to separate them, so the flag stays off by default until the cost
can be re-measured on a memory-pressured run.

## Concordance with STARsolo (mouse)

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

## Concordance with CellRanger (human)

The section above is mouse. On human, the comparator is CellRanger's own
published output for 10x PBMC 1k v3 (66.6M read pairs, GRCh38, `grch38_snp`
graph index). CellRanger uses STAR internally, so this stands in for a STARsolo
run. Annotation was restricted to CellRanger's own gene ids so that the gene
model is not a variable.

| Measure | Result |
|---|---|
| **Called-cell overlap** | **Jaccard 0.9185** (1,127 shared of 1,132 / 1,222) |
| **Per-cell total UMI** | **r = 0.9950** (Spearman 0.9941) |
| Per-gene total UMI | r = 0.9550 (Spearman 0.9552) |
| Genes within 2× | 94.4% (10,505 of 11,123 with ≥ 50 UMIs) |
| Total UMIs on shared cells | 8,618,514 vs 9,130,347 (0.944×) |
| Overall alignment rate | 87.20% |

Two sources of disagreement, neither a misassignment:

**Annotation drift.** CellRanger 3.0.0 is built on GENCODE v28. Gene ids are not
stable in meaning across seven years of releases — `CAST` (ENSG00000153113) is
the clearest case, where v50 carries an lncRNA at `5:95,962,001-96,631,085` and
the protein-coding gene at `96,247,756-96,779,595`. Matching gene *ids* does not
match their *coordinates*.

**Pseudogenes inside annotated exons.** `OLFM3` scores 12,992 UMIs here against
CellRanger's 0, which is impossible for a brain-specific gene in PBMCs. The
cause is `RPSAP19`, a processed pseudogene of the highly expressed ribosomal
protein gene RPSA, lying at `1:101,786,340-101,787,219` — entirely inside
OLFM3's first exon. It is in neither reference, so reads from RPSA that land on
it have only OLFM3 to be assigned to; STAR discards them as multimappers.

This is the same mechanism behind the paralogue families in the mouse table
above, and it is the strongest argument for `--solo-multi-mappers EM` on any
sample where ribosomal-protein genes matter.

### Memory, on human

| | HISAT2-solo (graph) | STAR (linear) |
|---|---|---|
| index on disk | **6.5 GB** | 29 GB |
| peak RSS while mapping | **9.0–9.4 GB** | 31.3 GB |

Note the direction of the trade: HISAT2's index is 4.5× smaller and needs ~3.4×
less RAM, but a human *graph* index cannot practically be built — use the
prebuilt one. STAR's is larger and slower to build, and you can rebuild it for
any assembly you like.

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
