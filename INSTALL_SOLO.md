# Install and run HISAT2-Solo

Solo lives on `solo/gene-model` in `iandriver/hisat2`. The Rust graph builder
lives in `iandriver/hisat2-solo-analysis` on `main`. They are separate source
trees: the builder writes HISAT2 index files, and Solo reads those files.
No merge is needed to use them together.

## 1. Build Solo

On Linux or macOS, install a C++17 compiler, GNU Make, zlib development headers,
Python 3 and Perl. On macOS, the compiler and Make come with the Xcode command
line tools. On Debian/Ubuntu the package names are `build-essential zlib1g-dev
python3 perl`. Samtools is needed below to index your reference FASTA.

```sh
git clone --branch solo/gene-model --single-branch https://github.com/iandriver/hisat2.git
cd hisat2
make -j 2
make install PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
hisat2 --version
hisat2 --help | grep -- --solo-out-dir
bash tests/run_tests.sh
bash tests/install_smoke.sh
```

Add the PATH line to your shell startup file for future sessions. Keep wrappers
and binaries together; copying only `hisat2` is insufficient. An existing
system installation can hide this fork: check `command -v hisat2`.

CMake is an alternative to Make:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build -j 2
cmake --install build
bash tests/install_smoke.sh "$HOME/.local/bin"
```

The installation smoke test uses generated reads with known alleles. It checks
the installed wrappers, gene/allele output and automatic helper lookup from an
unrelated directory with spaces and an apostrophe in the installation path.
The helper lookup test uses a spy, so it needs no NumPy and does not validate
the EmptyDrops statistical algorithm.
The current check results and their limits are recorded in
[SOLO_INSTALL_VALIDATION.md](docs/SOLO_INSTALL_VALIDATION.md).

## 2. Count a 10x library

You need a HISAT2 index, a GTF from the same assembly, the matching reference
FASTA and a chemistry-specific barcode whitelist. An ordinary linear index
works for gene counts; allele counting needs an index containing variants.

```sh
samtools faidx genome.fa
hisat2_extract_genes.py --fai genome.fa.fai -o genes.ht2gm annotation.gtf
hisat2 -x /path/to/index/genome \
  -1 cDNA_R2.fastq.gz -2 barcode_R1.fastq.gz --solo-barcode-mate 2 \
  --solo-cb-whitelist 3M-february-2018.txt \
  --solo-cb-len 16 --solo-umi-len 12 \
  --gene-annotation genes.ht2gm --gene-strand Reverse \
  --gene-feature Gene,GeneFull --solo-out-dir Solo.out \
  -p 4 -S /dev/null
```

This geometry is for 10x 3' v3. Choose barcode/UMI lengths and strand for your
library. `Gene` counts exons; `GeneFull` counts gene bodies, including introns,
which is useful for nuclei. Check `Solo.out/Gene/Summary.csv` and the filtered
matrix. A low assignment rate can mean the wrong strand; overlapping genes in
an unfiltered annotation can cause reads to be discarded as ambiguous.
See [README_SOLO.md](README_SOLO.md) for options and limitations.

The default cell filter has no Python package dependency. `EmptyDrops_CR`
requires NumPy in the `python3` environment on PATH. The `hisat2` wrapper
locates its installed helper automatically, including when run from another
directory. If NumPy is missing, it prints the error and keeps the knee-filtered
matrix. You can also refine a matrix explicitly:

```sh
hisat2_solo_filter.py --raw Solo.out/Gene/raw --expect-cells 3000
```

For example, install NumPy in an isolated environment and activate it before
running HISAT2:

```sh
python3 -m venv "$HOME/.venvs/hisat2-solo"
. "$HOME/.venvs/hisat2-solo/bin/activate"
python3 -m pip install numpy
```

Repeat the explicit refinement command for
`GeneFull/raw` if you also want to refine that feature's cell calls.
Summary statistics and clone-size reports describe the initial knee calls;
EmptyDrops subsequently rewrites the filtered matrix. Its cell set may differ
from those reports. If invoking `hisat2-align-s` or `hisat2-align-l` directly,
set `HISAT2_SOLO_FILTER` to the absolute helper path; prefer the wrapper.

## 3. Build a custom graph index with Rust

You can use an existing HISAT2 index without installing this builder. For a
custom graph, install a current stable Rust toolchain with Cargo, then:

```sh
git clone --branch main https://github.com/iandriver/hisat2-solo-analysis.git
cd hisat2-solo-analysis
cargo build --release --manifest-path design/rust/ht2fmt/Cargo.toml --bin ht2wg
mkdir -p "$HOME/.local/bin"
install -m 755 design/rust/ht2fmt/target/release/ht2wg "$HOME/.local/bin/ht2wg"
export PATH="$HOME/.local/bin:$PATH"
```

The builder takes a FASTA plus HISAT2-format `.snp` and `.haplotype` files,
not a VCF directly. The Solo checkout includes
`hisat2_extract_snps_haplotypes_VCF.py`; check its help and inspect the outputs
when preparing a panel. The existing extractor assumes diploid calls, so a
panel with haploid calls, such as male non-PAR chrX, needs separate handling.

For a whole-human-genome build:

```sh
mkdir -p out
HT2_LARGE=1 HT2_THREADS=8 \
  ht2wg genome.fa genome.snp genome.haplotype scratchdir out/genome 134217728
```

Plan for about 19 GB RAM at these settings and more than 325 GB free scratch
space, plus the inputs and final index. The recorded unannotated GRCh38 run
took 17.7 hours on an external SSD; RAM was sampled at roughly 11 GB, not
measured as a continuous peak. Rerun the same command and inputs to resume.
Use `hisat2 -x out/genome` to align against the resulting `.ht2l` files.

For annotation-aware construction, extract junctions and exons from the GTF:

```sh
hisat2_extract_splice_sites.py annotation.gtf > genome.ss
hisat2_extract_exons.py annotation.gtf > genome.exon
HT2_LARGE=1 HT2_THREADS=8 HT2_SS=genome.ss HT2_EXON=genome.exon \
  ht2wg genome.fa genome.snp genome.haplotype scratch-annotated out/genome-annotated 134217728
```

Build-time annotation helps splice-aware alignment; `genes.ht2gm` separately
controls gene assignment during Solo counting. A sidecar gene model alone does
not add junctions to the index. Annotation-aware byte equality is documented
on the example reference and chr22; do not read that as a verified whole-genome
annotated build. Builder design and validation are in the companion repository's
`design/rust/README.md`.
