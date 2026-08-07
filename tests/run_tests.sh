#!/bin/bash
#
# HISAT2 correctness tests.
#
# Runs from the repository root (or anywhere; paths are resolved relative to
# this script). Uses only the in-repo example/ data, so it needs no downloads
# and finishes in well under a minute.
#
#   tests/run_tests.sh
#
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
cd "$ROOT"

CXX="${CXX:-c++}"
IDX=example/index/22_20-21M_snp
R1=example/reads/reads_1.fa
R2=example/reads/reads_2.fa
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

pass=0; fail=0
ok()   { pass=$((pass+1)); echo "  ok    $1"; }
bad()  { fail=$((fail+1)); echo "  FAIL  $1"; }

# @PG carries the command line and @HD can vary; neither is part of the result.
canon() { grep -v '^@PG' "$1" | grep -v '^@HD'; }
hash_of() { canon "$1" | shasum | cut -d' ' -f1; }

echo "== 1. tools run =="
for t in ./hisat2 ./hisat2-build ./hisat2-inspect; do
    if $t --version > /dev/null 2>&1; then ok "$t --version"; else bad "$t --version"; fi
done

echo "== 2. alignment produces the expected rate =="
./hisat2 -x $IDX -f -1 $R1 -2 $R2 -p 1 --seed 0 --reorder \
    -S "$TMP/a.sam" --summary-file "$TMP/a.txt" > /dev/null 2>&1
rate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate/\1/p' "$TMP/a.txt")
if awk -v r="${rate:-0}" 'BEGIN { exit !(r+0 >= 90) }'; then
    ok "overall alignment rate ${rate}% (>= 90%)"
else
    bad "overall alignment rate ${rate}% (< 90%)"
fi
recs=$(grep -vc '^@' "$TMP/a.sam")
if [ "$recs" -gt 0 ]; then ok "$recs SAM records emitted"; else bad "no SAM records"; fi

echo "== 3. single-threaded output is reproducible =="
# Score ties are broken with a seeded PRNG, so a fixed --seed must reproduce.
h1=$(hash_of "$TMP/a.sam")
same=1
for i in 2 3; do
    ./hisat2 -x $IDX -f -1 $R1 -2 $R2 -p 1 --seed 0 --reorder \
        -S "$TMP/rep$i.sam" > /dev/null 2>&1
    [ "$(hash_of "$TMP/rep$i.sam")" = "$h1" ] || same=0
done
if [ $same -eq 1 ]; then ok "3 runs byte-identical"; else bad "repeat runs differ"; fi

echo "== 4. output does not depend on thread count =="
# --no-temp-splicesite is required, not incidental: novel splice sites are
# discovered in read order, so a read can only use junctions found by reads
# processed before it. With the feature on, -p 1 legitimately differs from
# -p > 1, and this test would fail for a reason unrelated to the change under
# test. See also --novel-splicesite-outfile/-infile for the two-pass form,
# which is both thread-invariant and more sensitive than either.
prev=""
same=1
for p in 1 2 4; do
    ./hisat2 -x $IDX -f -1 $R1 -2 $R2 -p $p --seed 0 --reorder --no-temp-splicesite \
        -S "$TMP/t$p.sam" > /dev/null 2>&1
    h=$(hash_of "$TMP/t$p.sam")
    [ -n "$prev" ] && [ "$h" != "$prev" ] && same=0
    prev="$h"
done
if [ $same -eq 1 ]; then ok "-p 1/2/4 byte-identical"; else bad "output varies with -p"; fi

echo "== 5. index round-trip =="
./hisat2-build -p 2 example/reference/22_20-21M.fa "$TMP/idx" > /dev/null 2>&1
n=$(ls "$TMP"/idx*.ht2 2>/dev/null | wc -l | tr -d ' ')
if [ "$n" -ge 6 ]; then ok "hisat2-build wrote $n index files"; else bad "hisat2-build produced $n files"; fi
if ./hisat2-inspect -s "$TMP/idx" 2>/dev/null | grep -q "Index version"; then
    ok "hisat2-inspect reads the index back"
else
    bad "hisat2-inspect failed"
fi
./hisat2 -x "$TMP/idx" -f -1 $R1 -2 $R2 -p 1 --seed 0 \
    -S /dev/null --summary-file "$TMP/f.txt" > /dev/null 2>&1
frate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate/\1/p' "$TMP/f.txt")
if awk -v r="${frate:-0}" 'BEGIN { exit !(r+0 >= 90) }'; then
    ok "alignment against the freshly built index: ${frate}%"
else
    bad "alignment against the freshly built index: ${frate}%"
fi

echo "== 6. gzipped input matches plain input =="
# Only meaningful when the binary was built with zlib; a WITH_ZLIB=0 build
# leaves .gz to the wrapper's named pipes, which this test does not exercise.
if ./hisat2-align-s --version 2>/dev/null | grep -q '\-DWITH_ZLIB'; then
    gzip -c $R1 > "$TMP/r1.fa.gz"
    gzip -c $R2 > "$TMP/r2.fa.gz"
    ./hisat2 -x $IDX -f -1 "$TMP/r1.fa.gz" -2 "$TMP/r2.fa.gz" -p 1 --seed 0 --reorder \
        -S "$TMP/gz.sam" > /dev/null 2>&1
    if [ "$(hash_of "$TMP/gz.sam")" = "$h1" ]; then
        ok "gzipped input is byte-identical to plain"
    else
        bad "gzipped input differs from plain"
    fi

    # A truncated file used to look like a short one: the wrapper's `gzip -dc`
    # died into a pipe and the aligner just saw EOF, reported fewer reads and
    # exited 0. Reading the stream natively means the error is not lost.
    sz=$(wc -c < "$TMP/r1.fa.gz")
    head -c $((sz / 2)) "$TMP/r1.fa.gz" > "$TMP/trunc.fa.gz"
    if ./hisat2-align-s -x $IDX -f -U "$TMP/trunc.fa.gz" -p 1 -S /dev/null > /dev/null 2>&1; then
        bad "truncated gzip input was accepted silently"
    else
        ok "truncated gzip input fails loudly"
    fi
else
    echo "  skip  built without zlib"
fi

echo "== 7. BAM output =="
# samtools is the oracle: converting our BAM back to SAM must reproduce what
# SAM mode wrote, record for record. @PG differs by construction (the command
# line contains --bam, and samtools appends its own line).
if ! ./hisat2-align-s --version 2>/dev/null | grep -q '\-DWITH_ZLIB'; then
    echo "  skip  built without zlib"
elif ! command -v samtools > /dev/null 2>&1; then
    echo "  skip  samtools not installed"
else
    ./hisat2 -x $IDX -f -1 $R1 -2 $R2 -p 1 --seed 0 --reorder --bam \
        -S "$TMP/a.bam" > /dev/null 2>&1
    if samtools quickcheck "$TMP/a.bam" 2>/dev/null; then
        ok "BAM passes samtools quickcheck (valid BGZF, EOF block present)"
    else
        bad "samtools quickcheck rejected the BAM"
    fi
    samtools view -h "$TMP/a.bam" 2>/dev/null | grep -v '^@PG' > "$TMP/rt.sam"
    grep -v '^@PG' "$TMP/a.sam" > "$TMP/ref.sam"
    if diff -q "$TMP/ref.sam" "$TMP/rt.sam" > /dev/null 2>&1; then
        ok "BAM round-trips to identical SAM ($(grep -vc '^@' "$TMP/ref.sam") records)"
    else
        bad "BAM round-trip differs from SAM"
    fi

    # Exercises what HISAT2 itself never emits: float and array tags, every
    # CIGAR operation, odd-length SEQ, '*' SEQ/QUAL, and the integer widths the
    # tag encoder switches between.
    if $CXX -std=c++17 -DWITH_ZLIB -iquote . -o "$TMP/unit_bam" \
            tests/unit_bam.cpp bam.cpp -lz > /dev/null 2>&1 \
       && "$TMP/unit_bam" "$TMP/unit.bam" > /dev/null 2>&1; then
        if samtools view "$TMP/unit.bam" > "$TMP/unit_rt.sam" 2>/dev/null; then
            # The test's own records are the expectation: every field must come
            # back byte for byte.
            if [ "$(wc -l < "$TMP/unit_rt.sam" | tr -d ' ')" = "8" ] && \
               grep -q 'Zh:f:3.5' "$TMP/unit_rt.sam" && \
               grep -q 'Zi:B:c,-1,2,-3' "$TMP/unit_rt.sam" && \
               grep -q 'Zj:B:i,100000,-100000' "$TMP/unit_rt.sam" && \
               grep -q 'Zg:i:-2147483648' "$TMP/unit_rt.sam" && \
               grep -q '1H2S3M1I1P1D3N2=1X1S1H' "$TMP/unit_rt.sam"; then
                ok "BAM encoder unit test: tags, CIGAR ops, odd SEQ, '*' fields"
            else
                bad "BAM encoder unit test produced unexpected records"
            fi
        else
            bad "samtools could not read the unit-test BAM"
        fi
    else
        bad "BAM encoder unit test failed to build or run"
    fi
fi

echo "== 8. single-cell counting =="
# Uses only example/ data: make_allelic_reads.py builds barcoded reads whose
# true REF/ALT allele is known by construction, plus a matching gene model and
# whitelist.
SC="$TMP/sc"
if python3 tests/make_allelic_reads.py --fasta example/reference/22_20-21M.fa \
        --snp example/reference/22_20-21M.snp --outdir "$SC" --nsnps 20 > /dev/null 2>&1; then
    ./hisat2 -x $IDX -U "$SC/reads.fq" --solo-cb-in-readname \
        --solo-cb-whitelist "$SC/whitelist.txt" --gene-annotation "$SC/model.ht2gm" \
        --gene-strand Unstranded --solo-out-dir "$SC/Solo.out" --solo-allelic \
        -p 4 -S /dev/null > /dev/null 2>&1
    if [ -f "$SC/Solo.out/Allelic/raw/ref.mtx" ]; then
        ok "--solo-allelic wrote per-cell REF/ALT matrices"
    else
        bad "--solo-allelic produced no Allelic output"
    fi

    # Counts must match the hand-computed expectation, and Summary.csv must
    # agree with the matrices it summarises. Those tallies are filled in while
    # the allelic matrices are written, so an ordering change can leave every
    # allelic line reading zero next to correct matrices -- which it once did.
    if python3 - "$SC" <<'PYEOF'
import sys, os
d = sys.argv[1]; D = os.path.join(d, "Solo.out", "Allelic", "raw")
feats = [l.split('\t')[0] for l in open(D + "/features.tsv")]
bcs = [l.strip() for l in open(D + "/barcodes.tsv")]
def load(p):
    m = {}
    for i, l in enumerate(open(p)):
        if i < 3: continue
        r, c, v = l.split(); m[(feats[int(r)-1], bcs[int(c)-1])] = int(v)
    return m
ref, alt = load(D + "/ref.mtx"), load(D + "/alt.mtx")
bad = 0
for i, l in enumerate(open(os.path.join(d, "expected.tsv"))):
    if i == 0: continue
    bc, rs, allele, n = l.rstrip().split('\t')
    if (ref if allele == 'ref' else alt).get((rs, bc), 0) != int(n): bad += 1
summ = {}
for l in open(os.path.join(d, "Solo.out", "Gene", "Summary.csv")):
    k, _, v = l.rstrip().partition(',')
    summ[k] = v
if bad: sys.exit(1)
if int(summ.get("Allelic UMIs: Reference", -1)) != sum(ref.values()): sys.exit(2)
if int(summ.get("Allelic UMIs: Alternate", -1)) != sum(alt.values()): sys.exit(3)
seen = set(k[0] for k in list(ref) + list(alt))
if int(summ.get("Variants Observed", -1)) != len(seen): sys.exit(4)
PYEOF
    then
        ok "allelic counts match expectation and Summary.csv"
    else
        bad "allelic counts or Summary.csv disagree (exit $?)"
    fi

    # Fractional counts are only legal in the UniqueAndMult-* matrices that
    # --solo-multi-mappers writes; everything else must stay integer.
    nonint=0
    for m in "$SC"/Solo.out/Gene/raw/matrix.mtx "$SC"/Solo.out/Allelic/raw/ref.mtx \
             "$SC"/Solo.out/Allelic/raw/alt.mtx; do
        head -1 "$m" | grep -q "coordinate integer" || nonint=$((nonint+1))
        c=$(tail -n +4 "$m" | awk '{if ($3 != int($3)) n++} END{print n+0}')
        nonint=$((nonint + c))
    done
    if [ $nonint -eq 0 ]; then ok "count matrices are integer"; else bad "$nonint non-integer matrices/values"; fi
else
    bad "could not generate single-cell test data"
fi

echo "== 9. paired and unpaired modes both work =="
./hisat2 -x $IDX -f -U $R1 -p 1 --seed 0 -S /dev/null --summary-file "$TMP/u.txt" > /dev/null 2>&1
urate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate/\1/p' "$TMP/u.txt")
if awk -v r="${urate:-0}" 'BEGIN { exit !(r+0 >= 90) }'; then
    ok "unpaired alignment: ${urate}%"
else
    bad "unpaired alignment: ${urate}%"
fi

echo
echo "$pass passed, $fail failed"
[ $fail -eq 0 ]
