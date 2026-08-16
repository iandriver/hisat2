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

echo "== 8. paired and unpaired modes both work =="
./hisat2 -x $IDX -f -U $R1 -p 1 --seed 0 -S /dev/null --summary-file "$TMP/u.txt" > /dev/null 2>&1
urate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate/\1/p' "$TMP/u.txt")
if awk -v r="${urate:-0}" 'BEGIN { exit !(r+0 >= 90) }'; then
    ok "unpaired alignment: ${urate}%"
else
    bad "unpaired alignment: ${urate}%"
fi

echo "== 9. hisat2_index_probe.py predicts what a build will do =="
# The probe exists because both of hisat2-build's failure modes are cheap to
# foresee and expensive to discover: overrunning the 2^32 node bound, and
# silently dropping variants from over-budget local graphs.
#
# The bounds it reads must stay in step with hier_idx_common.h -- if those
# constants change and the probe does not, it keeps answering confidently and
# wrongly, which is worse than not having it.
PSNP="$TMP/probe.snp"; PHAP="$TMP/probe.haplotype"; PFAI="$TMP/probe.fai"
printf 'chr1\t1000000\t0\t60\t61\n' > "$PFAI"
: > "$PSNP"; : > "$PHAP"
# 600 variants inside one 57,344 bp window: over any plausible capacity, so the
# probe must call it at risk. 20 more spread across the rest, which are not.
for i in $(seq 1 600); do
    printf 'rs%d\tsingle\tchr1\t%d\tA\n' "$i" "$((1000 + i * 20))" >> "$PSNP"
    printf 'ht%d\tchr1\t%d\t%d\trs%d\n' "$i" "$((1000 + i * 20))" "$((1000 + i * 20))" "$i" >> "$PHAP"
done
for i in $(seq 601 620); do
    printf 'rs%d\tsingle\tchr1\t%d\tA\n' "$i" "$((200000 + i * 500))" >> "$PSNP"
    printf 'ht%d\tchr1\t%d\t%d\trs%d\n' "$i" "$((200000 + i * 500))" "$((200000 + i * 500))" "$i" >> "$PHAP"
done
if python3 hisat2_index_probe.py --fai "$PFAI" --snp "$PSNP" --haplotype "$PHAP" \
        > "$TMP/probe.out" 2>&1; then
    ok "probe runs and reports a verdict"
else
    bad "probe exited non-zero on a small genome that should fit"
fi
if grep -q "windows over the budget      1" "$TMP/probe.out"; then
    ok "probe flags the one over-dense window and not the sparse ones"
else
    bad "probe did not identify the crowded window: $(grep -c . "$TMP/probe.out") lines"
fi
# A 1 Mb reference is 0.02% of the 32-bit ceiling, so it must not be called at risk.
if grep -q "fits -- chr1 built comfortably" "$TMP/probe.out"; then
    ok "probe passes a small reference on the global bound"
else
    bad "probe misjudged the global bound on a 1 Mb reference"
fi
# The whole point is that these track the builder. Compare against the header.
for c in "LOCAL_INDEX_SIZE:local_index_size" "LOCAL_MAX_GBWT:local_max_gbwt"; do
    pv=$(sed -n "s/^${c%%:*} = \(.*\)  *#.*/\1/p" hisat2_index_probe.py | head -1 | tr -d ' ')
    hv=$(sed -n "s/^static const uint32_t ${c##*:} = \(.*\);.*/\1/p" hier_idx_common.h | tr -d ' ')
    if [ -n "$pv" ] && [ "$pv" = "$hv" ]; then
        ok "${c%%:*} matches hier_idx_common.h ($hv)"
    else
        bad "${c%%:*} is '$pv' but hier_idx_common.h says '$hv'"
    fi
done

echo "== 10. the build reacts to header and flag changes =="
# Checked statically rather than by touching files: `touch Makefile` would force
# every object of all seven targets to recompile for whoever runs the suite next.
#
# Both halves have failed silently before. Objects that do not depend on their
# headers keep a stale binary through a source edit -- which can make an A/B
# comparison of two builds compare a binary against itself. Objects that do not
# depend on the Makefile survive a change to the index width, -fsigned-char or
# the assertion level, leaving a binary built half one way and half the other.
MK="$ROOT/Makefile"
objrules=$(grep -cE '^\.(obj/\$\(1\)|ht2lib-obj[a-z-]*)/%\.o:' "$MK")
withmk=$(grep -cE '^\.(obj/\$\(1\)|ht2lib-obj[a-z-]*)/%\.o:.*[^a-zA-Z]Makefile$' "$MK")
if [ "$objrules" -gt 0 ] && [ "$objrules" = "$withmk" ]; then
    ok "all $objrules object rules depend on the Makefile"
else
    bad "only $withmk of $objrules object rules depend on the Makefile"
fi
mmd=$(grep -E '\$\(CXX\)|\$\$\(CXX\)' "$MK" | grep -c '\-MMD -MP')
if [ "$mmd" = "$objrules" ]; then
    ok "all $objrules object rules generate .d files"
else
    bad "$mmd of $objrules compile rules generate .d files"
fi
if [ "$(grep -c '^-include .*\.o=\.d\|^-include \$\$(\$(1)_OBJS:\.o=\.d)' "$MK")" -ge 2 ]; then
    ok "generated .d files are included"
else
    bad "generated .d files are not included, so they have no effect"
fi

echo
echo "$pass passed, $fail failed"
[ $fail -eq 0 ]
