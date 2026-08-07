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

echo "== 7. paired and unpaired modes both work =="
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
