#!/usr/bin/env bash
# Exercise the installed wrappers and Solo from outside the source tree.
# Requires an existing Make build, or accepts a directory of installed binaries
# as its first argument. Copies that installation before testing it.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TMP=$(mktemp -d)
cleanup() {
    local result=$?
    if [ "$result" -ne 0 ]; then
        for log in "$TMP/install.log" "$TMP/alignment.log" "$TMP/fallback.log"; do
            [ ! -f "$log" ] || tail -20 "$log" >&2
        done
    fi
    rm -rf "$TMP"
}
trap cleanup EXIT
PREFIX="$TMP/install space's"
BIN="$PREFIX/bin"
if [ "$#" -gt 0 ]; then
    mkdir -p "$BIN"
    cp -R "$1/." "$BIN/"
else
    make -C "$ROOT" install PREFIX="$PREFIX" > "$TMP/install.log" 2>&1
fi
"$BIN/hisat2" --help > "$TMP/help.txt"
grep -q -- --solo-out-dir "$TMP/help.txt"
"$BIN/hisat2-build" --version > /dev/null
"$BIN/hisat2-inspect" --version > /dev/null
"$BIN/hisat2_extract_genes.py" --help > /dev/null
test -x "$BIN/hisat2_solo_filter.py"

python3 "$ROOT/tests/make_allelic_reads.py" \
  --fasta "$ROOT/example/reference/22_20-21M.fa" \
  --snp "$ROOT/example/reference/22_20-21M.snp" \
  --outdir "$TMP/input" --nsnps 20 > /dev/null
mkdir "$TMP/index space's"
cp "$ROOT"/example/index/22_20-21M_snp.*.ht2 "$TMP/index space's/"

# A helper spy verifies argv and path resolution independently of NumPy.
# Copy it over the helper in this disposable install, so the wrapper must
# resolve its own sibling rather than use a test-only environment override.
cat > "$BIN/hisat2_solo_filter.py" <<'PY'
import os, sys
assert os.path.isabs(__file__)
assert sys.argv[1] == '--raw'
assert os.path.isfile(os.path.join(sys.argv[2], 'matrix.mtx'))
assert sys.argv[3:] == ['--expect-cells', '3000']
with open(os.environ['SOLO_HELPER_MARKER'], 'a') as f:
    f.write(sys.argv[2] + '\n')
PY
mkdir "$TMP/unrelated"
cd "$TMP/unrelated"
export SOLO_HELPER_MARKER="$TMP/helper-ran"
unset HISAT2_SOLO_FILTER
run_solo() {
"$BIN/hisat2" -x "$TMP/index space's/22_20-21M_snp" \
  -U "$TMP/input/reads.fq" --solo-cb-in-readname \
  --solo-cb-whitelist "$TMP/input/whitelist.txt" \
  --gene-annotation "$TMP/input/model.ht2gm" --gene-strand Unstranded \
  --gene-feature Gene,GeneFull --solo-cell-filter EmptyDrops_CR \
  --solo-allelic --solo-out-dir "$1" \
  -p 2 -S /dev/null
}
run_solo "$TMP/output space's" > "$TMP/alignment.log" 2>&1
test "$(wc -l < "$SOLO_HELPER_MARKER" | tr -d ' ')" = 2
test -s "$TMP/output space's/Gene/raw/matrix.mtx"
test -s "$TMP/output space's/GeneFull/raw/matrix.mtx"
test -s "$TMP/output space's/Allelic/raw/alt.mtx"
test "$(grep -c 'EmptyDrops_CR refinement complete' "$TMP/alignment.log")" = 2

# A failing helper must leave usable knee matrices and expose its diagnostic.
cat > "$BIN/hisat2_solo_filter.py" <<'PY'
import sys
print('helper failure sentinel', file=sys.stderr)
sys.exit(2)
PY
run_solo "$TMP/fallback" > "$TMP/fallback.log" 2>&1
grep -q 'helper failure sentinel' "$TMP/fallback.log"
test "$(grep -c 'could not run EmptyDrops_CR automatically' "$TMP/fallback.log")" = 2
test -s "$TMP/fallback/Gene/filtered/matrix.mtx"
test -s "$TMP/fallback/GeneFull/filtered/matrix.mtx"
echo 'Installed wrappers, Solo matrices, helper lookup and fallback: PASS'
