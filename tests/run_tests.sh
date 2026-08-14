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

    # Gene and GeneFull from a single alignment pass.
    #
    # STARsolo and rustar take --soloFeatures Gene GeneFull and produce both
    # from one pass; doing it in two costs ~1.9x here (measured on 111,600
    # reads), because everything expensive -- alignment, CIGAR, reference
    # blocks -- is shared and only the interval query repeats.
    #
    # The check that matters is not that both directories appear but that each
    # matrix is byte-identical to the one a single-feature run produces. A
    # shared pass that quietly counted a read under the wrong feature, or let
    # one feature's state leak into the other, would still produce two
    # plausible matrices.
    # The generated model is one gene whose single exon spans the whole
    # reference, so Gene and GeneFull are equal there by construction and every
    # check below would pass even if the two features were swapped. Punch an
    # intron through the exon union, inside the ~11 kb the reads occupy, so the
    # features actually differ (69 vs 95 UMIs as written).
    awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} $1=="G"{print;next}
                $1=="E"{print "E",$2,$3,0,5000; print "E",$2,$3,8000,1000000; next} {print}' \
        "$SC/model.ht2gm" > "$SC/model_intron.ht2gm"
    for spec in Gene GeneFull Gene,GeneFull GeneFull,Gene; do
        ./hisat2 -x $IDX -U "$SC/reads.fq" --solo-cb-in-readname \
            --solo-cb-whitelist "$SC/whitelist.txt" --gene-annotation "$SC/model_intron.ht2gm" \
            --gene-strand Unstranded --gene-feature "$spec" \
            --solo-out-dir "$SC/ff_$(echo "$spec" | tr ',' '_')" \
            -p 4 -S /dev/null > /dev/null 2>&1
    done
    if [ -d "$SC/ff_Gene_GeneFull/Gene" ] && [ -d "$SC/ff_Gene_GeneFull/GeneFull" ]; then
        ok "one pass writes both Gene and GeneFull"
    else
        bad "--gene-feature Gene,GeneFull did not write both features"
    fi
    ffdiff=0
    for f in Gene GeneFull; do
        for m in matrix.mtx barcodes.tsv features.tsv; do
            cmp -s "$SC/ff_$f/$f/raw/$m" "$SC/ff_Gene_GeneFull/$f/raw/$m" || ffdiff=$((ffdiff+1))
        done
    done
    if [ $ffdiff -eq 0 ]; then
        ok "combined-pass matrices are identical to single-feature runs"
    else
        bad "$ffdiff combined-pass outputs differ from their single-feature run"
    fi
    if cmp -s "$SC/ff_Gene_GeneFull/Gene/raw/matrix.mtx" \
              "$SC/ff_GeneFull_Gene/Gene/raw/matrix.mtx"; then
        ok "output does not depend on the order features are listed"
    else
        bad "--gene-feature Gene,GeneFull differs from GeneFull,Gene"
    fi
    # sum(Gene) < sum(GeneFull), strictly: a read inside the exon union is
    # inside the body, and with the intron above some reads are body-only.
    # Requiring strict inequality is what keeps the comparison above honest --
    # if these two ever came out equal the identity checks would be vacuous.
    g=$(awk 'NR>3{s+=$3} END{print s+0}' "$SC/ff_Gene_GeneFull/Gene/raw/matrix.mtx")
    gf=$(awk 'NR>3{s+=$3} END{print s+0}' "$SC/ff_Gene_GeneFull/GeneFull/raw/matrix.mtx")
    if [ "$g" -lt "$gf" ]; then
        ok "sum(Gene) < sum(GeneFull) ($g < $gf), so the features are distinguishable"
    else
        bad "Gene and GeneFull agree ($g vs $gf) -- the identity checks above prove nothing"
    fi
else
    bad "could not generate single-cell test data"
fi

echo "== 9. hisat2_filter_snps.py masks pseudogenes but spares real exons =="
# A processed pseudogene lying across a real gene's exon (the DHFRP2/HLA-B case)
# must not cost that exon its variants.
cat > "$TMP/f.gtf" <<'EOF'
1	t	gene	100	200	.	+	.	gene_id "G1"; gene_biotype "protein_coding";
1	t	exon	100	200	.	+	.	gene_id "G1"; gene_biotype "protein_coding";
1	t	gene	500	600	.	+	.	gene_id "P1"; gene_biotype "processed_pseudogene";
1	t	gene	480	560	.	+	.	gene_id "G2"; gene_biotype "protein_coding";
1	t	exon	480	520	.	+	.	gene_id "G2"; gene_biotype "protein_coding";
2	t	gene	1000	1100	.	-	.	gene_id "P2"; gene_biotype "transcribed_processed_pseudogene";
2	t	gene	5000	5100	.	+	.	gene_id "U1"; gene_biotype "unprocessed_pseudogene";
EOF
# .snp positions are 0-based; the GTF above is 1-based inclusive.
printf 'a\tsingle\t1\t150\tA\nb\tsingle\t1\t499\tC\nc\tsingle\t1\t530\tG\nd\tsingle\t1\t600\tT\ne\tdeletion\t2\t995\t6\nf\tsingle\t2\t5050\tC\n' > "$TMP/f.snp"
printf 'ht0\t1\t150\t499\ta,b\nht1\t1\t530\t530\tc\nht2\t2\t995\t995\te\n' > "$TMP/f.hap"
if python3 hisat2_filter_snps.py --gtf "$TMP/f.gtf" \
        --haplotype "$TMP/f.hap" --haplotype-out "$TMP/f.out.hap" \
        "$TMP/f.snp" "$TMP/f.out.snp" > /dev/null 2>&1; then
    kept=$(cut -f1 "$TMP/f.out.snp" | tr '\n' ' ')
    # a: outside any pseudogene.  b: inside P1 but also inside G2's exon.
    # d: at 0-based 600 == 1-based 601, one past P1.  f: unprocessed, not masked.
    # c: inside P1, only in G2's intron.  e: deletion spanning into P2.
    if [ "$kept" = "a b d f " ]; then
        ok "keeps protected-exon and outside variants, drops pseudogene-only ones"
    else
        bad "filter kept '$kept', expected 'a b d f '"
    fi
    # --protect '' must drop b as well
    python3 hisat2_filter_snps.py --gtf "$TMP/f.gtf" --protect '' \
        "$TMP/f.snp" "$TMP/f.np.snp" > /dev/null 2>&1
    np=$(cut -f1 "$TMP/f.np.snp" | tr '\n' ' ')
    if [ "$np" = "a d f " ]; then ok "--protect '' drops the overlapped exon too"
    else bad "--protect '' kept '$np', expected 'a d f '"; fi
    # haplotypes: ht0 keeps a,b; ht1 loses its only variant; ht2 likewise
    hl=$(wc -l < "$TMP/f.out.hap" | tr -d ' ')
    hv=$(cut -f5 "$TMP/f.out.hap" | tr '\n' ' ')
    if [ "$hl" = "1" ] && [ "$hv" = "a,b " ]; then
        ok "haplotypes pruned to surviving variants"
    else
        bad "haplotype output was $hl line(s), variants '$hv'"
    fi
    # every id named by a haplotype must still exist in the .snp file
    miss=0
    for id in $(cut -f5 "$TMP/f.out.hap" | tr ',' '\n'); do
        cut -f1 "$TMP/f.out.snp" | grep -qx "$id" || miss=$((miss+1))
    done
    if [ $miss -eq 0 ]; then ok "no haplotype references a dropped variant"
    else bad "$miss haplotype ids missing from the .snp file"; fi
else
    bad "hisat2_filter_snps.py failed to run"
fi

echo "== 10. paired and unpaired modes both work =="
./hisat2 -x $IDX -f -U $R1 -p 1 --seed 0 -S /dev/null --summary-file "$TMP/u.txt" > /dev/null 2>&1
urate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate/\1/p' "$TMP/u.txt")
if awk -v r="${urate:-0}" 'BEGIN { exit !(r+0 >= 90) }'; then
    ok "unpaired alignment: ${urate}%"
else
    bad "unpaired alignment: ${urate}%"
fi

echo "== 11. the build reacts to header and flag changes =="
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
