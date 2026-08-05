/*
 * Unit tests for GeneModel / GeneIvIndex.
 *
 * Standalone: builds with
 *   c++ -std=c++11 -O2 -I. gene_model.cpp tests/unit_gene_model.cpp -o unit_gene_model
 *
 * The interval index is the crux of single-cell gene assignment, so it is
 * checked against a brute-force O(n) scan on randomized data rather than
 * against hand-picked cases only.
 */

#include "gene_model.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static int failures = 0;
static int checks = 0;

static void check(bool cond, const std::string& what) {
    checks++;
    if(!cond) { failures++; std::cerr << "  FAIL: " << what << "\n"; }
}

// ---------------------------------------------------------------- helpers

static std::vector<uint32_t> sortedUnique(std::vector<uint32_t> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

/** Deterministic PRNG so failures reproduce exactly. */
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
    int64_t range(int64_t lo, int64_t hi) { return lo + (int64_t)(next() % (uint64_t)(hi - lo)); }
};

// ------------------------------------------------- 1. index vs brute force

static void testIndexVsBruteForce() {
    std::cout << "1. interval index vs brute force (10k intervals, 10k queries)\n";
    Rng rng(0x9E3779B97F4A7C15ULL);

    std::vector<GeneIv> truth;
    GeneIvIndex idx;
    const int64_t GENOME = 5000000;
    for(int i = 0; i < 10000; i++) {
        GeneIv iv;
        iv.start = rng.range(0, GENOME);
        // Mostly short, but deliberately include some intervals past the
        // hoisting threshold so that path is exercised too.
        int64_t len = (i % 500 == 0) ? rng.range(1000000, 3000000) : rng.range(1, 20000);
        iv.end  = iv.start + len;
        iv.gene = (uint32_t)i;
        truth.push_back(iv);
        idx.add(iv);
    }
    idx.build();
    check(idx.size() == truth.size(), "index retains every interval");

    uint64_t maxScan = 0, totScan = 0;
    int nq = 10000;
    for(int q = 0; q < nq; q++) {
        int64_t s = rng.range(0, GENOME);
        int64_t e = s + rng.range(1, 500);

        std::vector<uint32_t> got;
        idx.query(s, e, got);
        got = sortedUnique(got);

        std::vector<uint32_t> want;
        for(size_t i = 0; i < truth.size(); i++) {
            if(truth[i].start < e && truth[i].end > s) want.push_back(truth[i].gene);
        }
        want = sortedUnique(want);

        if(got != want) {
            std::ostringstream os;
            os << "query [" << s << "," << e << ") got " << got.size()
               << " want " << want.size();
            check(false, os.str());
            if(failures > 3) return;
        } else {
            checks++;
        }
        maxScan = std::max(maxScan, idx.lastScanned);
        totScan += idx.lastScanned;
    }
    std::cout << "   scan length: mean " << (totScan / nq) << ", max " << maxScan
              << " (of " << truth.size() << " intervals)\n";
    // With hoisting, a point-ish query should touch a tiny fraction of the set.
    check(maxScan < truth.size() / 10, "hoisting keeps worst-case scan well below O(n)");
}

// ------------------------------------------------- 2. empty / edge cases

static void testEdgeCases() {
    std::cout << "2. edge cases\n";
    GeneIvIndex empty;
    empty.build();
    std::vector<uint32_t> out;
    empty.query(0, 100, out);
    check(out.empty(), "empty index yields nothing");

    GeneIvIndex one;
    GeneIv iv; iv.start = 100; iv.end = 200; iv.gene = 7;
    one.add(iv);
    one.build();

    out.clear(); one.query(150, 160, out);
    check(out.size() == 1 && out[0] == 7, "query inside interval hits");
    out.clear(); one.query(0, 100, out);
    check(out.empty(), "half-open: query ending exactly at start misses");
    out.clear(); one.query(200, 300, out);
    check(out.empty(), "half-open: query starting exactly at end misses");
    out.clear(); one.query(99, 101, out);
    check(out.size() == 1, "query straddling start hits");
    out.clear(); one.query(199, 201, out);
    check(out.size() == 1, "query straddling end hits");
    out.clear(); one.query(50, 50, out);
    check(out.empty(), "empty query range yields nothing");
}

// ------------------------------------------------- 3. assignGenes semantics

/** Writes a small hand-built model exercising the cases that break naive logic. */
static std::string writeMiniModel(const std::string& path) {
    std::ofstream o(path.c_str());
    o << "#hisat2-gene-model\tv1\n";
    o << "#ref\tchr1\t100000\n";
    // g0: two exons, + strand, intron 2000..3000
    o << "G\t0\tGENE0\tAlpha\tchr1\t+\t1000\t4000\n";
    o << "E\t0\tchr1\t1000\t2000\n";
    o << "E\t0\tchr1\t3000\t4000\n";
    o << "J\t0\tchr1\t1999\t3000\n";
    // g1: overlaps g0's second exon, opposite strand
    o << "G\t1\tGENE1\tBeta\tchr1\t-\t3500\t5000\n";
    o << "E\t1\tchr1\t3500\t5000\n";
    // g2: single exon lying entirely inside g0's intron
    o << "G\t2\tGENE2\tGamma\tchr1\t+\t2200\t2400\n";
    o << "E\t2\tchr1\t2200\t2400\n";
    o.close();
    return path;
}

static void testAssignGenes() {
    std::cout << "3. assignGenes containment semantics\n";
    std::string path = "/tmp/hisat2_unit_mini.ht2gm";
    writeMiniModel(path);

    std::vector<std::string> refnames; refnames.push_back("chr1");
    std::vector<int64_t> reflens; reflens.push_back(100000);
    GeneModel gm;
    std::string err;
    if(!gm.load(path, refnames, reflens, err)) {
        check(false, "load mini model: " + err);
        return;
    }
    check(gm.numGenes() == 3, "mini model has 3 genes");
    check(std::string(gm.geneId(0)) == "GENE0", "gene_id round-trips");
    check(std::string(gm.geneName(0)) == "Alpha", "gene_name round-trips");

    std::vector<uint32_t> out, scratch;
    std::vector<std::pair<int64_t, int64_t> > blocks;

    // A read wholly inside g0 exon 1.
    out.clear(); blocks.clear(); blocks.push_back(std::make_pair(1100, 1200));
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_UNSTRANDED, out, scratch);
    check(out.size() == 1 && out[0] == 0, "exonic read assigns to its gene");

    // A read in g0's intron: not exonic for g0, but inside its body.
    out.clear(); blocks.clear(); blocks.push_back(std::make_pair(2600, 2700));
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_UNSTRANDED, out, scratch);
    check(out.empty(), "intronic read is not counted under Gene");
    out.clear();
    gm.assignGenes(0, blocks, true, GENE_FEATURE_BODY, GENE_STRAND_UNSTRANDED, out, scratch);
    check(out.size() == 1 && out[0] == 0, "intronic read is counted under GeneFull");

    // A read straddling the exon/intron boundary belongs to neither.
    out.clear(); blocks.clear(); blocks.push_back(std::make_pair(1950, 2050));
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_UNSTRANDED, out, scratch);
    check(out.empty(), "read straddling an exon boundary is not exonic (containment, not overlap)");

    // A spliced read across g0's annotated junction: two blocks, both exonic.
    out.clear(); blocks.clear();
    blocks.push_back(std::make_pair(1900, 2000));
    blocks.push_back(std::make_pair(3000, 3100));
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_UNSTRANDED, out, scratch);
    check(out.size() == 1 && out[0] == 0, "spliced read assigns via block intersection");
    check(gm.junctionInGene(0, 1999, 3000), "annotated junction is recognised");
    check(!gm.junctionInGene(0, 1999, 3500), "unannotated junction is not");

    // A read inside g0's intron that lands on g2, which nests there.
    out.clear(); blocks.clear(); blocks.push_back(std::make_pair(2250, 2350));
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_UNSTRANDED, out, scratch);
    check(out.size() == 1 && out[0] == 2, "nested gene inside another's intron is found");

    // Overlap region of g0 exon2 and g1: ambiguous, both returned unstranded.
    out.clear(); blocks.clear(); blocks.push_back(std::make_pair(3600, 3700));
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_UNSTRANDED, out, scratch);
    check(sortedUnique(out).size() == 2, "overlapping genes both reported when unstranded");

    // ...but strandedness disambiguates them, which is the whole point of it.
    out.clear();
    gm.assignGenes(0, blocks, true, GENE_FEATURE_EXONIC, GENE_STRAND_FORWARD, out, scratch);
    check(out.size() == 1 && out[0] == 0, "forward-strand read resolves to the + gene");
    out.clear();
    gm.assignGenes(0, blocks, false, GENE_FEATURE_EXONIC, GENE_STRAND_FORWARD, out, scratch);
    check(out.size() == 1 && out[0] == 1, "reverse-strand read resolves to the - gene");
}

// ------------------------------------------------- 4. assembly-mismatch guard

static void testRefGuard() {
    std::cout << "4. assembly mismatch is a hard error\n";
    std::string path = "/tmp/hisat2_unit_mini.ht2gm";
    writeMiniModel(path);
    GeneModel gm;
    std::string err;

    std::vector<std::string> wrongName; wrongName.push_back("chrX");
    std::vector<int64_t> len1; len1.push_back(100000);
    check(!gm.load(path, wrongName, len1, err), "unknown sequence name is rejected");
    check(err.find("different assemblies") != std::string::npos,
          "error explains the assembly mismatch");

    GeneModel gm2;
    std::vector<std::string> rightName; rightName.push_back("chr1");
    std::vector<int64_t> wrongLen; wrongLen.push_back(999999);
    check(!gm2.load(path, rightName, wrongLen, err), "conflicting sequence length is rejected");
}

// ------------------------------------------------- 5. real model, if present

static void testRealModel(const char* modelPath, const char* faiPath) {
    std::cout << "5. real gene model: " << modelPath << "\n";
    std::vector<std::string> refnames;
    std::vector<int64_t> reflens;
    std::ifstream fai(faiPath);
    if(!fai.good()) { std::cout << "   (skipped: no .fai)\n"; return; }
    std::string line;
    while(std::getline(fai, line)) {
        std::istringstream ls(line);
        std::string name, len;
        std::getline(ls, name, '\t');
        std::getline(ls, len, '\t');
        refnames.push_back(name);
        reflens.push_back(strtoll(len.c_str(), NULL, 10));
    }

    GeneModel gm;
    std::string err;
    if(!gm.load(modelPath, refnames, reflens, err)) {
        std::cout << "   (skipped: " << err << ")\n";
        return;
    }
    std::cout << "   loaded " << gm.numGenes() << " genes across "
              << refnames.size() << " sequences\n";
    check(gm.numGenes() > 0, "real model loads");

    // Every gene's own body midpoint must find that gene under GeneFull.
    size_t selfHits = 0, tested = 0;
    std::vector<uint32_t> out, scratch;
    std::vector<std::pair<int64_t, int64_t> > blocks;
    for(uint32_t g = 0; g < gm.numGenes(); g += 97) {
        const GeneRec& rec = gm.gene(g);
        int64_t mid = (rec.bodyStart + rec.bodyEnd) / 2;
        out.clear(); blocks.clear();
        blocks.push_back(std::make_pair(mid, mid + 1));
        gm.assignGenes(rec.refid, blocks, rec.fw != 0, GENE_FEATURE_BODY,
                       GENE_STRAND_UNSTRANDED, out, scratch);
        tested++;
        if(std::find(out.begin(), out.end(), g) != out.end()) selfHits++;
    }
    std::cout << "   body midpoint self-hit: " << selfHits << "/" << tested << "\n";
    check(selfHits == tested, "every sampled gene is found at its own midpoint");
}

int main(int argc, char** argv) {
    testIndexVsBruteForce();
    testEdgeCases();
    testAssignGenes();
    testRefGuard();
    if(argc >= 3) testRealModel(argv[1], argv[2]);

    std::cout << "\n" << (failures ? "FAILED" : "OK") << ": "
              << (checks - failures) << "/" << checks << " checks passed\n";
    return failures ? 1 : 0;
}
