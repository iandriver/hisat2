/*
 * Unit tests for UMI collapsing and the reads-per-UMI (clone size) histogram.
 *
 *   c++ -std=c++11 -O2 -iquote . tests/unit_solo_umi.cpp -o unit_solo_umi
 *
 * The clone sizes are bookkeeping laid over deduplication, so the test that
 * matters most is that they change no decision: section 2 reimplements the
 * collapsing rules as they stood before the histogram existed and requires the
 * survivor count to agree on random blocks. A silent change here would move
 * every count in the matrix.
 */

#include "solo_umi.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

static int failures = 0, checks = 0;
static void check(bool cond, const std::string& what) {
    checks++;
    if(!cond) { failures++; std::cerr << "  FAIL: " << what << "\n"; }
}

struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
    uint32_t range(uint32_t lo, uint32_t hi) { return lo + (uint32_t)(next() % (uint64_t)(hi - lo)); }
};

/** Deduplication exactly as it was before clone sizes were recorded. */
static uint32_t referenceSurvivors(const std::vector<uint64_t>& umis,
                                   const std::vector<uint32_t>& reads,
                                   SoloUmiDedup dedup, int umiLen) {
    if(dedup == SOLO_UMI_EXACT) return (uint32_t)umis.size();
    std::vector<char> merged(umis.size(), 0);
    for(size_t a = 0; a < umis.size(); a++) {
        if(merged[a]) continue;
        for(size_t b = 0; b < umis.size(); b++) {
            if(a == b || merged[b]) continue;
            if(!soloUmiWithin1(umis[a], umis[b], umiLen)) continue;
            if(dedup == SOLO_UMI_1MM_ALL) {
                if(b > a) merged[b] = 1;
            } else {
                if(reads[b] < reads[a] || (reads[b] == reads[a] && b > a)) merged[b] = 1;
            }
        }
    }
    uint32_t n = 0;
    for(size_t a = 0; a < umis.size(); a++) if(!merged[a]) n++;
    return n;
}

// ------------------------------------------------------- pass-through modes

static void testTrivialModes() {
    std::cout << "1. Exact and NoDedup pass clone sizes through\n";
    std::vector<uint64_t> umis;
    std::vector<uint32_t> reads, clones;
    umis.push_back(0x111); reads.push_back(4);
    umis.push_back(0x112); reads.push_back(1);   // 1MM neighbour of the first
    umis.push_back(0xabc); reads.push_back(7);

    soloUmiClones(umis, reads, SOLO_UMI_EXACT, 12, clones);
    check(clones.size() == 3, "Exact keeps every distinct UMI");
    check(clones[0] == 4 && clones[1] == 1 && clones[2] == 7, "Exact clone sizes are the read counts");

    soloUmiClones(umis, reads, SOLO_UMI_NODEDUP, 12, clones);
    check(clones.size() == 3, "NoDedup keeps every distinct UMI");
}

// ------------------------------------------------------------- no drift

static void testMatchesReference(SoloUmiDedup dedup, const char* name) {
    std::cout << "2. " << name << " survivor count is unchanged\n";
    Rng rng(0x5eed1234u);
    const int umiLen = 12;
    bool agree = true;
    size_t merges = 0;
    for(int t = 0; t < 4000 && agree; t++) {
        // Seed a few UMIs, then add 1MM neighbours of them so blocks actually
        // merge; wholly random 12-mers almost never collide.
        size_t nSeed = rng.range(1, 5);
        std::vector<uint64_t> umis;
        std::vector<uint32_t> reads;
        for(size_t i = 0; i < nSeed; i++) {
            umis.push_back(rng.next() & 0xffffffull);
            reads.push_back(rng.range(1, 30));
        }
        size_t nNbr = rng.range(0, 6);
        for(size_t i = 0; i < nNbr && !umis.empty(); i++) {
            uint64_t base = umis[rng.range(0, (uint32_t)umis.size())];
            int pos = (int)rng.range(0, (uint32_t)umiLen);
            uint64_t cur = (base >> (2 * pos)) & 3ull;
            uint64_t alt = (cur + 1 + rng.range(0, 3)) & 3ull;
            if(alt == cur) continue;
            uint64_t nb = (base & ~(3ull << (2 * pos))) | (alt << (2 * pos));
            bool dup = false;
            for(size_t k = 0; k < umis.size(); k++) if(umis[k] == nb) dup = true;
            if(dup) continue;
            umis.push_back(nb);
            reads.push_back(rng.range(1, 30));
        }

        std::vector<uint32_t> clones;
        soloUmiClones(umis, reads, dedup, umiLen, clones);
        uint32_t want = referenceSurvivors(umis, reads, dedup, umiLen);
        if(clones.size() != want) {
            agree = false;
            std::cerr << "    block of " << umis.size() << ": got " << clones.size()
                      << " survivors, reference says " << want << "\n";
        }
        if(clones.size() < umis.size()) merges++;

        // Reads are conserved: no molecule's support is invented or lost.
        uint64_t inSum = 0, outSum = 0;
        for(size_t k = 0; k < reads.size(); k++) inSum += reads[k];
        for(size_t k = 0; k < clones.size(); k++) outSum += clones[k];
        if(inSum != outSum) {
            agree = false;
            std::cerr << "    reads not conserved: " << inSum << " in, " << outSum << " out\n";
        }
    }
    check(agree, std::string(name) + " agrees with the pre-histogram rule and conserves reads");
    check(merges > 200, std::string(name) + " test blocks actually merge (saw " +
                        std::to_string(merges) + ")");
}

// --------------------------------------------------------------- absorbing

static void testAbsorption() {
    std::cout << "3. a merged UMI's reads join the molecule that absorbed it\n";
    std::vector<uint64_t> umis;
    std::vector<uint32_t> reads, clones;
    umis.push_back(0x000000); reads.push_back(10);
    umis.push_back(0x000001); reads.push_back(1);   // 1MM error of the first
    soloUmiClones(umis, reads, SOLO_UMI_1MM_CR, 12, clones);
    check(clones.size() == 1, "the error is merged away");
    check(clones.size() == 1 && clones[0] == 11, "its read joins the survivor's clone");

    // A chain: c is a 1MM error of b, b of a. Whichever way the chain is
    // resolved, all 12 reads must end up on the surviving molecules.
    umis.clear(); reads.clear();
    umis.push_back(0x000000); reads.push_back(9);
    umis.push_back(0x000001); reads.push_back(2);
    umis.push_back(0x000002); reads.push_back(1);
    soloUmiClones(umis, reads, SOLO_UMI_1MM_CR, 12, clones);
    uint64_t tot = 0;
    for(size_t k = 0; k < clones.size(); k++) tot += clones[k];
    check(tot == 12, "a chain of errors conserves reads");
}

// --------------------------------------------------------------- histogram

static void testHistogram() {
    std::cout << "4. histogram totals and percentile\n";
    SoloCloneHist h;
    for(int i = 0; i < 90; i++) h.add(1);
    for(int i = 0; i < 9; i++)  h.add(5);
    h.add(400);

    check(h.umis() == 100, "molecule count");
    check(h.reads() == 90 + 45 + 400, "read count");
    check(h.maxClone == 400, "max clone size");
    check(h.percentile(0.50) == 1, "median sits in the singleton bin");
    check(h.percentile(0.99) == 400, "P99 reaches the tail");

    // Clone sizes past the last exact bin are counted, not resolved, and must
    // not be lost from the totals.
    SoloCloneHist big;
    big.add(SoloCloneHist::kMaxBin + 7);
    big.add(3);
    check(big.nOver == 1, "oversized clone goes to the overflow counter");
    check(big.umis() == 2, "overflow still counts as a molecule");
    check(big.reads() == SoloCloneHist::kMaxBin + 10, "overflow still counts its reads");
    check(big.maxClone == SoloCloneHist::kMaxBin + 7, "overflow is reflected in the max");

    SoloCloneHist empty;
    check(empty.umis() == 0 && empty.reads() == 0 && empty.percentile(0.99) == 0,
          "an empty histogram reports zeros rather than dividing by zero");
}

int main() {
    testTrivialModes();
    testMatchesReference(SOLO_UMI_1MM_CR,  "1MM_CR");
    testMatchesReference(SOLO_UMI_1MM_ALL, "1MM_All");
    testAbsorption();
    testHistogram();
    std::cout << (failures ? "FAILED " : "passed ") << (checks - failures) << "/" << checks << "\n";
    return failures ? 1 : 0;
}
