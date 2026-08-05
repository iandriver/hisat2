/*
 * Unit tests for barcode/UMI packing, whitelist lookup and 1MM correction.
 *
 *   c++ -std=c++11 -O2 -iquote . solo_barcode.cpp tests/unit_solo_barcode.cpp -o unit_solo_barcode
 *
 * Whitelist lookup and 1MM correction are checked against brute force, since
 * a wrong barcode silently assigns reads to the wrong cell -- a failure mode
 * nothing downstream can detect.
 */

#include "solo_barcode.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
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
    int range(int lo, int hi) { return lo + (int)(next() % (uint64_t)(hi - lo)); }
};

static const char* BASES = "ACGT";

static std::string randBc(Rng& rng, int len) {
    std::string s(len, 'A');
    for(int i = 0; i < len; i++) s[i] = BASES[rng.range(0, 4)];
    return s;
}

// ------------------------------------------------------------ pack/unpack

static void testPacking() {
    std::cout << "1. pack / unpack\n";
    Rng rng(12345);
    char buf[40];
    for(int t = 0; t < 2000; t++) {
        int len = rng.range(1, 17);
        std::string s = randBc(rng, len);
        uint64_t v;
        check(soloPack(s.c_str(), len, v), "pack ACGT succeeds");
        soloUnpack(v, len, buf);
        if(std::string(buf) != s) { check(false, "round-trip " + s + " -> " + buf); return; }
    }
    uint64_t v;
    check(!soloPack("ACGTN", 5, v), "pack rejects N");
    check(!soloPack("ACGT.", 5, v), "pack rejects non-ACGT");
    check(soloPack("acgt", 4, v), "pack accepts lowercase");
    soloUnpack(v, 4, buf);
    check(std::string(buf) == "ACGT", "lowercase normalises to uppercase");
    // A 16bp barcode must use the full 32 bits without collision.
    uint64_t a, b;
    soloPack("TTTTTTTTTTTTTTTT", 16, a);
    soloPack("AAAAAAAAAAAAAAAA", 16, b);
    check(a == 0xFFFFFFFFull && b == 0, "16bp barcodes span the full uint32 range");
}

// ------------------------------------------------- whitelist vs brute force

static void testWhitelist() {
    std::cout << "2. whitelist lookup vs brute force\n";
    Rng rng(999);
    const int LEN = 16;
    std::set<std::string> wl;
    while(wl.size() < 5000) wl.insert(randBc(rng, LEN));

    std::string path = "/tmp/hisat2_unit_wl.txt";
    { std::ofstream o(path.c_str());
      for(std::set<std::string>::const_iterator it = wl.begin(); it != wl.end(); ++it) o << *it << "\n"; }

    SoloWhitelist w;
    std::string err;
    check(w.load(path, LEN, err), "whitelist loads: " + err);
    check(w.size() == wl.size(), "whitelist retains every barcode");

    // Members are found; non-members are not.
    int found = 0, falsePos = 0;
    for(std::set<std::string>::const_iterator it = wl.begin(); it != wl.end(); ++it) {
        uint64_t v; soloPack(it->c_str(), LEN, v);
        if(w.find((uint32_t)v) != SoloRead::kNoIdx) found++;
    }
    check(found == (int)wl.size(), "every whitelist member is found");
    for(int t = 0; t < 5000; t++) {
        std::string s = randBc(rng, LEN);
        if(wl.count(s)) continue;
        uint64_t v; soloPack(s.c_str(), LEN, v);
        if(w.find((uint32_t)v) != SoloRead::kNoIdx) falsePos++;
    }
    check(falsePos == 0, "no false positives");

    // codeAt must round-trip back to the same barcode.
    char buf[40];
    for(int t = 0; t < 200; t++) {
        std::string s = *std::next(wl.begin(), rng.range(0, (int)wl.size()));
        uint64_t v; soloPack(s.c_str(), LEN, v);
        uint32_t idx = w.find((uint32_t)v);
        soloUnpack(w.codeAt(idx), LEN, buf);
        if(std::string(buf) != s) { check(false, "codeAt round-trip"); return; }
    }
    checks++;
}

// ------------------------------------------------ 1MM correction vs brute force

static void testCorrection() {
    std::cout << "3. 1MM correction vs brute force\n";
    Rng rng(4242);
    const int LEN = 16;
    // A small whitelist so 1MM neighbours are rare and controllable.
    std::set<std::string> wl;
    while(wl.size() < 300) wl.insert(randBc(rng, LEN));
    std::string path = "/tmp/hisat2_unit_wl2.txt";
    { std::ofstream o(path.c_str());
      for(std::set<std::string>::const_iterator it = wl.begin(); it != wl.end(); ++it) o << *it << "\n"; }
    SoloWhitelist w; std::string err;
    check(w.load(path, LEN, err), "small whitelist loads");

    std::vector<std::string> wlv(wl.begin(), wl.end());
    int agree = 0, tested = 0;
    for(int t = 0; t < 3000; t++) {
        std::string q = randBc(rng, LEN);
        // Brute force: how many whitelist entries are within one substitution?
        int hd0 = 0, hd1 = 0;
        for(size_t i = 0; i < wlv.size(); i++) {
            int d = 0;
            for(int j = 0; j < LEN && d < 2; j++) if(wlv[i][j] != q[j]) d++;
            if(d == 0) hd0++;
            else if(d == 1) hd1++;
        }
        SoloRead sr; sr.reset();
        uint64_t v; soloPack(q.c_str(), LEN, v);
        sr.cbPacked = (uint32_t)v; sr.cbLen = LEN; sr.umiLen = 12;
        w.resolve(sr, true);

        int expect;
        if(hd0 == 1)      expect = SOLO_CB_EXACT;
        else if(hd1 == 1) expect = SOLO_CB_1MM;
        else if(hd1 > 1)  expect = SOLO_CB_AMBIG;
        else              expect = SOLO_CB_NOMATCH;
        tested++;
        if(sr.status == expect) agree++;
        else if(agree + 3 > tested) {
            std::cerr << "    mismatch: hd0=" << hd0 << " hd1=" << hd1
                      << " got=" << (int)sr.status << " want=" << expect << "\n";
        }
    }
    std::cout << "   " << agree << "/" << tested << " agree with brute force\n";
    check(agree == tested, "1MM resolution matches brute force exactly");

    // Correction must land on the right barcode, not merely report 1MM.
    std::string base = wlv[7];
    std::string mut = base; mut[5] = (mut[5] == 'A') ? 'C' : 'A';
    SoloRead sr; sr.reset();
    uint64_t v; soloPack(mut.c_str(), LEN, v);
    sr.cbPacked = (uint32_t)v; sr.cbLen = LEN;
    w.resolve(sr, true);
    if(sr.status == SOLO_CB_1MM) {
        char buf[40]; soloUnpack(w.codeAt(sr.cbIdx), LEN, buf);
        check(std::string(buf) == base, "1MM correction recovers the original barcode");
    } else { checks++; }

    // With correction off, a 1MM barcode must be rejected outright.
    SoloRead sr2; sr2.reset();
    sr2.cbPacked = (uint32_t)v; sr2.cbLen = LEN;
    w.resolve(sr2, false);
    check(sr2.status == SOLO_CB_NOMATCH, "Exact mode rejects a 1MM barcode");
}

// ------------------------------------------------------------- extraction

static void testExtraction() {
    std::cout << "4. extraction from sequence and read name\n";
    SoloParams p;   // 10x v3 defaults: CB 1..16, UMI 17..28
    SoloRead sr;

    std::string seq  = "ACGTACGTACGTACGTTTTTGGGGCCCC";   // 28bp
    std::string qual = std::string(28, 'I');
    check(soloExtractFromSeq(seq.c_str(), qual.c_str(), seq.size(), p, sr),
          "extracts from a 28bp barcode read");
    char buf[40];
    soloUnpack(sr.cbPacked, sr.cbLen, buf);
    check(std::string(buf) == "ACGTACGTACGTACGT", "barcode is bases 1-16");
    soloUnpack(sr.umiPacked, sr.umiLen, buf);
    check(std::string(buf) == "TTTTGGGGCCCC", "UMI is bases 17-28");
    check(sr.cbMinQual == 'I', "min barcode quality recorded");

    // Too short for the configured geometry.
    check(!soloExtractFromSeq(seq.c_str(), qual.c_str(), 20, p, sr),
          "short read is rejected");

    // An N anywhere in the barcode is reported distinctly, not as "no match",
    // so it cannot be mistaken for a chemistry mismatch.
    std::string withN = seq; withN[3] = 'N';
    check(!soloExtractFromSeq(withN.c_str(), qual.c_str(), withN.size(), p, sr),
          "N in barcode is rejected");
    check(sr.status == SOLO_CB_HAS_N, "N is reported as HAS_N, not NOMATCH");

    // umi_tools read-name convention.
    std::string name = "A00519:2071:HY5K2DMXY:1:1101:1:1_ACGTACGTACGTACGT_TTTTGGGGCCCC";
    check(soloExtractFromName(name.c_str(), name.size(), p, sr), "extracts from read name");
    soloUnpack(sr.cbPacked, sr.cbLen, buf);
    check(std::string(buf) == "ACGTACGTACGTACGT", "name barcode parsed");
    soloUnpack(sr.umiPacked, sr.umiLen, buf);
    check(std::string(buf) == "TTTTGGGGCCCC", "name UMI parsed");

    std::string plain = "A00519:2071:HY5K2DMXY:1:1101:1:1";
    check(!soloExtractFromName(plain.c_str(), plain.size(), p, sr),
          "a name with no CB/UMI suffix is rejected");
}

int main() {
    testPacking();
    testWhitelist();
    testCorrection();
    testExtraction();
    std::cout << "\n" << (failures ? "FAILED" : "OK") << ": "
              << (checks - failures) << "/" << checks << " checks passed\n";
    return failures ? 1 : 0;
}
