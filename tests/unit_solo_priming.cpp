/*
 * Unit tests for the internal-priming window rule.
 *
 *   c++ -std=c++11 -O2 -iquote . tests/unit_solo_priming.cpp -o unit_solo_priming
 *
 * The rule decides whether a molecule gets flagged as an artifact, so both its
 * directions matter: a real polyadenylation site must not trip it (the A's are
 * on the transcript, not in the reference) and a genomic A-stretch must.
 */

#include "solo_priming.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

static int failures = 0, checks = 0;
static void check(bool cond, const std::string& what) {
    checks++;
    if(!cond) { failures++; std::cerr << "  FAIL: " << what << "\n"; }
}

/** ACGT text -> the 2-bit codes BitPairReference hands back. */
static std::vector<char> enc(const std::string& s) {
    std::vector<char> v;
    for(size_t i = 0; i < s.size(); i++) {
        switch(s[i]) {
            case 'A': v.push_back(0); break;
            case 'C': v.push_back(1); break;
            case 'G': v.push_back(2); break;
            case 'T': v.push_back(3); break;
            default:  v.push_back(4); break;
        }
    }
    return v;
}
static bool primed(const std::string& s, int base = 0) {
    std::vector<char> v = enc(s);
    return soloIsPrimingWindow(v.empty() ? NULL : &v[0], v.size(), base);
}

int main() {
    std::cout << "1. the two rules\n";
    check( primed("AAAAAAAAAAAAAAAAAAAA"), "a full window of A is primed");
    check(!primed("CGCGCGCGCGCGCGCGCGCG"), "no A at all is not primed");
    check( primed("AAAAAACGCGCGCGCGCGCG"), "a run of 6 is primed on the run rule alone");
    check(!primed("AAAAACGCGCGCGCGCGCGC"), "a run of 5 with 5 A total is not");
    // 12 A, longest run 3 -- the count rule has to carry this one.
    check( primed("AAACAAACAAACAAACGCGC"), "12 scattered A is primed on count");
    check(!primed("AAACAAACAAACGCGCGCGC"), "9 scattered A is not");

    std::cout << "2. strand\n";
    check( primed("TTTTTTTTTTTTTTTTTTTT", 3), "a window of T is primed when counting T");
    check(!primed("TTTTTTTTTTTTTTTTTTTT", 0), "the same window is not primed counting A");
    check( primed("AAAAAAAAAAAAAAAAAAAA", 0) && !primed("AAAAAAAAAAAAAAAAAAAA", 3),
           "and the mirror holds");

    std::cout << "3. windows truncated by a contig end\n";
    check( primed("AAAAAACGCG"), "6 of 10 trips the run rule");
    check( primed("AACAACAACA"), "6 of 10 scattered meets the scaled count (need 6)");
    check(!primed("AACAACGCGC"), "4 of 10 does not");
    check(!primed(""), "an empty window is not primed");
    check(!primed("A"), "a single base is not primed");

    std::cout << "4. sequences that must not trip it\n";
    // A real poly(A) site: the tail is on the transcript, so the genome past a
    // genuine 3' end looks like ordinary sequence.
    check(!primed("GCTAGCTAGCTAGCTAGCTA"), "ordinary sequence past a real 3' end");
    check(!primed("ACGTACGTACGTACGTACGT"), "even 25% A, evenly spaced, is not");
    check(!primed("NNNNNNNNNNNNNNNNNNNN"), "an N gap is not primed");

    std::cout << (failures ? "FAILED " : "passed ") << (checks - failures) << "/" << checks << "\n";
    return failures ? 1 : 0;
}
