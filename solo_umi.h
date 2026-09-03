/*
 * Copyright 2026, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * HISAT 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SOLO_UMI_H_
#define SOLO_UMI_H_

#include <stdint.h>
#include <vector>

/** UMI collapsing rule. Names match STARsolo so results are comparable. */
enum SoloUmiDedup {
    SOLO_UMI_EXACT = 0,   // distinct UMI sequences
    SOLO_UMI_1MM_CR,      // CellRanger: merge a UMI into a more abundant 1MM neighbour
    SOLO_UMI_1MM_ALL,     // merge all UMIs within one mismatch (connected components)
    SOLO_UMI_NODEDUP      // count reads, not UMIs
};

/** True if two packed UMIs of the given length differ at exactly one base. */
inline bool soloUmiWithin1(uint64_t a, uint64_t b, int len) {
    uint64_t x = a ^ b;
    if(x == 0) return false;
    int diffs = 0;
    for(int i = 0; i < len; i++) {
        if((x >> (2 * i)) & 3ull) {
            if(++diffs > 1) return false;
        }
    }
    return diffs == 1;
}

/**
 * Collapses one (cell, gene) block of UMIs and reports each surviving
 * molecule's clone size -- the number of reads supporting it, including the
 * reads of any UMI merged into it as a sequencing error.
 *
 * `umis` holds the distinct packed UMIs of the block and `reads` their read
 * counts; both are indexed alike. On return `clones` holds one entry per
 * surviving molecule, so `clones.size()` is the deduplicated count for every
 * mode except NODEDUP, where the count is the sum instead. `survivors`, when
 * given, receives the index in `umis` each surviving entry came from, so a
 * caller carrying per-UMI attributes can follow them through the merge.
 *
 * The merge rules are exactly those of the deduplication they serve; the clone
 * sizes are bookkeeping laid over the top and change no decision. Read totals
 * are conserved: sum(clones) == sum(reads).
 */
inline void soloUmiClones(const std::vector<uint64_t>& umis,
                          const std::vector<uint32_t>& reads,
                          SoloUmiDedup dedup,
                          int umiLen,
                          std::vector<uint32_t>& clones,
                          std::vector<uint32_t>* survivors = NULL)
{
    clones.clear();
    if(survivors != NULL) survivors->clear();
    if(dedup == SOLO_UMI_EXACT || dedup == SOLO_UMI_NODEDUP) {
        clones.assign(reads.begin(), reads.end());
        if(survivors != NULL)
            for(size_t a = 0; a < umis.size(); a++) survivors->push_back((uint32_t)a);
        return;
    }

    // Blocks are tiny (a handful of UMIs per gene per cell), so the quadratic
    // neighbour search here is cheaper than building an index.
    const size_t n = umis.size();
    std::vector<char> merged(n, 0);
    std::vector<uint32_t> parent(n);
    for(size_t a = 0; a < n; a++) parent[a] = (uint32_t)a;

    for(size_t a = 0; a < n; a++) {
        if(merged[a]) continue;
        for(size_t b = 0; b < n; b++) {
            if(a == b || merged[b]) continue;
            if(!soloUmiWithin1(umis[a], umis[b], umiLen)) continue;
            if(dedup == SOLO_UMI_1MM_ALL) {
                if(b > a) { merged[b] = 1; parent[b] = (uint32_t)a; }
            } else {
                // 1MM_CR: the less-supported UMI is assumed to be a sequencing
                // error of the more-supported one.
                if(reads[b] < reads[a] || (reads[b] == reads[a] && b > a)) {
                    merged[b] = 1; parent[b] = (uint32_t)a;
                }
            }
        }
    }

    // Give each merged UMI's reads to the survivor that absorbed it. An
    // absorber is unmerged when it takes a child but may itself be absorbed
    // later, so follow the chain to its root. The chain cannot cycle: a UMI
    // that is already merged is skipped by the outer loop and so never becomes
    // an absorber.
    std::vector<uint32_t> clone(reads.begin(), reads.end());
    for(size_t b = 0; b < n; b++) {
        if(!merged[b]) continue;
        size_t r = parent[b], guard = 0;
        while(merged[r] && guard++ < n) r = parent[r];
        clone[r] += reads[b];
    }
    for(size_t a = 0; a < n; a++) {
        if(merged[a]) continue;
        clones.push_back(clone[a]);
        if(survivors != NULL) survivors->push_back((uint32_t)a);
    }
}

/**
 * Distribution of clone sizes over every molecule in a run.
 *
 * A UMI is meant to be attached to one molecule during reverse transcription,
 * which makes clone size a property of sequencing depth alone. Residual
 * UMI-bearing oligos can instead reprime during preamplification PCR and put a
 * second, valid UMI on a molecule that already has one; because the sequence
 * is intact, no amount of error-correcting deduplication can see it. What such
 * "phantom" UMIs do leave behind is a heavier tail in this distribution than
 * depth alone produces, which is why the histogram is emitted rather than the
 * single mean that Sequencing Saturation already reports.
 *
 * See Sugino & Lee, "The Phantom of the PCR" (bioRxiv 2026,
 * doi:10.64898/2026.08.01.742199) for the two-state branching model that
 * turns this shape into a contamination estimate.
 */
struct SoloCloneHist {
    // Exact bins for clone sizes 1..kMaxBin; anything larger is rare and only
    // needs to be counted, not resolved. A run saturated enough to exceed this
    // has other problems.
    static const uint32_t kMaxBin = 65536;

    std::vector<uint64_t> bins;      // bins[k] = molecules with k reads
    uint64_t nOver = 0;              // molecules with more than kMaxBin reads
    uint64_t readsOver = 0;          // and their reads, so totals still balance
    uint32_t maxClone = 0;

    void clear() { bins.clear(); nOver = 0; readsOver = 0; maxClone = 0; }

    void add(uint32_t cloneSize) {
        if(cloneSize == 0) return;
        if(cloneSize > maxClone) maxClone = cloneSize;
        if(cloneSize > kMaxBin) { nOver++; readsOver += cloneSize; return; }
        if(bins.size() <= cloneSize) bins.resize(cloneSize + 1, 0);
        bins[cloneSize]++;
    }

    /** Molecules counted. */
    uint64_t umis() const {
        uint64_t t = nOver;
        for(size_t k = 1; k < bins.size(); k++) t += bins[k];
        return t;
    }

    /** Reads they account for. */
    uint64_t reads() const {
        uint64_t t = readsOver;
        for(size_t k = 1; k < bins.size(); k++) t += bins[k] * (uint64_t)k;
        return t;
    }

    /**
     * Smallest clone size at or below which the given fraction of molecules
     * falls. Reported rather than a mean because the mean is just Sequencing
     * Saturation restated, while the tail is where phantoms show up.
     */
    uint32_t percentile(double p) const {
        const uint64_t tot = umis();
        if(tot == 0) return 0;
        uint64_t want = (uint64_t)(p * (double)tot);
        if(want >= tot) want = tot - 1;
        uint64_t seen = 0;
        for(size_t k = 1; k < bins.size(); k++) {
            seen += bins[k];
            if(seen > want) return (uint32_t)k;
        }
        return maxClone;
    }
};

#endif /* SOLO_UMI_H_ */
