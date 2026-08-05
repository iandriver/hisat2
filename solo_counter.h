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

#ifndef SOLO_COUNTER_H_
#define SOLO_COUNTER_H_

#include <stdint.h>
#include <string>
#include <vector>

#include "gene_model.h"
#include "solo_barcode.h"

#include "ds.h"
#include "read.h"
#include "aligner_result.h"

/** UMI collapsing rule. Names match STARsolo so results are comparable. */
enum SoloUmiDedup {
    SOLO_UMI_EXACT = 0,   // distinct UMI sequences
    SOLO_UMI_1MM_CR,      // CellRanger: merge a UMI into a more abundant 1MM neighbour
    SOLO_UMI_1MM_ALL,     // merge all UMIs within one mismatch (connected components)
    SOLO_UMI_NODEDUP      // count reads, not UMIs
};

/**
 * One (cell, gene, UMI) observation. Exactly 16 bytes with no padding.
 *
 * Gene indices are far below 2^28 for any real annotation (GRCm39 has 33,696),
 * so four bits are free for flags without growing the record. Size matters:
 * a saturated run produces hundreds of millions of these.
 */
struct SoloRec {
    uint32_t cb;             // whitelist index
    uint32_t geneAndFlags;   // gene index in the low 28 bits
    uint64_t umi;            // 2-bit packed

    static const uint32_t kGeneMask = 0x0fffffffu;
    uint32_t gene() const { return geneAndFlags & kGeneMask; }
    void setGene(uint32_t g) { geneAndFlags = (geneAndFlags & ~kGeneMask) | (g & kGeneMask); }
};

/**
 * A read whose barcode had several 1-mismatch whitelist neighbours.
 *
 * CellRanger breaks the tie with a posterior weighted by how often each
 * candidate was seen as an exact match, which is unknowable while streaming.
 * These are parked until finalize(), by which point the exact-match counts
 * exist, and resolved then.
 */
struct SoloAmbigRec {
    uint32_t cbRaw;          // packed barcode as sequenced
    uint32_t geneAndFlags;
    uint64_t umi;
};

class SoloCounter;

/**
 * Per-thread accumulator. Threads never share these, so recording a read
 * takes no lock; everything is merged after the worker threads join.
 */
class SoloCounterThread {
public:
    SoloCounterThread(SoloCounter* parent) : parent_(parent) {}

    /**
     * Records one read.
     *
     * `results` holds every candidate alignment, so this must be called before
     * the aligner picks a primary -- that choice breaks score ties at random,
     * and letting it decide would make counts non-deterministic. When all
     * candidates land on a single gene the read counts as unique-gene even
     * though it is multi-locus, which is what STARsolo does.
     */
    void addRead(const Read& rd, const EList<AlnRes>* results, size_t nresults);

    friend class SoloCounter;

private:
    SoloCounter* parent_;
    std::vector<SoloRec>      recs_;
    std::vector<SoloAmbigRec> ambig_;
    // Per-thread tallies, summed at finalize.
    uint64_t nReads_ = 0, nValidCB_ = 0, nAmbigCB_ = 0, nNoCB_ = 0;
    uint64_t nUnmapped_ = 0, nNoGene_ = 0, nMultiGene_ = 0, nCounted_ = 0;
};

/** Owns configuration, merges the per-thread arenas, and writes the matrices. */
class SoloCounter {
public:
    SoloCounter()
        : gm_(NULL), wl_(NULL), params_(NULL),
          feature_(GENE_FEATURE_EXONIC), strand_(GENE_STRAND_UNSTRANDED),
          dedup_(SOLO_UMI_1MM_CR) {}
    ~SoloCounter();

    void init(const GeneModel* gm,
              const SoloWhitelist* wl,
              const SoloParams* params,
              GeneFeature feature,
              GeneStrand strand,
              SoloUmiDedup dedup,
              const std::string& outDir);

    bool enabled() const { return gm_ != NULL && wl_ != NULL; }

    /**
     * Pre-allocates one accumulator per worker.  Must be called before threads
     * are spawned: allocating lazily from inside the workers would race on the
     * vector itself.
     */
    void reserveThreads(size_t n);
    /** Accumulator for worker i. Read-only once reserveThreads() has run. */
    SoloCounterThread* threadFor(size_t i) const { return threads_[i]; }

    /** Merges, resolves ambiguous barcodes, deduplicates and writes output. */
    bool finalize(std::string& err);

    const GeneModel*     geneModel() const { return gm_; }
    const SoloWhitelist* whitelist() const { return wl_; }
    GeneFeature feature() const { return feature_; }
    GeneStrand  strand()  const { return strand_; }

private:
    bool writeMatrix(const std::vector<SoloRec>& counts,
                     const std::vector<uint32_t>& tally,
                     std::string& err) const;
    bool writeSummary(std::string& err) const;

    const GeneModel*     gm_;
    const SoloWhitelist* wl_;
    const SoloParams*    params_;
    GeneFeature feature_;
    GeneStrand  strand_;
    SoloUmiDedup dedup_;
    std::string outDir_;
    std::vector<SoloCounterThread*> threads_;

    // Totals, filled by finalize().
    uint64_t nReads_ = 0, nValidCB_ = 0, nAmbigCB_ = 0, nAmbigResolved_ = 0;
    uint64_t nNoCB_ = 0, nUnmapped_ = 0, nNoGene_ = 0, nMultiGene_ = 0;
    uint64_t nCounted_ = 0, nUMIs_ = 0, nCells_ = 0, nGenesDetected_ = 0;
};

#endif /* SOLO_COUNTER_H_ */
