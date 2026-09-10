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
#include <utility>
#include <vector>

#include "gene_model.h"
#include "solo_barcode.h"
#include "solo_umi.h"
#include "solo_priming.h"

#include "ds.h"
#include "read.h"
#include "aligner_result.h"

/** Cell-calling rule. Names match STARsolo. */
enum SoloCellFilter {
    SOLO_FILTER_NONE = 0,
    SOLO_FILTER_CELLRANGER22,   // the classic knee: percentile of the top N, over a ratio
    SOLO_FILTER_TOPCELLS,       // fixed number of highest-count barcodes
    SOLO_FILTER_EMPTYDROPS      // deferred to hisat2_solo_filter.py
};

/** What to do with reads whose alignments span several genes. */
enum SoloMultiMapper {
    SOLO_MULTI_UNIQUE = 0,   // discard them (STARsolo's default)
    SOLO_MULTI_UNIFORM,      // split each molecule evenly across its genes
    SOLO_MULTI_EM            // distribute in proportion to per-cell abundance
};

/** How a molecule's splicing state is recorded, for RNA velocity. */
enum SoloVeloClass {
    SOLO_VELO_SPLICED = 0,
    SOLO_VELO_UNSPLICED = 1,
    SOLO_VELO_AMBIGUOUS = 2
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

    // Velocyto class lives in the two bits above the gene index; internal
    // priming takes the next one up. Both ride along for free.
    static const uint32_t kGeneMask   = 0x0fffffffu;
    static const uint32_t kVeloShift  = 28;
    static const uint32_t kVeloMask   = 0x30000000u;
    static const uint32_t kPrimedBit  = 0x40000000u;
    bool primed() const { return (geneAndFlags & kPrimedBit) != 0; }
    void setPrimed(bool v) {
        if(v) geneAndFlags |= kPrimedBit; else geneAndFlags &= ~kPrimedBit;
    }
    uint32_t veloClass() const { return (geneAndFlags & kVeloMask) >> kVeloShift; }
    void setVeloClass(uint32_t c) {
        geneAndFlags = (geneAndFlags & ~kVeloMask) | ((c << kVeloShift) & kVeloMask);
    }

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
/** One allele observation: which variant, and whether the read carried ALT. */
struct SoloAllelicRec {
    uint32_t cb;
    uint32_t variantAndAllele;   // slot in the low 31 bits, ALT flag in bit 31
    uint64_t umi;

    static const uint32_t kAltBit = 0x80000000u;
    uint32_t variant() const { return variantAndAllele & ~kAltBit; }
    bool     isAlt()   const { return (variantAndAllele & kAltBit) != 0; }
};

/**
 * A molecule compatible with several genes.
 *
 * The gene list lives in a separate flat arena rather than inline, so the
 * overwhelmingly common single-gene record stays exactly 16 bytes.
 */
struct SoloMultiRec {
    uint32_t cb;
    uint32_t off;    // start of this record's genes in the arena
    uint32_t n;      // how many genes
    uint64_t umi;
};

struct SoloAmbigRec {
    uint32_t cbRaw;          // packed barcode as sequenced
    uint32_t geneAndFlags;
    uint64_t umi;
};

/**
 * Positions of the index's SNPs, per reference sequence.
 *
 * HISAT2 stores variants in joined-genome coordinates and keeps their real
 * rsIDs (unlike genes, whose identity the index discards), so this converts
 * them once to per-chromosome coordinates and sorts them for range queries.
 * Immutable after build(), so lookups need no locking.
 *
 * One site and allele can have several ALTDB records -- the loader lists every
 * deletion twice, and hisat2-build keeps a .snp line written twice -- and a
 * read's edit names only one of them. build() therefore merges records with
 * the same refid, pos, type, len and seq into one variant, so the read is ALT
 * at that variant rather than ALT at one record and REF at its twin.
 */
class SoloVariantIndex {
public:
    /**
     * pos is 0-based, per-chromosome, and for a deletion its first deleted
     * base. type, len and seq are the ALT's, with seq 0 for a deletion.
     * altIdx indexes ALTDB::alts().
     */
    void add(int32_t refid, int64_t pos, uint32_t type, uint32_t len, uint64_t seq,
             uint32_t altIdx, const std::string& name);
    void build();
    size_t size() const { return nameOff_.size(); }
    bool empty() const { return nameOff_.empty(); }

    /** Variant slots whose position falls in [s, e) on refid. */
    void query(int32_t refid, int64_t s, int64_t e, std::vector<uint32_t>& out) const;
    /** Slot for an ALTDB index, or 0xffffffff if that alt is not a tracked SNP. */
    uint32_t slotForAlt(uint32_t altIdx) const;
    const char* name(uint32_t slot) const { return &nameBlob_[nameOff_[slot]]; }

private:
    struct Entry { int64_t pos; uint32_t slot; };
    std::vector<std::vector<Entry> > byRef_;   // sorted by pos
    std::vector<uint32_t> nameOff_;
    std::vector<char>     nameBlob_;
    std::vector<uint32_t> altToSlot_;          // ALTDB index -> slot
    std::vector<uint64_t> kind_, seq_;         // allele per added record; freed by build()
};

class SoloCounter;

/**
 * Per-thread accumulator. Threads never share these, so recording a read
 * takes no lock; everything is merged after the worker threads join.
 */
class SoloCounterThread {
public:
    SoloCounterThread(SoloCounter* parent);

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

    /** True if the genome just past an alignment's 3' end is A-rich. */
    bool primedAt(const BitPairReference& ref,
                  int32_t refid, int64_t leftmost, int64_t rightmost, bool fw);

    /**
     * Gene-assignment state for one requested feature (Gene, GeneFull, ...).
     *
     * Everything expensive about a read -- alignment, CIGAR, reference blocks --
     * is independent of which feature it is counted under, so a run asking for
     * both features does that work once and only repeats the interval query.
     */
    struct FeatureArena {
        std::vector<SoloRec>      recs;
        std::vector<SoloAmbigRec> ambig;
        std::vector<SoloMultiRec> multi;
        std::vector<uint32_t>     multiGenes;  // arena for the gene lists
        uint64_t nNoGene = 0, nMultiGene = 0, nCounted = 0;
    };

private:
    SoloCounter* parent_;
    std::vector<FeatureArena> feats_;
    std::vector<SoloAllelicRec> allelic_;
    std::vector<SoloRec> velo_;
    // Per-thread tallies, summed at finalize. These count reads, not
    // (read, feature) pairs, so they stay outside FeatureArena.
    uint64_t nReads_ = 0, nValidCB_ = 0, nAmbigCB_ = 0, nNoCB_ = 0;
    uint64_t nUnmapped_ = 0;
    uint64_t nRefObs_ = 0, nAltObs_ = 0;
    uint64_t nPrimed_ = 0;          // reads whose 3' end sits on an A-rich window
    uint64_t nPrimeChecked_ = 0;    // reads that reached a count, so were checked
    // getStretch writes through these; they must not be shared between threads.
    EList<uint32_t> primeBuf_;
    ASSERT_ONLY(SStringExpandable<uint32_t> primeScratch_);
};

/** Owns configuration, merges the per-thread arenas, and writes the matrices. */
class SoloCounter {
public:
    SoloCounter()
        : gm_(NULL), wl_(NULL), params_(NULL), vi_(NULL),
          strand_(GENE_STRAND_UNSTRANDED),
          dedup_(SOLO_UMI_1MM_CR) {}
    ~SoloCounter();

    void init(const GeneModel* gm,
              const SoloWhitelist* wl,
              const SoloParams* params,
              const std::vector<GeneFeature>& features,
              GeneStrand strand,
              SoloUmiDedup dedup,
              const std::string& outDir);

    void setCellFilter(SoloCellFilter f, int expectedCells, double maxPercentile,
                       int maxMinRatio, int topCells) {
        filter_ = f; expectedCells_ = expectedCells;
        maxPercentile_ = maxPercentile; maxMinRatio_ = maxMinRatio; topCells_ = topCells;
    }
    SoloCellFilter cellFilter() const { return filter_; }

    /** Enables allelic output. The index must outlive the counter. */
    void setVariantIndex(const SoloVariantIndex* vi) { vi_ = vi; }
    /** Enables spliced/unspliced/ambiguous output. */
    void setVelocyto(bool v) { velocyto_ = v; }
    /**
     * Enables internal-priming detection. The reference must outlive the
     * counter; without one the check is skipped rather than guessed at.
     */
    void setInternalPriming(bool on, const BitPairReference* ref) {
        primingOn_ = on; ref_ = ref;
    }
    bool internalPriming() const { return primingOn_ && ref_ != NULL; }
    const BitPairReference* reference() const { return ref_; }
    void setMultiMapper(SoloMultiMapper m) { multiMode_ = m; }
    SoloMultiMapper multiMapper() const { return multiMode_; }
    bool velocyto() const { return velocyto_; }
    const SoloVariantIndex* variantIndex() const { return vi_; }

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
    /**
     * The feature currently being written. finalize() walks the requested
     * features one at a time and the output helpers read this, so they need no
     * per-feature plumbing of their own.
     */
    GeneFeature feature() const { return features_[curFeat_]; }
    size_t numFeatures() const { return features_.size(); }
    GeneFeature featureAt(size_t i) const { return features_[i]; }
    GeneStrand  strand()  const { return strand_; }

private:
    /** Counts and writes one requested feature. */
    bool finalizeFeature(size_t fi, std::string& err);
    bool writeMatrix(const std::vector<SoloRec>& counts,
                     const std::vector<uint32_t>& tally,
                     std::string& err) const;
    bool writeSummary(std::string& err) const;
    /** Writes the reads-per-UMI histogram for the current feature. */
    bool writeCloneHist(std::string& err) const;
    /** Writes the same histogram split by gene-expression decile. */
    bool writeCloneHistByExpr(std::string& err) const;
    /** Writes the per-gene internal-priming table. */
    bool writePriming(std::string& err) const;
    /** Writes each called cell's clone-size histogram in the detector's bins. */
    bool writeCloneHistByCell(std::string& err) const;
    bool writeVelocyto(std::vector<SoloRec>& velo, std::string& err);
    bool writeMultiMatrix(const std::vector<SoloRec>& counts,
                          const std::vector<uint32_t>& tally,
                          std::vector<SoloMultiRec>& multi,
                          const std::vector<uint32_t>& multiGenes,
                          std::string& err);
    /** Chooses called barcodes from per-barcode UMI totals. */
    void callCells(const std::vector<SoloRec>& counts,
                   const std::vector<uint32_t>& tally,
                   std::vector<uint32_t>& called) const;
    bool writeFiltered(const std::vector<SoloRec>& counts,
                       const std::vector<uint32_t>& tally,
                       const std::vector<uint32_t>& called,
                       std::string& err) const;
    bool writeAllelic(std::vector<SoloAllelicRec>& allelic, std::string& err);

    const GeneModel*     gm_;
    const SoloWhitelist* wl_;
    const SoloParams*    params_;
    const SoloVariantIndex* vi_;
    bool velocyto_ = false;
    bool primingOn_ = false;
    const BitPairReference* ref_ = NULL;
    SoloMultiMapper multiMode_ = SOLO_MULTI_UNIQUE;
    std::vector<GeneFeature> features_;
    size_t      curFeat_ = 0;
    GeneStrand  strand_;
    SoloUmiDedup dedup_;
    SoloCellFilter filter_ = SOLO_FILTER_CELLRANGER22;
    int    expectedCells_ = 3000;
    double maxPercentile_ = 0.99;
    int    maxMinRatio_   = 10;
    int    topCells_      = 3000;
    std::string outDir_;
    std::vector<SoloCounterThread*> threads_;

    // Totals, filled by finalize().
    uint64_t nReads_ = 0, nValidCB_ = 0, nAmbigCB_ = 0, nAmbigResolved_ = 0;
    uint64_t nNoCB_ = 0, nUnmapped_ = 0, nNoGene_ = 0, nMultiGene_ = 0;
    uint64_t nCounted_ = 0, nUMIs_ = 0, nCells_ = 0, nGenesDetected_ = 0;
    uint64_t nRefUMIs_ = 0, nAltUMIs_ = 0, nVariantsSeen_ = 0;
    uint64_t nRefObs_ = 0, nAltObs_ = 0;
    uint64_t nCalledCells_ = 0, nUMIsInCells_ = 0;
    uint64_t nSpliced_ = 0, nUnspliced_ = 0, nAmbiguous_ = 0;
    uint64_t nMultiUMIs_ = 0;

    // Reads-per-UMI distributions for the feature being written: every
    // barcode, and called cells alone. Ambient barcodes are overwhelmingly
    // single-read molecules and would bury the tail that matters, so the two
    // are kept apart rather than summed.
    SoloCloneHist cloneAll_, cloneCell_;

    // The same distribution for called cells, split by how much of the
    // library a molecule's gene accounts for. A heavy tail confined to the
    // first decile is a handful of very highly expressed genes; one spread
    // across every decile is a property of the chemistry.
    static const size_t kDeciles = 10;
    std::vector<SoloCloneHist> cloneByDecile_;
    std::vector<uint64_t> decileGenes_, decileUmis_;

    // Internal priming, per feature: reads seen by the check, reads called,
    // and the same at the molecule level plus a per-gene breakdown.
    uint64_t nPrimedReads_ = 0, nPrimeChecked_ = 0;
    uint64_t nPrimedUmis_ = 0;
    std::vector<uint64_t> genePrimedUmis_, geneUmis_;

    // One clone-size histogram per called cell, in the detector's bins. The
    // phantom detector is a per-cell method, so a pooled histogram -- which
    // averages over cells of very different depth -- is the wrong input.
    std::vector<uint32_t> cellIdx_;
    std::vector<SoloCloneHist> cellHist_;
};

#endif /* SOLO_COUNTER_H_ */
