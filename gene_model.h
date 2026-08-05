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

#ifndef GENE_MODEL_H_
#define GENE_MODEL_H_

#include <stdint.h>
#include <string>
#include <vector>

/**
 * A gene model loaded from a .ht2gm sidecar (see hisat2_extract_genes.py),
 * with an interval index for "which genes overlap this range".
 *
 * HISAT2's index deliberately has no notion of a gene: annotation stored in
 * the .ht2 files is anonymous (gfm.h names every exon record "exon") and its
 * exon coordinates are fuzzed by +/-10bp as a pseudogene heuristic
 * (splice_site.cpp). Neither is usable for counting, so the gene model is a
 * separate sidecar loaded at align time.
 *
 * Everything here is immutable after load(), so queries need no locking. This
 * is deliberately unlike SpliceSiteDB, which needs a sharded mutex only
 * because it is mutated during alignment.
 *
 * Uses std::vector rather than EList so it can be unit-tested standalone, and
 * to avoid ds.h's push_back_array(), which memcpys non-trivially-copyable
 * types.
 */

/** Which annotation a query is against. */
enum GeneFeature {
    GENE_FEATURE_EXONIC = 0,   // Gene: read must fall within the exon union
    GENE_FEATURE_BODY   = 1    // GeneFull: read must fall within the gene body
};

/** Library strandedness, matching STARsolo's --soloStrand. */
enum GeneStrand {
    GENE_STRAND_UNSTRANDED = 0,
    GENE_STRAND_FORWARD    = 1,   // read strand matches gene strand (10x 3'/5')
    GENE_STRAND_REVERSE    = 2
};

/** One gene. Strings live in a single arena, referenced by offset. */
struct GeneRec {
    uint32_t idOff;        // offset into nameBlob_ of gene_id
    uint32_t nameOff;      // offset into nameBlob_ of gene_name
    int32_t  refid;        // index into the reference-name list
    int8_t   fw;           // 1 if on the + strand
    int64_t  bodyStart;    // 0-based, half-open
    int64_t  bodyEnd;
};

/** An interval carrying the gene it belongs to. */
struct GeneIv {
    int64_t  start;        // 0-based, half-open
    int64_t  end;
    uint32_t gene;
};

/** An annotated junction: 0-based last base of donor exon / first base of acceptor. */
struct GeneJunc {
    int64_t  donor;
    int64_t  acceptor;
    uint32_t gene;
};

/**
 * Per-reference augmented sorted interval array (the cgranges layout).
 *
 * Intervals are sorted by start, with a running maximum of end alongside. A
 * query binary-searches for the first interval starting at or after the query
 * end, then scans backwards while maxEnd still permits an overlap. That gives
 * O(log n + k) with two flat arrays and no pointer chasing.
 *
 * Genes vary enormously in length -- the longest on GRCm39 is ~3 Mb -- and a
 * single long interval raises maxEnd for everything before it, lengthening the
 * backward scan. Intervals at or above kLongIv are therefore hoisted into a
 * short separate list scanned unconditionally, which keeps the main scan tight.
 */
class GeneIvIndex {
public:
    static const int64_t kLongIv = 1000000;

    void build();
    void add(const GeneIv& iv) { ivs_.push_back(iv); }
    /** Appends the genes of all intervals overlapping [s, e). May repeat a gene. */
    void query(int64_t s, int64_t e, std::vector<uint32_t>& out) const;
    size_t size() const { return ivs_.size() + longIvs_.size(); }
    /** Intervals visited by the last query; for tuning, not correctness. */
    mutable uint64_t lastScanned;

    GeneIvIndex() : lastScanned(0) {}

private:
    std::vector<GeneIv>  ivs_;      // sorted by start, short intervals only
    std::vector<int64_t> maxEnd_;   // maxEnd_[i] = max end over ivs_[0..i]
    std::vector<GeneIv>  longIvs_;  // >= kLongIv, scanned every query
};

class GeneModel {
public:
    GeneModel() : loaded_(false) {}

    /**
     * Loads a .ht2gm sidecar, resolving chromosome names against the index.
     *
     * refnames/reflens come from the loaded GFM. Any #ref line naming a
     * sequence the index does not have, or giving a conflicting length, is a
     * hard error: a model built against a different assembly would otherwise
     * produce a plausible-looking but entirely wrong matrix, with nothing
     * downstream able to detect it.
     *
     * Returns false and fills err on failure.
     */
    bool load(const std::string& path,
              const std::vector<std::string>& refnames,
              const std::vector<int64_t>& reflens,
              std::string& err);

    bool loaded() const { return loaded_; }
    size_t numGenes() const { return genes_.size(); }
    const char* geneId(uint32_t g)   const { return &nameBlob_[genes_[g].idOff]; }
    const char* geneName(uint32_t g) const { return &nameBlob_[genes_[g].nameOff]; }
    const GeneRec& gene(uint32_t g)  const { return genes_[g]; }

    /** Genes whose exon union overlaps [s, e) on refid. Appends; may repeat. */
    void overlapExon(int32_t refid, int64_t s, int64_t e,
                     std::vector<uint32_t>& out) const;
    /** Genes whose body overlaps [s, e) on refid. Appends; may repeat. */
    void overlapBody(int32_t refid, int64_t s, int64_t e,
                     std::vector<uint32_t>& out) const;

    /** True if [s, e) lies entirely within one exon of gene g. */
    bool blockInExonsOf(uint32_t g, int64_t s, int64_t e) const;
    /** True if the junction (donor, acceptor) is annotated in gene g. */
    bool junctionInGene(uint32_t g, int64_t donor, int64_t acceptor) const;

    /**
     * Assigns an alignment to gene(s).
     *
     * Follows STARsolo's rule: containment, not mere overlap. Every reference
     * block must lie within the same gene's model, so the result is the
     * intersection of the per-block candidate sets, not their union. A read
     * straddling a gene boundary belongs to neither gene.
     *
     * blocks must be sorted and non-overlapping (see AlnRes::refBlocks).
     * scratch is caller-owned so that steady-state queries do not allocate.
     * Returns the number of genes appended to out.
     */
    size_t assignGenes(int32_t refid,
                       const std::vector<std::pair<int64_t, int64_t> >& blocks,
                       bool readFw,
                       GeneFeature feat,
                       GeneStrand strand,
                       std::vector<uint32_t>& out,
                       std::vector<uint32_t>& scratch) const;

private:
    bool strandOk(uint32_t g, bool readFw, GeneStrand strand) const;

    std::vector<GeneRec>     genes_;
    std::vector<char>        nameBlob_;   // all gene_id/gene_name text
    std::vector<GeneIv>      exons_;      // grouped by gene, CSR-indexed
    std::vector<uint32_t>    exonOff_;    // gene g owns exons_[exonOff_[g]..[g+1])
    std::vector<GeneJunc>    juncs_;      // grouped by gene, CSR-indexed
    std::vector<uint32_t>    juncOff_;
    std::vector<GeneIvIndex> exonIdx_;    // per refid: Gene queries
    std::vector<GeneIvIndex> bodyIdx_;    // per refid: GeneFull queries
    bool loaded_;
};

#endif /* GENE_MODEL_H_ */
