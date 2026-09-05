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

#include "solo_counter.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <sys/stat.h>

#include "read.h"
#include "aligner_result.h"
#include <limits>

namespace {

bool byCbUmiGene(const SoloRec& a, const SoloRec& b) {
    if(a.cb != b.cb)   return a.cb < b.cb;
    if(a.umi != b.umi) return a.umi < b.umi;
    return a.gene() < b.gene();
}

/**
 * One molecule after multi-gene resolution, with the reads behind it.
 *
 * Resolution collapses a (cell, UMI) run to a single record, so without this
 * the read counts would be gone by the time deduplication needs them -- and
 * 1MM_CR, whose whole rule is "merge the less-supported UMI into the better
 * supported one", would silently decide every case on array position instead.
 * There is one of these per molecule rather than per read, so the extra field
 * costs far less than it would inside SoloRec.
 */
struct ResolvedRec {
    SoloRec  rec;
    uint32_t reads;
};

bool byCbGeneUmiR(const ResolvedRec& a, const ResolvedRec& b) {
    if(a.rec.cb != b.rec.cb)         return a.rec.cb < b.rec.cb;
    if(a.rec.gene() != b.rec.gene()) return a.rec.gene() < b.rec.gene();
    return a.rec.umi < b.rec.umi;
}

bool byCbGeneUmi(const SoloRec& a, const SoloRec& b) {
    if(a.cb != b.cb)             return a.cb < b.cb;
    if(a.gene() != b.gene())     return a.gene() < b.gene();
    return a.umi < b.umi;
}

bool makeDir(const std::string& path) {
    struct stat st;
    if(stat(path.c_str(), &st) == 0) return S_ISDIR(st.st_mode);
    return mkdir(path.c_str(), 0755) == 0;
}

bool makeDirs(const std::string& path) {
    // Create each component in turn; the tree is only two or three deep.
    std::string cur;
    for(size_t i = 0; i < path.size(); i++) {
        cur += path[i];
        if(path[i] == '/' && cur.size() > 1) {
            if(!makeDir(cur.substr(0, cur.size() - 1))) return false;
        }
    }
    return makeDir(path);
}

} // namespace

void SoloVariantIndex::add(int32_t refid, int64_t pos, uint32_t altIdx,
                           const std::string& name) {
    if(refid < 0) return;
    uint32_t slot = (uint32_t)nameOff_.size();
    nameOff_.push_back((uint32_t)nameBlob_.size());
    nameBlob_.insert(nameBlob_.end(), name.begin(), name.end());
    nameBlob_.push_back('\0');
    if((size_t)refid >= byRef_.size()) byRef_.resize(refid + 1);
    Entry e; e.pos = pos; e.slot = slot;
    byRef_[refid].push_back(e);
    if(altIdx >= altToSlot_.size()) altToSlot_.resize(altIdx + 1, 0xffffffffu);
    altToSlot_[altIdx] = slot;
}

namespace {
bool entryLess(const std::pair<int64_t, uint32_t>& a, const std::pair<int64_t, uint32_t>& b) {
    return a.first < b.first;
}
}

void SoloVariantIndex::build() {
    for(size_t r = 0; r < byRef_.size(); r++) {
        std::sort(byRef_[r].begin(), byRef_[r].end(),
                  [](const Entry& a, const Entry& b) { return a.pos < b.pos; });
    }
}

uint32_t SoloVariantIndex::slotForAlt(uint32_t altIdx) const {
    if(altIdx >= altToSlot_.size()) return 0xffffffffu;
    return altToSlot_[altIdx];
}

void SoloVariantIndex::query(int32_t refid, int64_t s, int64_t e,
                             std::vector<uint32_t>& out) const {
    if(refid < 0 || (size_t)refid >= byRef_.size()) return;
    const std::vector<Entry>& v = byRef_[refid];
    if(v.empty()) return;
    // Variants are points, so a plain lower_bound on position suffices.
    size_t lo = 0, hi = v.size();
    while(lo < hi) { size_t mid = (lo + hi) / 2; if(v[mid].pos < s) lo = mid + 1; else hi = mid; }
    for(size_t i = lo; i < v.size() && v[i].pos < e; i++) out.push_back(v[i].slot);
}

SoloCounter::~SoloCounter() {
    for(size_t i = 0; i < threads_.size(); i++) delete threads_[i];
}

void SoloCounter::init(const GeneModel* gm,
                       const SoloWhitelist* wl,
                       const SoloParams* params,
                       const std::vector<GeneFeature>& features,
                       GeneStrand strand,
                       SoloUmiDedup dedup,
                       const std::string& outDir) {
    gm_ = gm; wl_ = wl; params_ = params;
    features_ = features; curFeat_ = 0; strand_ = strand; dedup_ = dedup;
    outDir_ = outDir;
}

SoloCounterThread::SoloCounterThread(SoloCounter* parent) : parent_(parent) {
    // Sized here rather than on first use: the workers run concurrently, so a
    // lazy resize would be a data race on the arena vector.
    feats_.resize(parent->numFeatures());
}

void SoloCounter::reserveThreads(size_t n) {
    for(size_t i = threads_.size(); i < n; i++) threads_.push_back(new SoloCounterThread(this));
}

/**
 * Reads the genome just past an alignment's 3' end and asks whether it looks
 * like an oligo(dT) mispriming site.
 *
 * "Past the 3' end" is in the read's own orientation: forward alignments look
 * rightward and count A, reverse ones look leftward and count T, because the
 * transcript's A's are the reference's T's on that strand. A window running off
 * the end of a contig is judged on what is there.
 */
bool SoloCounterThread::primedAt(const BitPairReference& ref,
                                 int32_t refid, int64_t leftmost, int64_t rightmost,
                                 bool fw) {
    if(refid < 0) return false;
    const int64_t reflen = (int64_t)ref.approxLen((TIndexOffU)refid);
    int64_t from, to;
    if(fw) { from = rightmost;          to = from + kPrimeWindow; }
    else   { to   = leftmost;           from = to - kPrimeWindow; }
    if(from < 0) from = 0;
    if(to > reflen) to = reflen;
    const int64_t n = to - from;
    if(n <= 0) return false;

    primeBuf_.resize((size_t)((n + 16 + 3) / 4) + 1);
    int off = ref.getStretch(
        primeBuf_.ptr(),
        (size_t)refid,
        (size_t)from,
        (size_t)n
        ASSERT_ONLY(, primeScratch_));
    const char* w = (const char*)primeBuf_.ptr() + off;
    return soloIsPrimingWindow(w, (size_t)n, fw ? 0 : 3);
}

void SoloCounterThread::addRead(const Read& rd, const EList<AlnRes>* results,
                                size_t nresults) {
    const SoloCounter* p = parent_;
    nReads_++;

    const SoloRead& sr = rd.solo;
    if(!sr.hasBarcode()) { nNoCB_++; return; }
    if(sr.status == SOLO_CB_AMBIG) {
        nAmbigCB_++;
    } else if(!sr.corrected()) {
        nNoCB_++;
        return;
    } else {
        nValidCB_++;
    }

    if(results == NULL || nresults == 0) { nUnmapped_++; return; }

    // Union the gene assignments over the candidate alignments the aligner
    // would actually report. If they all agree on one gene the read is
    // unique-gene even when it is multi-locus, which is what STARsolo counts.
    //
    // "Would actually report" is load-bearing. finishRead hands us every
    // candidate, including ones selectByScore then discards, and ranking by
    // the raw score alone does not reproduce its choice: a read spanning a
    // junction aligns to its parent gene spliced AND to a processed retrogene
    // contiguously, both at raw score 0. Ranking by hisat2_score() -- which
    // packs the raw score above the repeat / known-transcript / splice-site /
    // trim preferences (aligner_result.h calculate_hisat2_score) -- is what
    // makes the spliced parent win, exactly as selectByScore does.
    //
    // Without this filter every ribosomal-protein read with a retrogene copy
    // looked multi-gene and was discarded: Rps27 counted 258 UMIs against
    // rustar's 9,293 on 10M mouse reads. We keep ALL candidates tied at the
    // top rather than picking one, so selectByScore's random tiebreak still
    // never influences a count.
    static thread_local std::vector<uint32_t> scratch, blockGenes;
    static thread_local std::vector<std::vector<uint32_t> > featGenes;
    static thread_local std::vector<uint32_t> bodyGenes, veloGenes, veloClasses;
    static thread_local std::vector<std::pair<int64_t, int64_t> > blocks;
    static thread_local StackedAln staln;
    const size_t nfeat = feats_.size();
    if(featGenes.size() < nfeat) featGenes.resize(nfeat);
    for(size_t fi = 0; fi < nfeat; fi++) featGenes[fi].clear();
    veloGenes.clear(); veloClasses.clear();
    // -1 until the genome has been consulted; see the deferred check below.
    int  primedRead = -1;
    int32_t primeRef = -1;
    bool primeFw = false;
    int64_t primeL = 0, primeR = 0;

    // Rank exactly as selectByScore does, then keep only the top group.
    TAlScore bestScore = std::numeric_limits<TAlScore>::min();
    for(size_t i = 0; i < nresults; i++) {
        TAlScore h = (*results)[i].score().hisat2_score();
        if(h > bestScore) bestScore = h;
    }

    for(size_t i = 0; i < nresults; i++) {
        const AlnRes& rs = (*results)[i];
        if(rs.score().hisat2_score() != bestScore) continue;
        // The CIGAR is the authoritative description of which reference bases
        // the alignment covers, and reusing it keeps the counting path exactly
        // consistent with the emitted GX:Z tag.
        staln.reset();
        const_cast<AlnRes&>(rs).initStacked(rd, staln);
        staln.leftAlign(false);
        staln.buildCigar(false);
        if(!staln.cigarBuilt()) continue;

        blocks.clear();
        int64_t pos = (int64_t)rs.refoff();
        int64_t blockStart = pos;
        bool haveBlock = false;
        const EList<char>&   ops  = staln.cigarOps();
        const EList<size_t>& runs = staln.cigarRuns();
        for(size_t k = 0; k < ops.size(); k++) {
            char op = ops[k];
            int64_t run = (int64_t)runs[k];
            if(op == 'M' || op == '=' || op == 'X' || op == 'D') {
                if(!haveBlock) { blockStart = pos; haveBlock = true; }
                pos += run;
            } else if(op == 'N') {
                if(haveBlock) { blocks.push_back(std::make_pair(blockStart, pos)); haveBlock = false; }
                pos += run;
            }
        }
        if(haveBlock) blocks.push_back(std::make_pair(blockStart, pos));
        if(blocks.empty()) continue;

        // Remember where the first alignment the aligner would report ends --
        // the first to survive the top-score filter above, which is not
        // necessarily candidate 0. The check itself is deferred: reading the
        // genome there costs a random access into three gigabytes, and better
        // than half of all reads never reach a count.
        if(primeRef < 0 && p->internalPriming()) {
            primeRef = (int32_t)rs.refid();
            primeFw  = rs.fw();
            primeL   = blocks.front().first;
            primeR   = blocks.back().second;
        }

        // One interval query per requested feature. The blocks above are what
        // cost anything to produce, and they are shared.
        for(size_t fi = 0; fi < nfeat; fi++) {
            blockGenes.clear();
            p->geneModel()->assignGenes((int32_t)rs.refid(), blocks, rs.fw(),
                                        p->featureAt(fi), p->strand(), blockGenes, scratch);
            for(size_t g = 0; g < blockGenes.size(); g++) featGenes[fi].push_back(blockGenes[g]);
        }

        // Velocyto works off the gene body, not the exon union: an intronic
        // read is exactly the signal of interest, and the exonic assignment
        // would discard it.
        if(p->velocyto()) {
            bodyGenes.clear();
            p->geneModel()->assignGenes((int32_t)rs.refid(), blocks, rs.fw(),
                                        GENE_FEATURE_BODY, p->strand(), bodyGenes, scratch);
            const GeneModel* gm = p->geneModel();
            for(size_t bg = 0; bg < bodyGenes.size(); bg++) {
                uint32_t g = bodyGenes[bg];
                bool intronic = false, splicedJn = false;
                for(size_t b = 0; b < blocks.size(); b++) {
                    // A block not contained in a single exon either lies in an
                    // intron or straddles a boundary; both are unspliced evidence.
                    if(!gm->blockInExonsOf(g, blocks[b].first, blocks[b].second)) intronic = true;
                    if(b + 1 < blocks.size() &&
                       gm->junctionInGene(g, blocks[b].second - 1, blocks[b + 1].first)) {
                        splicedJn = true;
                    }
                }
                uint32_t cls = intronic ? (splicedJn ? SOLO_VELO_AMBIGUOUS : SOLO_VELO_UNSPLICED)
                                        : SOLO_VELO_SPLICED;
                veloGenes.push_back(g);
                veloClasses.push_back(cls);
            }
        }
    }

    // Allele observations.  A read that took an alternate path through the
    // graph carries an edit tagged with that variant's ALTDB index; any tracked
    // variant inside the alignment without such an edit was seen as reference.
    // Restricted to uniquely-aligned reads: for a multi-locus read the alleles
    // it supports are ambiguous.
    const SoloVariantIndex* vi = p->variantIndex();
    if(vi != NULL && !vi->empty() && nresults == 1 && sr.corrected()) {
        static thread_local std::vector<uint32_t> altSeen, spanned;
        const AlnRes& rs = (*results)[0];
        altSeen.clear(); spanned.clear();
        const EList<Edit>& eds = rs.ned();
        for(size_t k = 0; k < eds.size(); k++) {
            uint32_t slot = vi->slotForAlt((uint32_t)eds[k].snpID);
            if(slot != 0xffffffffu) altSeen.push_back(slot);
        }
        std::sort(altSeen.begin(), altSeen.end());
        for(size_t b = 0; b < blocks.size(); b++) {
            vi->query((int32_t)rs.refid(), blocks[b].first, blocks[b].second, spanned);
        }
        std::sort(spanned.begin(), spanned.end());
        spanned.erase(std::unique(spanned.begin(), spanned.end()), spanned.end());
        for(size_t k = 0; k < spanned.size(); k++) {
            bool isAlt = std::binary_search(altSeen.begin(), altSeen.end(), spanned[k]);
            SoloAllelicRec ar;
            ar.cb = sr.cbIdx;
            ar.variantAndAllele = spanned[k] | (isAlt ? SoloAllelicRec::kAltBit : 0u);
            ar.umi = sr.umiPacked;
            allelic_.push_back(ar);
            if(isAlt) nAltObs_++; else nRefObs_++;
        }
    }


    if(p->velocyto() && !veloGenes.empty() && sr.corrected()) {
        // Collapse across alignments: a molecule showing both intronic and
        // spliced evidence is ambiguous, matching how the per-UMI merge below
        // resolves multiple reads.
        uint32_t g0 = veloGenes[0];
        bool single = true, anyUnspliced = false, anySpliced = false;
        for(size_t k = 0; k < veloGenes.size(); k++) {
            if(veloGenes[k] != g0) { single = false; break; }
            if(veloClasses[k] == SOLO_VELO_UNSPLICED) anyUnspliced = true;
            else if(veloClasses[k] == SOLO_VELO_SPLICED) anySpliced = true;
            else { anyUnspliced = true; anySpliced = true; }
        }
        if(single) {
            SoloRec vr;
            vr.cb = sr.cbIdx;
            vr.geneAndFlags = g0;
            vr.umi = sr.umiPacked;
            vr.setVeloClass(anyUnspliced && anySpliced ? SOLO_VELO_AMBIGUOUS
                            : (anyUnspliced ? SOLO_VELO_UNSPLICED : SOLO_VELO_SPLICED));
            velo_.push_back(vr);
        }
    }

    for(size_t fi = 0; fi < nfeat; fi++) {
        std::vector<uint32_t>& genes = featGenes[fi];
        FeatureArena& fa = feats_[fi];
        std::sort(genes.begin(), genes.end());
        genes.erase(std::unique(genes.begin(), genes.end()), genes.end());

        if(genes.empty())      { fa.nNoGene++;    continue; }
        if(genes.size() > 1) {
            fa.nMultiGene++;
            // Keep the gene set when multimapper resolution is on: EM needs only
            // this, not the alignments, which is what lets resolution happen at
            // finalize without holding anything from the streaming pass.
            if(p->multiMapper() != SOLO_MULTI_UNIQUE && sr.corrected()) {
                SoloMultiRec mr;
                mr.cb  = sr.cbIdx;
                mr.off = (uint32_t)fa.multiGenes.size();
                mr.n   = (uint32_t)genes.size();
                mr.umi = sr.umiPacked;
                for(size_t k = 0; k < genes.size(); k++) fa.multiGenes.push_back(genes[k]);
                fa.multi.push_back(mr);
            }
            continue;
        }

        fa.nCounted++;
        if(primedRead < 0 && p->internalPriming() && primeRef >= 0) {
            nPrimeChecked_++;
            primedRead = primedAt(*p->reference(), primeRef, primeL, primeR, primeFw) ? 1 : 0;
            if(primedRead) nPrimed_++;
        }
        if(sr.status == SOLO_CB_AMBIG) {
            SoloAmbigRec ar;
            ar.cbRaw = sr.cbPacked;
            ar.geneAndFlags = genes[0] | (primedRead > 0 ? SoloRec::kPrimedBit : 0u);
            ar.umi = sr.umiPacked;
            fa.ambig.push_back(ar);
        } else {
            SoloRec r;
            r.cb = sr.cbIdx;
            r.geneAndFlags = genes[0];
            r.umi = sr.umiPacked;
            r.setPrimed(primedRead > 0);
            fa.recs.push_back(r);
        }
    }
}

bool SoloCounter::finalize(std::string& err) {
    if(!enabled()) { err = "solo counter not initialised"; return false; }

    // Streams that do not depend on the feature are merged once. Velocyto
    // already queries the gene body explicitly and allele observations are
    // positional, so neither is repeated per feature.
    std::vector<SoloAllelicRec> allelic;
    std::vector<SoloRec> velo;
    for(size_t i = 0; i < threads_.size(); i++) {
        SoloCounterThread* t = threads_[i];
        allelic.insert(allelic.end(), t->allelic_.begin(), t->allelic_.end());
        velo.insert(velo.end(), t->velo_.begin(), t->velo_.end());
        nRefObs_ += t->nRefObs_; nAltObs_ += t->nAltObs_;
        nReads_ += t->nReads_; nValidCB_ += t->nValidCB_; nAmbigCB_ += t->nAmbigCB_;
        nPrimedReads_ += t->nPrimed_; nPrimeChecked_ += t->nPrimeChecked_;
        nNoCB_ += t->nNoCB_;   nUnmapped_ += t->nUnmapped_;
        std::vector<SoloAllelicRec>().swap(t->allelic_);
        std::vector<SoloRec>().swap(t->velo_);
    }

    // Before the per-feature loop, not after: these populate the velocyto and
    // allelic tallies that each feature's Summary.csv reports, and
    // finalizeFeature() writes that summary. Running them afterwards leaves
    // every one of those lines reading zero while the matrices themselves are
    // correct -- which is exactly what it did until this was caught.
    if(!velo.empty() && !writeVelocyto(velo, err)) return false;
    if(!allelic.empty() && !writeAllelic(allelic, err)) return false;

    for(size_t fi = 0; fi < features_.size(); fi++) {
        curFeat_ = fi;
        if(!finalizeFeature(fi, err)) return false;
    }
    return true;
}

/**
 * Counts and writes one feature. Called once per --gene-feature entry, with
 * curFeat_ set so the output helpers pick up the right name and directory.
 */
bool SoloCounter::finalizeFeature(size_t fi, std::string& err) {
    // Per-feature totals. A run counting both Gene and GeneFull writes two
    // Summary.csv files, and each must describe its own feature.
    nNoGene_ = nMultiGene_ = nCounted_ = 0;
    nUMIs_ = nCells_ = nGenesDetected_ = 0;
    nCalledCells_ = nUMIsInCells_ = nMultiUMIs_ = nAmbigResolved_ = 0;
    cloneAll_.clear(); cloneCell_.clear();
    cloneByDecile_.assign(kDeciles, SoloCloneHist());
    cellIdx_.clear(); cellHist_.clear();
    nPrimedUmis_ = 0;
    genePrimedUmis_.assign(gm_ != NULL ? gm_->numGenes() : 0, 0);
    geneUmis_.assign(gm_ != NULL ? gm_->numGenes() : 0, 0);
    decileGenes_.assign(kDeciles, 0);
    decileUmis_.assign(kDeciles, 0);

    std::vector<SoloRec> all;
    std::vector<SoloAmbigRec> ambig;
    std::vector<SoloMultiRec> multi;
    std::vector<uint32_t> multiGenes;
    size_t total = 0;
    for(size_t i = 0; i < threads_.size(); i++) total += threads_[i]->feats_[fi].recs.size();
    all.reserve(total);
    for(size_t i = 0; i < threads_.size(); i++) {
        SoloCounterThread::FeatureArena& fa = threads_[i]->feats_[fi];
        all.insert(all.end(), fa.recs.begin(), fa.recs.end());
        ambig.insert(ambig.end(), fa.ambig.begin(), fa.ambig.end());
        // Gene-list offsets are thread-local, so rebase them onto the merged arena.
        {
            uint32_t base = (uint32_t)multiGenes.size();
            for(size_t k = 0; k < fa.multi.size(); k++) {
                SoloMultiRec mr = fa.multi[k];
                mr.off += base;
                multi.push_back(mr);
            }
            multiGenes.insert(multiGenes.end(), fa.multiGenes.begin(), fa.multiGenes.end());
        }
        nNoGene_ += fa.nNoGene; nMultiGene_ += fa.nMultiGene; nCounted_ += fa.nCounted;
        // Release each arena as it is consumed, so peak memory is roughly one
        // copy of the records rather than two.
        std::vector<SoloRec>().swap(fa.recs);
        std::vector<SoloAmbigRec>().swap(fa.ambig);
        std::vector<SoloMultiRec>().swap(fa.multi);
        std::vector<uint32_t>().swap(fa.multiGenes);
    }

    // Ambiguous barcodes: now that the pass is complete, the exact-match count
    // of every whitelist barcode is known, so a 1MM barcode with several
    // candidates can go to whichever candidate was actually observed most.
    if(!ambig.empty()) {
        std::map<uint32_t, uint32_t> obs;
        for(size_t i = 0; i < all.size(); i++) obs[all[i].cb]++;
        int cbLen = wl_->cbLen();
        for(size_t i = 0; i < ambig.size(); i++) {
            uint32_t best = SoloRead::kNoIdx, bestCount = 0;
            bool tie = false;
            for(int pos = 0; pos < cbLen; pos++) {
                int shift = 2 * (cbLen - 1 - pos);
                uint32_t cur = (ambig[i].cbRaw >> shift) & 3u;
                for(uint32_t alt = 0; alt < 4; alt++) {
                    if(alt == cur) continue;
                    uint32_t cand = (ambig[i].cbRaw & ~(3u << shift)) | (alt << shift);
                    uint32_t idx = wl_->find(cand);
                    if(idx == SoloRead::kNoIdx) continue;
                    std::map<uint32_t, uint32_t>::const_iterator it = obs.find(idx);
                    uint32_t c = (it == obs.end()) ? 0 : it->second;
                    if(best == SoloRead::kNoIdx || c > bestCount) { best = idx; bestCount = c; tie = false; }
                    else if(c == bestCount) tie = true;
                }
            }
            if(best != SoloRead::kNoIdx && !tie) {
                SoloRec r;
                r.cb = best;
                r.geneAndFlags = ambig[i].geneAndFlags;
                r.umi = ambig[i].umi;
                all.push_back(r);
                nAmbigResolved_++;
            }
        }
    }

    if(all.empty()) { err = "no reads were assigned to a cell and gene"; return false; }

    // Multi-gene UMIs: the same UMI in the same cell hitting two genes is a
    // collision. CellRanger gives it to whichever gene has more reads and
    // discards it on a tie.
    std::sort(all.begin(), all.end(), byCbUmiGene);
    std::vector<ResolvedRec> resolved;
    resolved.reserve(all.size());
    size_t i = 0;
    while(i < all.size()) {
        size_t j = i;
        while(j < all.size() && all[j].cb == all[i].cb && all[j].umi == all[i].umi) j++;
        // Within this (cell, UMI) run, find the most-supported gene.
        uint32_t bestGene = all[i].gene(), bestN = 0;
        bool tie = false;
        size_t k = i;
        while(k < j) {
            size_t m = k;
            while(m < j && all[m].gene() == all[k].gene()) m++;
            uint32_t n = (uint32_t)(m - k);
            if(n > bestN)      { bestN = n; bestGene = all[k].gene(); tie = false; }
            else if(n == bestN && all[k].gene() != bestGene) tie = true;
            k = m;
        }
        if(!tie) {
            ResolvedRec rr;
            rr.rec = all[i];
            rr.rec.setGene(bestGene);
            // Every read in the run came off the same molecule, including any
            // that favoured a losing gene, so they all count towards its depth.
            rr.reads = (uint32_t)(j - i);
            resolved.push_back(rr);
        } else {
            nMultiGene_++;
        }
        i = j;
    }
    std::vector<SoloRec>().swap(all);

    // UMI collapsing within each (cell, gene).
    std::sort(resolved.begin(), resolved.end(), byCbGeneUmiR);
    std::vector<SoloRec> counts;    // one entry per (cell, gene)
    std::vector<uint32_t> tally;    // its UMI count
    std::vector<uint32_t> clones;   // reads per surviving molecule, per block
    std::vector<uint32_t> survivors;
    std::vector<uint32_t> cloneArena, cloneOff, cloneN;   // the same, kept
    int umiLen = params_ != NULL ? params_->umiLen : 12;

    i = 0;
    while(i < resolved.size()) {
        size_t j = i;
        while(j < resolved.size() && resolved[j].rec.cb == resolved[i].rec.cb &&
              resolved[j].rec.gene() == resolved[i].rec.gene()) j++;

        // Distinct UMIs in this block, with read support.
        std::vector<uint64_t> umis;
        std::vector<uint32_t> reads;
        static thread_local std::vector<char> primedFlag;
        primedFlag.clear();
        uint64_t blockReads = 0;
        for(size_t k = i; k < j; ) {
            size_t m = k;
            uint32_t nr = 0;
            bool pr = false;
            while(m < j && resolved[m].rec.umi == resolved[k].rec.umi) {
                nr += resolved[m].reads; pr = pr || resolved[m].rec.primed(); m++;
            }
            umis.push_back(resolved[k].rec.umi);
            reads.push_back(nr);
            primedFlag.push_back(pr ? 1 : 0);
            blockReads += nr;
            k = m;
        }

        soloUmiClones(umis, reads, dedup_, umiLen, clones, &survivors);
        if(primingOn_) {
            const uint32_t g = resolved[i].rec.gene();
            for(size_t x = 0; x < survivors.size(); x++) {
                if(!primedFlag[survivors[x]]) continue;
                nPrimedUmis_++;
                if(g < genePrimedUmis_.size()) genePrimedUmis_[g]++;
            }
            if(g < geneUmis_.size()) geneUmis_[g] += clones.size();
        }

        // Every mode but NODEDUP counts molecules; NODEDUP counts their reads.
        uint32_t n = 0;
        if(dedup_ == SOLO_UMI_NODEDUP) {
            n = (uint32_t)blockReads;
        } else {
            n = (uint32_t)clones.size();
        }

        if(n > 0) {
            counts.push_back(resolved[i].rec);
            tally.push_back(n);
            nUMIs_ += n;
            // Clone sizes are binned once cell calling has run, so park them
            // in a flat arena keyed by position in counts. One uint32 per
            // molecule, which is an order of magnitude smaller than the
            // per-read records already held.
            cloneOff.push_back((uint32_t)cloneArena.size());
            cloneN.push_back((uint32_t)clones.size());
            cloneArena.insert(cloneArena.end(), clones.begin(), clones.end());
        }
        i = j;
    }

    // Cells and genes with any signal, for the summary.
    {
        std::vector<uint32_t> cbs, gs;
        cbs.reserve(counts.size()); gs.reserve(counts.size());
        for(size_t k = 0; k < counts.size(); k++) { cbs.push_back(counts[k].cb); gs.push_back(counts[k].gene()); }
        std::sort(cbs.begin(), cbs.end());
        cbs.erase(std::unique(cbs.begin(), cbs.end()), cbs.end());
        nCells_ = cbs.size();
        std::sort(gs.begin(), gs.end());
        gs.erase(std::unique(gs.begin(), gs.end()), gs.end());
        nGenesDetected_ = gs.size();
    }

    // A very low assignment rate almost always means --gene-strand is wrong;
    // the failure is otherwise silent and produces a plausible but tiny matrix.
    if(nReads_ > 0 && strand_ != GENE_STRAND_UNSTRANDED &&
       (double)nCounted_ / (double)nReads_ < 0.10) {
        fprintf(stderr,
                "Warning: only %.1f%% of reads were assigned to a gene. This "
                "usually means --gene-strand is inverted for this library; try "
                "the opposite setting or Unstranded.\n",
                100.0 * (double)nCounted_ / (double)nReads_);
    }

    if(!writeMatrix(counts, tally, err)) return false;

    // Which barcodes are cells is needed both by the filtered matrix and by
    // the clone-size histogram below, so it is decided once here.
    std::vector<char> isCell(wl_->size(), 0);
    if(filter_ != SOLO_FILTER_NONE) {
        std::vector<uint32_t> called;
        callCells(counts, tally, called);
        nCalledCells_ = called.size();
        for(size_t i = 0; i < called.size(); i++) isCell[called[i]] = 1;
        for(size_t k = 0; k < counts.size(); k++)
            if(isCell[counts[k].cb]) nUMIsInCells_ += tally[k];
        if(!writeFiltered(counts, tally, called, err)) return false;

        // EmptyDrops_CR needs an ambient-profile estimate, a Monte-Carlo
        // simulation and an FDR correction, which live in a bundled Python
        // script. The knee result above is already on disk, so if the script
        // cannot run the user still has a usable filtered matrix.
        if(filter_ == SOLO_FILTER_EMPTYDROPS) {
            const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
            std::string raw = outDir_ + "/" + featName + "/raw";
            std::string script = "hisat2_solo_filter.py";
            std::string cmd = "python3 " + script + " --raw '" + raw +
                              "' --expect-cells " + std::to_string(expectedCells_);
            if(system((cmd + " 2>/dev/null >/dev/null").c_str()) != 0) {
                fprintf(stderr,
                        "Note: could not run EmptyDrops_CR automatically. The knee-filtered "
                        "matrix has been written; to refine it run:\n  %s\n", cmd.c_str());
            } else {
                fprintf(stderr, "EmptyDrops_CR refinement complete.\n");
            }
        }
    }

    // Rank genes by how many molecules they hold in called cells, then cut
    // them into deciles of molecule mass -- not of gene count, which would put
    // nearly every molecule in one bin. The first decile is therefore the
    // handful of genes carrying a tenth of the library.
    const size_t nGenes = gm_ != NULL ? gm_->numGenes() : 0;
    std::vector<uint8_t> geneDecile(nGenes, (uint8_t)(kDeciles - 1));
    {
        std::vector<uint64_t> geneUmis(nGenes, 0);
        uint64_t totalUmis = 0;
        for(size_t k = 0; k < counts.size(); k++) {
            if(!isCell[counts[k].cb] && filter_ != SOLO_FILTER_NONE) continue;
            geneUmis[counts[k].gene()] += tally[k];
            totalUmis += tally[k];
        }
        std::vector<uint32_t> order;
        order.reserve(nGenes);
        for(size_t g = 0; g < nGenes; g++) if(geneUmis[g]) order.push_back((uint32_t)g);
        std::sort(order.begin(), order.end(),
                  [&geneUmis](uint32_t a, uint32_t b) {
                      if(geneUmis[a] != geneUmis[b]) return geneUmis[a] > geneUmis[b];
                      return a < b;   // stable across runs
                  });
        uint64_t seen = 0;
        for(size_t idx = 0; idx < order.size(); idx++) {
            size_t d = totalUmis ? (size_t)((double)seen * kDeciles / (double)totalUmis) : 0;
            if(d >= kDeciles) d = kDeciles - 1;
            geneDecile[order[idx]] = (uint8_t)d;
            decileGenes_[d]++;
            decileUmis_[d] += geneUmis[order[idx]];
            seen += geneUmis[order[idx]];
        }
    }

    // Bin the parked clone sizes now that cells are known, then release the
    // arena before the multimapper pass allocates.
    // Per-cell histograms are indexed through a sparse map so only called
    // cells occupy anything; a whitelist-sized vector of histograms would be
    // millions of empty objects.
    std::map<uint32_t, size_t> cellSlot;
    for(size_t k = 0; k < cloneOff.size(); k++) {
        const bool cell = isCell[counts[k].cb] != 0;
        const uint32_t g = counts[k].gene();
        const uint8_t d = g < geneDecile.size() ? geneDecile[g] : (uint8_t)(kDeciles - 1);
        SoloCloneHist* ch = NULL;
        if(cell) {
            std::map<uint32_t, size_t>::iterator it = cellSlot.find(counts[k].cb);
            if(it == cellSlot.end()) {
                it = cellSlot.insert(std::make_pair(counts[k].cb, cellHist_.size())).first;
                cellHist_.push_back(SoloCloneHist());
                cellIdx_.push_back(counts[k].cb);
            }
            ch = &cellHist_[it->second];
        }
        for(uint32_t x = 0; x < cloneN[k]; x++) {
            const uint32_t c = cloneArena[cloneOff[k] + x];
            cloneAll_.add(c);
            if(cell) { cloneCell_.add(c); cloneByDecile_[d].add(c); ch->add(c); }
        }
    }
    std::vector<uint32_t>().swap(cloneArena);
    std::vector<uint32_t>().swap(cloneOff);
    std::vector<uint32_t>().swap(cloneN);
    if(!writeCloneHist(err)) return false;
    if(!writeCloneHistByExpr(err)) return false;
    if(!writePriming(err)) return false;
    if(!writeCloneHistByCell(err)) return false;

    if(!multi.empty() && !writeMultiMatrix(counts, tally, multi, multiGenes, err)) return false;
    if(!writeSummary(err)) return false;
    return true;
}

void SoloCounter::callCells(const std::vector<SoloRec>& counts,
                            const std::vector<uint32_t>& tally,
                            std::vector<uint32_t>& called) const {
    called.clear();
    if(filter_ == SOLO_FILTER_NONE) return;

    // Total UMIs per barcode. counts is already sorted by (cb, gene), so
    // barcodes arrive in contiguous runs.
    std::vector<std::pair<uint64_t, uint32_t> > perCell;   // (umis, cb)
    for(size_t i = 0; i < counts.size(); ) {
        size_t j = i;
        uint64_t tot = 0;
        while(j < counts.size() && counts[j].cb == counts[i].cb) { tot += tally[j]; j++; }
        perCell.push_back(std::make_pair(tot, counts[i].cb));
        i = j;
    }
    // Descending by count; ties broken by barcode index so the result does not
    // depend on sort stability.
    std::sort(perCell.begin(), perCell.end(),
              [](const std::pair<uint64_t, uint32_t>& a, const std::pair<uint64_t, uint32_t>& b) {
                  if(a.first != b.first) return a.first > b.first;
                  return a.second < b.second;
              });
    if(perCell.empty()) return;

    if(filter_ == SOLO_FILTER_TOPCELLS) {
        size_t n = std::min((size_t)topCells_, perCell.size());
        for(size_t i = 0; i < n; i++) called.push_back(perCell[i].second);
    } else {
        // CellRanger 2.2 knee: take the maxPercentile-th count among the top
        // nExpectedCells barcodes, divide by maxMinRatio, keep everything at or
        // above that. EmptyDrops_CR uses this as its starting point too, so the
        // same code serves both until the Python pass refines it.
        size_t idx = (size_t)((double)expectedCells_ * (1.0 - maxPercentile_));
        if(idx >= perCell.size()) idx = perCell.size() - 1;
        uint64_t umiMax = perCell[idx].first;
        uint64_t umiMin = umiMax / (uint64_t)(maxMinRatio_ > 0 ? maxMinRatio_ : 1);
        if(umiMin < 1) umiMin = 1;
        for(size_t i = 0; i < perCell.size(); i++) {
            if(perCell[i].first >= umiMin) called.push_back(perCell[i].second);
            else break;   // sorted descending
        }
    }
    std::sort(called.begin(), called.end());
}

bool SoloCounter::writeFiltered(const std::vector<SoloRec>& counts,
                                const std::vector<uint32_t>& tally,
                                const std::vector<uint32_t>& called,
                                std::string& err) const {
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string dir = outDir_ + "/" + std::string(featName) + "/filtered";
    if(!makeDirs(dir)) { err = "could not create output directory: " + dir; return false; }

    // Called barcodes become columns 1..N, so the filtered matrix is dense in
    // cells the way CellRanger and STARsolo emit it.
    std::vector<uint32_t> colOf(wl_->size(), 0xffffffffu);
    for(size_t i = 0; i < called.size(); i++) colOf[called[i]] = (uint32_t)i;

    {
        std::ofstream o((dir + "/features.tsv").c_str());
        if(!o.good()) { err = "could not write filtered features.tsv"; return false; }
        for(size_t g = 0; g < gm_->numGenes(); g++)
            o << gm_->geneId(g) << "\t" << gm_->geneName(g) << "\tGene Expression\n";
    }
    {
        std::ofstream o((dir + "/barcodes.tsv").c_str());
        if(!o.good()) { err = "could not write filtered barcodes.tsv"; return false; }
        char buf[40];
        for(size_t i = 0; i < called.size(); i++) {
            soloUnpack(wl_->codeAt(called[i]), wl_->cbLen(), buf);
            o << buf << "\n";
        }
    }
    size_t nz = 0;
    for(size_t k = 0; k < counts.size(); k++)
        if(colOf[counts[k].cb] != 0xffffffffu) nz++;
    {
        std::ofstream o((dir + "/matrix.mtx").c_str());
        if(!o.good()) { err = "could not write filtered matrix.mtx"; return false; }
        o << "%%MatrixMarket matrix coordinate integer general\n%\n";
        o << gm_->numGenes() << " " << called.size() << " " << nz << "\n";
        for(size_t k = 0; k < counts.size(); k++) {
            uint32_t c = colOf[counts[k].cb];
            if(c == 0xffffffffu) continue;
            o << (counts[k].gene() + 1) << " " << (c + 1) << " " << tally[k] << "\n";
        }
    }
    return true;
}

bool SoloCounter::writeMultiMatrix(const std::vector<SoloRec>& counts,
                                   const std::vector<uint32_t>& tally,
                                   std::vector<SoloMultiRec>& multi,
                                   const std::vector<uint32_t>& multiGenes,
                                   std::string& err) {
    // Collapse multi-gene reads to molecules first: the same UMI seen several
    // times is one molecule, and distributing it once per read would inflate
    // multimapper genes relative to unique ones.
    // The gene list must participate in the ordering, not just its length.
    // Two reads of the same molecule can carry different compatible gene sets,
    // and picking whichever landed first would make the result depend on the
    // thread count. Ordering by (fewest genes, then lexicographically) makes
    // the choice deterministic and prefers the most specific assignment.
    std::sort(multi.begin(), multi.end(),
              [&](const SoloMultiRec& a, const SoloMultiRec& b) {
                  if(a.cb != b.cb) return a.cb < b.cb;
                  if(a.umi != b.umi) return a.umi < b.umi;
                  if(a.n != b.n) return a.n < b.n;
                  for(uint32_t k = 0; k < a.n; k++) {
                      uint32_t ga = multiGenes[a.off + k], gb = multiGenes[b.off + k];
                      if(ga != gb) return ga < gb;
                  }
                  return false;
              });
    std::vector<SoloMultiRec> mol;
    for(size_t i = 0; i < multi.size(); ) {
        size_t j = i;
        while(j < multi.size() && multi[j].cb == multi[i].cb && multi[j].umi == multi[i].umi) j++;
        mol.push_back(multi[i]);
        i = j;
    }
    nMultiUMIs_ = mol.size();

    // Unique counts, indexed for fast per-cell lookup.
    std::map<std::pair<uint32_t, uint32_t>, double> acc;
    for(size_t k = 0; k < counts.size(); k++)
        acc[std::make_pair(counts[k].cb, counts[k].gene())] += (double)tally[k];

    // Resolve per cell: EM over one cell's molecules, which is how STARsolo
    // does it and keeps each problem small.
    std::sort(mol.begin(), mol.end(),
              [](const SoloMultiRec& a, const SoloMultiRec& b) { return a.cb < b.cb; });
    for(size_t i = 0; i < mol.size(); ) {
        size_t j = i;
        while(j < mol.size() && mol[j].cb == mol[i].cb) j++;
        uint32_t cb = mol[i].cb;

        if(multiMode_ == SOLO_MULTI_UNIFORM) {
            for(size_t k = i; k < j; k++)
                for(uint32_t g = 0; g < mol[k].n; g++)
                    acc[std::make_pair(cb, multiGenes[mol[k].off + g])] += 1.0 / (double)mol[k].n;
        } else {
            // Gene set for this cell: those with unique support plus any
            // reachable from a multi-gene molecule.
            std::map<uint32_t, double> theta;
            for(size_t k = i; k < j; k++)
                for(uint32_t g = 0; g < mol[k].n; g++)
                    theta[multiGenes[mol[k].off + g]] = 0.0;
            for(std::map<uint32_t, double>::iterator it = theta.begin(); it != theta.end(); ++it) {
                std::map<std::pair<uint32_t, uint32_t>, double>::const_iterator u =
                    acc.find(std::make_pair(cb, it->first));
                // A pseudocount keeps a gene with no unique support reachable;
                // initialising it at zero would pin it there forever.
                it->second = (u == acc.end() ? 0.0 : u->second) + 1e-3;
            }
            for(int iter = 0; iter < 100; iter++) {
                std::map<uint32_t, double> next;
                for(std::map<uint32_t, double>::iterator it = theta.begin(); it != theta.end(); ++it)
                    next[it->first] = 0.0;
                for(size_t k = i; k < j; k++) {
                    double sum = 0.0;
                    for(uint32_t g = 0; g < mol[k].n; g++) sum += theta[multiGenes[mol[k].off + g]];
                    if(sum <= 0.0) continue;
                    for(uint32_t g = 0; g < mol[k].n; g++) {
                        uint32_t gi = multiGenes[mol[k].off + g];
                        next[gi] += theta[gi] / sum;
                    }
                }
                double delta = 0.0;
                for(std::map<uint32_t, double>::iterator it = next.begin(); it != next.end(); ++it) {
                    std::map<std::pair<uint32_t, uint32_t>, double>::const_iterator u =
                        acc.find(std::make_pair(cb, it->first));
                    double base = (u == acc.end() ? 0.0 : u->second);
                    double val = base + it->second + 1e-3;
                    delta += std::abs(val - theta[it->first]);
                    theta[it->first] = val;
                }
                if(delta < 1e-6) break;
            }
            // Fold the converged multimapper mass back onto the unique counts.
            for(size_t k = i; k < j; k++) {
                double sum = 0.0;
                for(uint32_t g = 0; g < mol[k].n; g++) sum += theta[multiGenes[mol[k].off + g]];
                if(sum <= 0.0) continue;
                for(uint32_t g = 0; g < mol[k].n; g++) {
                    uint32_t gi = multiGenes[mol[k].off + g];
                    acc[std::make_pair(cb, gi)] += theta[gi] / sum;
                }
            }
        }
        i = j;
    }

    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string raw = outDir_ + "/" + featName + "/raw";
    if(!makeDirs(raw)) { err = "could not create output directory: " + raw; return false; }
    std::string fn = raw + (multiMode_ == SOLO_MULTI_UNIFORM
                            ? "/UniqueAndMult-Uniform.mtx" : "/UniqueAndMult-EM.mtx");
    std::ofstream o(fn.c_str());
    if(!o.good()) { err = "could not write " + fn; return false; }
    // "real", not "integer": distributing a molecule across genes yields
    // fractional counts, and rounding here would lose the point of doing it.
    o << "%%MatrixMarket matrix coordinate real general\n%\n";
    o << gm_->numGenes() << " " << wl_->size() << " " << acc.size() << "\n";
    o.setf(std::ios::fixed); o.precision(5);
    for(std::map<std::pair<uint32_t, uint32_t>, double>::const_iterator it = acc.begin();
        it != acc.end(); ++it) {
        o << (it->first.second + 1) << " " << (it->first.first + 1) << " " << it->second << "\n";
    }
    return true;
}

bool SoloCounter::writeVelocyto(std::vector<SoloRec>& velo, std::string& err) {
    // Collapse to one class per molecule. STARsolo's rule: a UMI is spliced
    // only if every read supporting it is spliced, unspliced only if every
    // read is unspliced, and ambiguous otherwise.
    std::sort(velo.begin(), velo.end(),
              [](const SoloRec& a, const SoloRec& b) {
                  if(a.cb != b.cb) return a.cb < b.cb;
                  if(a.gene() != b.gene()) return a.gene() < b.gene();
                  return a.umi < b.umi;
              });

    std::vector<SoloRec> keys;
    std::vector<uint32_t> nS, nU, nA;
    size_t i = 0;
    while(i < velo.size()) {
        size_t j = i;
        while(j < velo.size() && velo[j].cb == velo[i].cb && velo[j].gene() == velo[i].gene()) j++;
        uint32_t cS = 0, cU = 0, cA = 0;
        for(size_t k = i; k < j; ) {
            size_t m = k;
            bool anyS = false, anyU = false, anyA = false;
            while(m < j && velo[m].umi == velo[k].umi) {
                uint32_t c = velo[m].veloClass();
                if(c == SOLO_VELO_SPLICED) anyS = true;
                else if(c == SOLO_VELO_UNSPLICED) anyU = true;
                else anyA = true;
                m++;
            }
            if(anyA || (anyS && anyU)) cA++;
            else if(anyU) cU++;
            else cS++;
            k = m;
        }
        SoloRec r = velo[i];
        keys.push_back(r); nS.push_back(cS); nU.push_back(cU); nA.push_back(cA);
        nSpliced_ += cS; nUnspliced_ += cU; nAmbiguous_ += cA;
        i = j;
    }

    std::string raw = outDir_ + "/Velocyto/raw";
    if(!makeDirs(raw)) { err = "could not create output directory: " + raw; return false; }
    {
        std::ofstream o((raw + "/features.tsv").c_str());
        if(!o.good()) { err = "could not write Velocyto features.tsv"; return false; }
        for(size_t g = 0; g < gm_->numGenes(); g++)
            o << gm_->geneId(g) << "\t" << gm_->geneName(g) << "\tGene Expression\n";
    }
    {
        std::ofstream o((raw + "/barcodes.tsv").c_str());
        if(!o.good()) { err = "could not write Velocyto barcodes.tsv"; return false; }
        char buf[40];
        for(size_t b = 0; b < wl_->size(); b++) {
            soloUnpack(wl_->codeAt((uint32_t)b), wl_->cbLen(), buf);
            o << buf << "\n";
        }
    }
    const char* fname[3] = { "spliced.mtx", "unspliced.mtx", "ambiguous.mtx" };
    const std::vector<uint32_t>* vecs[3] = { &nS, &nU, &nA };
    for(int a = 0; a < 3; a++) {
        const std::vector<uint32_t>& N = *vecs[a];
        size_t nz = 0;
        for(size_t k = 0; k < N.size(); k++) if(N[k] > 0) nz++;
        std::ofstream o((raw + "/" + fname[a]).c_str());
        if(!o.good()) { err = std::string("could not write ") + fname[a]; return false; }
        o << "%%MatrixMarket matrix coordinate integer general\n%\n";
        o << gm_->numGenes() << " " << wl_->size() << " " << nz << "\n";
        for(size_t k = 0; k < N.size(); k++) {
            if(N[k] == 0) continue;
            o << (keys[k].gene() + 1) << " " << (keys[k].cb + 1) << " " << N[k] << "\n";
        }
    }
    return true;
}

bool SoloCounter::writeAllelic(std::vector<SoloAllelicRec>& allelic, std::string& err) {
    // Deduplicate by UMI within each (cell, variant, allele), the same rule the
    // gene matrix uses, so an allele is counted once per molecule rather than
    // once per read.
    std::sort(allelic.begin(), allelic.end(),
              [](const SoloAllelicRec& a, const SoloAllelicRec& b) {
                  if(a.cb != b.cb) return a.cb < b.cb;
                  if(a.variantAndAllele != b.variantAndAllele)
                      return a.variantAndAllele < b.variantAndAllele;
                  return a.umi < b.umi;
              });

    std::vector<SoloAllelicRec> keys;   // one per (cell, variant, allele)
    std::vector<uint32_t> refN, altN;
    std::vector<SoloAllelicRec> uniq;
    size_t i = 0;
    while(i < allelic.size()) {
        size_t j = i;
        while(j < allelic.size() && allelic[j].cb == allelic[i].cb &&
              allelic[j].variantAndAllele == allelic[i].variantAndAllele) j++;
        uint32_t n = 0;
        for(size_t k = i; k < j; ) {
            size_t m = k;
            while(m < j && allelic[m].umi == allelic[k].umi) m++;
            n++;
            k = m;
        }
        uniq.push_back(allelic[i]);
        refN.push_back(allelic[i].isAlt() ? 0 : n);
        altN.push_back(allelic[i].isAlt() ? n : 0);
        if(allelic[i].isAlt()) nAltUMIs_ += n; else nRefUMIs_ += n;
        i = j;
    }

    {
        std::vector<uint32_t> vs;
        for(size_t k = 0; k < uniq.size(); k++) vs.push_back(uniq[k].variant());
        std::sort(vs.begin(), vs.end());
        vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
        nVariantsSeen_ = vs.size();
    }

    std::string raw = outDir_ + "/Allelic/raw";
    if(!makeDirs(raw)) { err = "could not create output directory: " + raw; return false; }
    {
        std::ofstream o((raw + "/features.tsv").c_str());
        if(!o.good()) { err = "could not write Allelic features.tsv"; return false; }
        for(size_t v = 0; v < vi_->size(); v++) o << vi_->name((uint32_t)v) << "\tVariant\n";
    }
    {
        std::ofstream o((raw + "/barcodes.tsv").c_str());
        if(!o.good()) { err = "could not write Allelic barcodes.tsv"; return false; }
        char buf[40];
        for(size_t b = 0; b < wl_->size(); b++) {
            soloUnpack(wl_->codeAt((uint32_t)b), wl_->cbLen(), buf);
            o << buf << "\n";
        }
    }
    // Two matrices rather than one with an allele axis, so each can be read by
    // the ordinary 10x readers and subtracted or ratioed directly.
    const char* fname[2] = { "ref.mtx", "alt.mtx" };
    for(int a = 0; a < 2; a++) {
        const std::vector<uint32_t>& N = (a == 0) ? refN : altN;
        size_t nz = 0;
        for(size_t k = 0; k < N.size(); k++) if(N[k] > 0) nz++;
        std::ofstream o((raw + "/" + fname[a]).c_str());
        if(!o.good()) { err = std::string("could not write ") + fname[a]; return false; }
        o << "%%MatrixMarket matrix coordinate integer general\n%\n";
        o << vi_->size() << " " << wl_->size() << " " << nz << "\n";
        for(size_t k = 0; k < N.size(); k++) {
            if(N[k] == 0) continue;
            o << (uniq[k].variant() + 1) << " " << (uniq[k].cb + 1) << " " << N[k] << "\n";
        }
    }
    return true;
}

bool SoloCounter::writeMatrix(const std::vector<SoloRec>& counts,
                              const std::vector<uint32_t>& tally,
                              std::string& err) const {
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string base = outDir_ + "/" + featName;
    std::string raw  = base + "/raw";
    if(!makeDirs(raw)) { err = "could not create output directory: " + raw; return false; }

    // features.tsv -- gene order is the model's order, so matrices from
    // different runs against the same annotation are directly comparable.
    {
        std::ofstream o((raw + "/features.tsv").c_str());
        if(!o.good()) { err = "could not write features.tsv"; return false; }
        for(size_t g = 0; g < gm_->numGenes(); g++) {
            o << gm_->geneId(g) << "\t" << gm_->geneName(g) << "\tGene Expression\n";
        }
    }
    // barcodes.tsv -- the full whitelist, matching STARsolo's raw layout.
    {
        std::ofstream o((raw + "/barcodes.tsv").c_str());
        if(!o.good()) { err = "could not write barcodes.tsv"; return false; }
        char buf[40];
        for(size_t b = 0; b < wl_->size(); b++) {
            soloUnpack(wl_->codeAt((uint32_t)b), wl_->cbLen(), buf);
            o << buf << "\n";
        }
    }
    // matrix.mtx -- rows are features, columns barcodes (CellRanger orientation).
    {
        std::ofstream o((raw + "/matrix.mtx").c_str());
        if(!o.good()) { err = "could not write matrix.mtx"; return false; }
        o << "%%MatrixMarket matrix coordinate integer general\n%\n";
        o << gm_->numGenes() << " " << wl_->size() << " " << counts.size() << "\n";
        for(size_t k = 0; k < counts.size(); k++) {
            o << (counts[k].gene() + 1) << " " << (counts[k].cb + 1) << " " << tally[k] << "\n";
        }
    }
    return true;
}

bool SoloCounter::writeCloneHist(std::string& err) const {
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string path = outDir_ + "/" + featName + "/UMIcloneSize.tsv";
    std::ofstream o(path.c_str());
    if(!o.good()) { err = "could not write " + path; return false; }

    o << "# reads per UMI, over molecules surviving deduplication\n";
    o << "# umis_all: every barcode.  umis_in_cells: called cells only"
      << (filter_ == SOLO_FILTER_NONE ? " (no cell filter ran; column is zero)" : "")
      << "\n";
    o << "reads_per_umi\tumis_all\tumis_in_cells\n";

    size_t upto = cloneAll_.bins.size();
    if(cloneCell_.bins.size() > upto) upto = cloneCell_.bins.size();
    for(size_t k = 1; k < upto; k++) {
        uint64_t a = k < cloneAll_.bins.size()  ? cloneAll_.bins[k]  : 0;
        uint64_t c = k < cloneCell_.bins.size() ? cloneCell_.bins[k] : 0;
        if(a == 0 && c == 0) continue;   // the tail is sparse; skip empty bins
        o << k << "\t" << a << "\t" << c << "\n";
    }
    // Clone sizes past the last exact bin are counted but not resolved.
    if(cloneAll_.nOver > 0 || cloneCell_.nOver > 0) {
        o << ">" << SoloCloneHist::kMaxBin << "\t"
          << cloneAll_.nOver << "\t" << cloneCell_.nOver << "\n";
    }
    return true;
}

bool SoloCounter::writeCloneHistByExpr(std::string& err) const {
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string path = outDir_ + "/" + featName + "/UMIcloneSizeByExpr.tsv";
    std::ofstream o(path.c_str());
    if(!o.good()) { err = "could not write " + path; return false; }

    o << "# reads per UMI in called cells, split by gene-expression decile\n";
    o << "# d1 holds the genes accounting for the first tenth of all molecules\n";
    o << "# genes per decile:";
    for(size_t d = 0; d < kDeciles; d++) o << " d" << (d + 1) << "=" << decileGenes_[d];
    o << "\n# molecules per decile:";
    for(size_t d = 0; d < kDeciles; d++) o << " d" << (d + 1) << "=" << decileUmis_[d];
    o << "\nreads_per_umi";
    for(size_t d = 0; d < kDeciles; d++) o << "\td" << (d + 1);
    o << "\n";

    size_t upto = 0;
    for(size_t d = 0; d < kDeciles; d++)
        if(cloneByDecile_[d].bins.size() > upto) upto = cloneByDecile_[d].bins.size();
    for(size_t k = 1; k < upto; k++) {
        uint64_t row = 0;
        for(size_t d = 0; d < kDeciles; d++)
            if(k < cloneByDecile_[d].bins.size()) row += cloneByDecile_[d].bins[k];
        if(row == 0) continue;
        o << k;
        for(size_t d = 0; d < kDeciles; d++)
            o << "\t" << (k < cloneByDecile_[d].bins.size() ? cloneByDecile_[d].bins[k] : 0);
        o << "\n";
    }
    return true;
}

bool SoloCounter::writePriming(std::string& err) const {
    if(!internalPriming()) return true;
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string path = outDir_ + "/" + featName + "/InternalPriming.tsv";
    std::ofstream o(path.c_str());
    if(!o.good()) { err = "could not write " + path; return false; }

    o << "# molecules whose 3' end sits on an A-rich genomic window ("
      << kPrimeMinCount << "+ A in " << kPrimeWindow << " nt, or a run of "
      << kPrimeMinRun << ")\n";
    o << "# an oligo(dT) primer that annealed inside a transcript rather than to"
         " its poly(A) tail\n";
    o << "gene_id\tgene_name\tumis\tprimed_umis\tfraction\n";

    // Ordered by molecule count: the genes worth looking at first are the ones
    // carrying enough signal for the fraction to mean anything.
    std::vector<uint32_t> order;
    for(size_t g = 0; g < geneUmis_.size(); g++) if(geneUmis_[g]) order.push_back((uint32_t)g);
    std::sort(order.begin(), order.end(), [this](uint32_t a, uint32_t b) {
        if(geneUmis_[a] != geneUmis_[b]) return geneUmis_[a] > geneUmis_[b];
        return a < b;
    });
    o.setf(std::ios::fixed); o.precision(4);
    for(size_t k = 0; k < order.size(); k++) {
        const uint32_t g = order[k];
        o << gm_->geneId(g) << "\t" << gm_->geneName(g) << "\t"
          << geneUmis_[g] << "\t" << genePrimedUmis_[g] << "\t"
          << (double)genePrimedUmis_[g] / (double)geneUmis_[g] << "\n";
    }
    return true;
}

bool SoloCounter::writeCloneHistByCell(std::string& err) const {
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string path = outDir_ + "/" + featName + "/UMIcloneSizeByCell.tsv";
    std::ofstream o(path.c_str());
    if(!o.good()) { err = "could not write " + path; return false; }

    o << "# per-cell clone-size histogram, called cells, in the bins the\n"
         "# PhantomUMI detector consumes: k=1..19 with a k>=20 catch-all.\n"
         "# umi_reads is the cell's total reads over surviving molecules; the\n"
         "# detector's own threshold is >= 10000.\n";
    o << "barcode\tumis\tumi_reads";
    for(int b = 1; b < kCloneBins; b++) o << "\tk" << b;
    o << "\tk" << kCloneBins << "plus\n";

    std::vector<uint64_t> bins(kCloneBins, 0);
    char buf[64];
    for(size_t i = 0; i < cellHist_.size(); i++) {
        const SoloCloneHist& h = cellHist_[i];
        soloCloneBins(h, &bins[0]);
        soloUnpack(wl_->codeAt(cellIdx_[i]), wl_->cbLen(), buf);
        o << buf << "\t" << h.umis() << "\t" << h.reads();
        for(int b = 0; b < kCloneBins; b++) o << "\t" << bins[b];
        o << "\n";
    }
    return true;
}

bool SoloCounter::writeSummary(std::string& err) const {
    const char* featName = (feature() == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
    std::string path = outDir_ + "/" + featName + "/Summary.csv";
    std::ofstream o(path.c_str());
    if(!o.good()) { err = "could not write " + path; return false; }
    double reads = (double)(nReads_ > 0 ? nReads_ : 1);
    o.setf(std::ios::fixed); o.precision(6);
    o << "Number of Reads," << nReads_ << "\n";
    o << "Reads With Valid Barcodes," << (double)(nValidCB_ + nAmbigResolved_) / reads << "\n";
    o << "Reads Mapped to Gene: Unique Gene," << (double)nCounted_ / reads << "\n";
    o << "Reads Unmapped," << (double)nUnmapped_ / reads << "\n";
    o << "Reads With No Gene," << (double)nNoGene_ / reads << "\n";
    o << "Reads Multi-Gene," << (double)nMultiGene_ / reads << "\n";
    o << "Barcodes Ambiguous 1MM," << nAmbigCB_ << "\n";
    o << "Barcodes Ambiguous Resolved," << nAmbigResolved_ << "\n";
    o << "Sequencing Saturation,"
      << (nCounted_ > 0 ? 1.0 - (double)nUMIs_ / (double)nCounted_ : 0.0) << "\n";
    o.unsetf(std::ios::fixed);
    // Saturation is the mean of the reads-per-UMI distribution restated. These
    // two describe its tail, which is where PCR-repriming artefacts show; the
    // shape itself is in UMIcloneSize.tsv.
    o << "UMI Clone Size P99," << cloneAll_.percentile(0.99) << "\n";
    o << "UMI Clone Size Max," << cloneAll_.maxClone << "\n";
    if(internalPriming()) {
        o.setf(std::ios::fixed); o.precision(6);
        o << "Reads Internally Primed,"
          << (nPrimeChecked_ ? (double)nPrimedReads_ / (double)nPrimeChecked_ : 0.0) << "\n";
        o.unsetf(std::ios::fixed);
        o << "UMIs Internally Primed," << nPrimedUmis_ << "\n";
        o.setf(std::ios::fixed); o.precision(6);
        o << "Fraction of UMIs Internally Primed,"
          << (nUMIs_ ? (double)nPrimedUmis_ / (double)nUMIs_ : 0.0) << "\n";
        o.unsetf(std::ios::fixed);
    }
    // Not "Estimated Number of Cells": no cell calling has happened, this is
    // every barcode with at least one UMI. Naming it otherwise would invite
    // comparison against CellRanger's filtered cell count, which is a
    // different quantity entirely.
    o << "Barcodes With UMIs," << nCells_ << "\n";
    o << "Total Genes Detected," << nGenesDetected_ << "\n";
    o << "Total UMIs," << nUMIs_ << "\n";
    if(multiMode_ != SOLO_MULTI_UNIQUE) {
        o << "Multi-Gene UMIs Distributed," << nMultiUMIs_ << "\n";
    }
    if(velocyto_) {
        o << "Velocyto Spliced," << nSpliced_ << "\n";
        o << "Velocyto Unspliced," << nUnspliced_ << "\n";
        o << "Velocyto Ambiguous," << nAmbiguous_ << "\n";
    }
    if(filter_ != SOLO_FILTER_NONE) {
        o << "Estimated Number of Cells," << nCalledCells_ << "\n";
        o << "UMIs in Cells," << nUMIsInCells_ << "\n";
        o.setf(std::ios::fixed); o.precision(6);
        o << "Fraction of UMIs in Cells,"
          << (nUMIs_ > 0 ? (double)nUMIsInCells_ / (double)nUMIs_ : 0.0) << "\n";
        o.unsetf(std::ios::fixed);
    }
    if(vi_ != NULL && !vi_->empty()) {
        o << "Variants Observed," << nVariantsSeen_ << "\n";
        o << "Allelic UMIs: Reference," << nRefUMIs_ << "\n";
        o << "Allelic UMIs: Alternate," << nAltUMIs_ << "\n";
    }
    return true;
}
