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

namespace {

bool byCbUmiGene(const SoloRec& a, const SoloRec& b) {
    if(a.cb != b.cb)   return a.cb < b.cb;
    if(a.umi != b.umi) return a.umi < b.umi;
    return a.gene() < b.gene();
}

bool byCbGeneUmi(const SoloRec& a, const SoloRec& b) {
    if(a.cb != b.cb)             return a.cb < b.cb;
    if(a.gene() != b.gene())     return a.gene() < b.gene();
    return a.umi < b.umi;
}

/** True if two packed UMIs of the given length differ at exactly one base. */
inline bool umiWithin1(uint64_t a, uint64_t b, int len) {
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

SoloCounter::~SoloCounter() {
    for(size_t i = 0; i < threads_.size(); i++) delete threads_[i];
}

void SoloCounter::init(const GeneModel* gm,
                       const SoloWhitelist* wl,
                       const SoloParams* params,
                       GeneFeature feature,
                       GeneStrand strand,
                       SoloUmiDedup dedup,
                       const std::string& outDir) {
    gm_ = gm; wl_ = wl; params_ = params;
    feature_ = feature; strand_ = strand; dedup_ = dedup;
    outDir_ = outDir;
}

void SoloCounter::reserveThreads(size_t n) {
    for(size_t i = threads_.size(); i < n; i++) threads_.push_back(new SoloCounterThread(this));
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

    // Union the gene assignments over every candidate alignment. If they all
    // agree on one gene the read is unique-gene even when it is multi-locus,
    // which is what STARsolo counts by default.
    static thread_local std::vector<uint32_t> genes, scratch, blockGenes;
    static thread_local std::vector<std::pair<int64_t, int64_t> > blocks;
    static thread_local StackedAln staln;
    genes.clear();

    for(size_t i = 0; i < nresults; i++) {
        const AlnRes& rs = (*results)[i];
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

        blockGenes.clear();
        p->geneModel()->assignGenes((int32_t)rs.refid(), blocks, rs.fw(),
                                    p->feature(), p->strand(), blockGenes, scratch);
        for(size_t g = 0; g < blockGenes.size(); g++) genes.push_back(blockGenes[g]);
    }

    std::sort(genes.begin(), genes.end());
    genes.erase(std::unique(genes.begin(), genes.end()), genes.end());

    if(genes.empty())      { nNoGene_++;    return; }
    if(genes.size() > 1)   { nMultiGene_++; return; }   // EM handling is F6

    nCounted_++;
    if(sr.status == SOLO_CB_AMBIG) {
        SoloAmbigRec ar;
        ar.cbRaw = sr.cbPacked;
        ar.geneAndFlags = genes[0];
        ar.umi = sr.umiPacked;
        ambig_.push_back(ar);
    } else {
        SoloRec r;
        r.cb = sr.cbIdx;
        r.geneAndFlags = genes[0];
        r.umi = sr.umiPacked;
        recs_.push_back(r);
    }
}

bool SoloCounter::finalize(std::string& err) {
    if(!enabled()) { err = "solo counter not initialised"; return false; }

    std::vector<SoloRec> all;
    std::vector<SoloAmbigRec> ambig;
    size_t total = 0;
    for(size_t i = 0; i < threads_.size(); i++) total += threads_[i]->recs_.size();
    all.reserve(total);
    for(size_t i = 0; i < threads_.size(); i++) {
        SoloCounterThread* t = threads_[i];
        all.insert(all.end(), t->recs_.begin(), t->recs_.end());
        ambig.insert(ambig.end(), t->ambig_.begin(), t->ambig_.end());
        nReads_ += t->nReads_; nValidCB_ += t->nValidCB_; nAmbigCB_ += t->nAmbigCB_;
        nNoCB_ += t->nNoCB_;   nUnmapped_ += t->nUnmapped_;
        nNoGene_ += t->nNoGene_; nMultiGene_ += t->nMultiGene_; nCounted_ += t->nCounted_;
        // Release each arena as it is consumed, so peak memory is roughly one
        // copy of the records rather than two.
        std::vector<SoloRec>().swap(t->recs_);
        std::vector<SoloAmbigRec>().swap(t->ambig_);
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
    std::vector<SoloRec> resolved;
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
            SoloRec r = all[i];
            r.setGene(bestGene);
            resolved.push_back(r);
        } else {
            nMultiGene_++;
        }
        i = j;
    }
    std::vector<SoloRec>().swap(all);

    // UMI collapsing within each (cell, gene).
    std::sort(resolved.begin(), resolved.end(), byCbGeneUmi);
    std::vector<SoloRec> counts;    // one entry per (cell, gene)
    std::vector<uint32_t> tally;    // its UMI count
    int umiLen = params_ != NULL ? params_->umiLen : 12;

    i = 0;
    while(i < resolved.size()) {
        size_t j = i;
        while(j < resolved.size() && resolved[j].cb == resolved[i].cb &&
              resolved[j].gene() == resolved[i].gene()) j++;

        // Distinct UMIs in this block, with read support.
        std::vector<uint64_t> umis;
        std::vector<uint32_t> reads;
        for(size_t k = i; k < j; ) {
            size_t m = k;
            while(m < j && resolved[m].umi == resolved[k].umi) m++;
            umis.push_back(resolved[k].umi);
            reads.push_back((uint32_t)(m - k));
            k = m;
        }

        uint32_t n = 0;
        if(dedup_ == SOLO_UMI_NODEDUP) {
            n = (uint32_t)(j - i);
        } else if(dedup_ == SOLO_UMI_EXACT) {
            n = (uint32_t)umis.size();
        } else {
            // Blocks are tiny (a handful of UMIs per gene per cell), so the
            // quadratic neighbour search here is cheaper than building an index.
            std::vector<char> merged(umis.size(), 0);
            for(size_t a = 0; a < umis.size(); a++) {
                if(merged[a]) continue;
                for(size_t b = 0; b < umis.size(); b++) {
                    if(a == b || merged[b]) continue;
                    if(!umiWithin1(umis[a], umis[b], umiLen)) continue;
                    if(dedup_ == SOLO_UMI_1MM_ALL) {
                        if(b > a) merged[b] = 1;
                    } else {
                        // 1MM_CR: the less-supported UMI is assumed to be a
                        // sequencing error of the more-supported one.
                        if(reads[b] < reads[a] || (reads[b] == reads[a] && b > a)) merged[b] = 1;
                    }
                }
            }
            for(size_t a = 0; a < umis.size(); a++) if(!merged[a]) n++;
        }

        if(n > 0) {
            SoloRec r = resolved[i];
            counts.push_back(r);
            tally.push_back(n);
            nUMIs_ += n;
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
    if(!writeSummary(err)) return false;
    return true;
}

bool SoloCounter::writeMatrix(const std::vector<SoloRec>& counts,
                              const std::vector<uint32_t>& tally,
                              std::string& err) const {
    const char* featName = (feature_ == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
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

bool SoloCounter::writeSummary(std::string& err) const {
    const char* featName = (feature_ == GENE_FEATURE_BODY) ? "GeneFull" : "Gene";
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
    // Not "Estimated Number of Cells": no cell calling has happened, this is
    // every barcode with at least one UMI. Naming it otherwise would invite
    // comparison against CellRanger's filtered cell count, which is a
    // different quantity entirely.
    o << "Barcodes With UMIs," << nCells_ << "\n";
    o << "Total Genes Detected," << nGenesDetected_ << "\n";
    o << "Total UMIs," << nUMIs_ << "\n";
    return true;
}
