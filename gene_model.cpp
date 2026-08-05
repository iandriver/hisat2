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

#include "gene_model.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>

namespace {

bool ivStartLess(const GeneIv& a, const GeneIv& b) {
    if(a.start != b.start) return a.start < b.start;
    return a.end < b.end;
}

std::string joinErr(const std::string& path, size_t lineno,
                    const std::string& msg) {
    std::ostringstream os;
    os << path << ":" << lineno << ": " << msg;
    return os.str();
}

} // namespace

void GeneIvIndex::build() {
    // Hoist the few very long intervals so they do not inflate maxEnd_ for
    // everything preceding them.
    std::vector<GeneIv> shortIvs;
    shortIvs.reserve(ivs_.size());
    for(size_t i = 0; i < ivs_.size(); i++) {
        if(ivs_[i].end - ivs_[i].start >= kLongIv) longIvs_.push_back(ivs_[i]);
        else shortIvs.push_back(ivs_[i]);
    }
    ivs_.swap(shortIvs);

    std::sort(ivs_.begin(), ivs_.end(), ivStartLess);
    maxEnd_.resize(ivs_.size());
    int64_t running = INT64_MIN;
    for(size_t i = 0; i < ivs_.size(); i++) {
        if(ivs_[i].end > running) running = ivs_[i].end;
        maxEnd_[i] = running;
    }
}

void GeneIvIndex::query(int64_t s, int64_t e, std::vector<uint32_t>& out) const {
    lastScanned = 0;
    if(e <= s) return;

    for(size_t i = 0; i < longIvs_.size(); i++) {
        lastScanned++;
        if(longIvs_[i].start < e && longIvs_[i].end > s) out.push_back(longIvs_[i].gene);
    }
    if(ivs_.empty()) return;

    // First index whose start >= e; nothing at or after it can overlap.
    GeneIv probe;
    probe.start = e; probe.end = INT64_MIN; probe.gene = 0;
    size_t hi = (size_t)(std::lower_bound(ivs_.begin(), ivs_.end(), probe, ivStartLess)
                         - ivs_.begin());

    // Walk backwards while some interval at or before i can still reach past s.
    for(size_t i = hi; i-- > 0; ) {
        lastScanned++;
        if(maxEnd_[i] <= s) break;              // nothing earlier can overlap
        if(ivs_[i].end > s) out.push_back(ivs_[i].gene);
    }
}

bool GeneModel::load(const std::string& path,
                     const std::vector<std::string>& refnames,
                     const std::vector<int64_t>& reflens,
                     std::string& err) {
    std::ifstream in(path.c_str());
    if(!in.good()) { err = "could not open gene model: " + path; return false; }

    std::map<std::string, int32_t> refmap;
    for(size_t i = 0; i < refnames.size(); i++) {
        // Index reference names may carry a description after whitespace; the
        // GTF only ever names the first token, so key on that.
        const std::string& rn = refnames[i];
        size_t sp = rn.find_first_of(" \t");
        refmap[sp == std::string::npos ? rn : rn.substr(0, sp)] = (int32_t)i;
    }

    // Gene records are read in file order, which hisat2_extract_genes.py emits
    // sorted, so gene index == file order == the order used for features.tsv.
    std::vector<std::vector<GeneIv> >   pendingExons;
    std::vector<std::vector<GeneJunc> > pendingJuncs;
    bool sawMagic = false;
    size_t lineno = 0;
    std::string line;

    while(std::getline(in, line)) {
        lineno++;
        if(line.empty()) continue;

        if(line[0] == '#') {
            if(!sawMagic) {
                if(line.compare(0, 19, "#hisat2-gene-model\t") != 0) {
                    err = joinErr(path, lineno, "not a .ht2gm file (bad magic)");
                    return false;
                }
                sawMagic = true;
                continue;
            }
            if(line.compare(0, 5, "#ref\t") == 0) {
                std::istringstream ls(line);
                std::string tag, name, lenStr;
                std::getline(ls, tag, '\t');
                std::getline(ls, name, '\t');
                std::getline(ls, lenStr, '\t');
                std::map<std::string, int32_t>::const_iterator it = refmap.find(name);
                if(it == refmap.end()) {
                    err = joinErr(path, lineno,
                        "gene model references sequence '" + name +
                        "' which is not in the index -- the model and the index "
                        "were built from different assemblies");
                    return false;
                }
                if(lenStr != "-" && !lenStr.empty() && it->second < (int32_t)reflens.size()) {
                    int64_t want = strtoll(lenStr.c_str(), NULL, 10);
                    int64_t have = reflens[it->second];
                    if(want != have) {
                        std::ostringstream os;
                        os << "sequence '" << name << "' is " << want
                           << " bp in the gene model but " << have
                           << " bp in the index -- different assemblies";
                        err = joinErr(path, lineno, os.str());
                        return false;
                    }
                }
            }
            continue;
        }
        if(!sawMagic) {
            err = joinErr(path, lineno, "not a .ht2gm file (missing header)");
            return false;
        }

        std::istringstream ls(line);
        std::string tag;
        std::getline(ls, tag, '\t');

        if(tag == "G") {
            std::string idxStr, gid, gname, chrom, strand, s, e;
            std::getline(ls, idxStr, '\t'); std::getline(ls, gid,    '\t');
            std::getline(ls, gname,  '\t'); std::getline(ls, chrom,  '\t');
            std::getline(ls, strand, '\t'); std::getline(ls, s,      '\t');
            std::getline(ls, e,      '\t');
            if(e.empty()) { err = joinErr(path, lineno, "malformed G record"); return false; }

            uint32_t idx = (uint32_t)strtoul(idxStr.c_str(), NULL, 10);
            if(idx != genes_.size()) {
                err = joinErr(path, lineno, "G records must be numbered consecutively from 0");
                return false;
            }
            std::map<std::string, int32_t>::const_iterator it = refmap.find(chrom);
            if(it == refmap.end()) {
                err = joinErr(path, lineno, "unknown sequence '" + chrom + "'");
                return false;
            }
            GeneRec g;
            g.idOff   = (uint32_t)nameBlob_.size();
            nameBlob_.insert(nameBlob_.end(), gid.begin(), gid.end());
            nameBlob_.push_back('\0');
            g.nameOff = (uint32_t)nameBlob_.size();
            nameBlob_.insert(nameBlob_.end(), gname.begin(), gname.end());
            nameBlob_.push_back('\0');
            g.refid     = it->second;
            g.fw        = (strand == "+") ? 1 : 0;
            g.bodyStart = strtoll(s.c_str(), NULL, 10);
            g.bodyEnd   = strtoll(e.c_str(), NULL, 10);
            genes_.push_back(g);
            pendingExons.push_back(std::vector<GeneIv>());
            pendingJuncs.push_back(std::vector<GeneJunc>());
        } else if(tag == "E" || tag == "J") {
            std::string idxStr, chrom, a, b;
            std::getline(ls, idxStr, '\t'); std::getline(ls, chrom, '\t');
            std::getline(ls, a, '\t');      std::getline(ls, b, '\t');
            if(b.empty()) { err = joinErr(path, lineno, "malformed " + tag + " record"); return false; }
            uint32_t idx = (uint32_t)strtoul(idxStr.c_str(), NULL, 10);
            if(idx >= genes_.size()) {
                err = joinErr(path, lineno, tag + " record precedes its G record");
                return false;
            }
            if(tag == "E") {
                GeneIv iv;
                iv.start = strtoll(a.c_str(), NULL, 10);
                iv.end   = strtoll(b.c_str(), NULL, 10);
                iv.gene  = idx;
                pendingExons[idx].push_back(iv);
            } else {
                GeneJunc j;
                j.donor    = strtoll(a.c_str(), NULL, 10);
                j.acceptor = strtoll(b.c_str(), NULL, 10);
                j.gene     = idx;
                pendingJuncs[idx].push_back(j);
            }
        }
        // Unknown record types are ignored so the format can gain new ones.
    }

    if(genes_.empty()) { err = "gene model contains no genes: " + path; return false; }

    // Flatten into CSR, and feed the per-reference interval indexes.
    size_t nref = refnames.size();
    exonIdx_.resize(nref);
    bodyIdx_.resize(nref);
    exonOff_.reserve(genes_.size() + 1);
    juncOff_.reserve(genes_.size() + 1);

    for(size_t g = 0; g < genes_.size(); g++) {
        exonOff_.push_back((uint32_t)exons_.size());
        std::vector<GeneIv>& ge = pendingExons[g];
        std::sort(ge.begin(), ge.end(), ivStartLess);
        int32_t refid = genes_[g].refid;
        for(size_t i = 0; i < ge.size(); i++) {
            exons_.push_back(ge[i]);
            exonIdx_[refid].add(ge[i]);
        }
        juncOff_.push_back((uint32_t)juncs_.size());
        std::vector<GeneJunc>& gj = pendingJuncs[g];
        for(size_t i = 0; i < gj.size(); i++) juncs_.push_back(gj[i]);

        GeneIv body;
        body.start = genes_[g].bodyStart;
        body.end   = genes_[g].bodyEnd;
        body.gene  = (uint32_t)g;
        bodyIdx_[refid].add(body);
    }
    exonOff_.push_back((uint32_t)exons_.size());
    juncOff_.push_back((uint32_t)juncs_.size());

    for(size_t i = 0; i < nref; i++) { exonIdx_[i].build(); bodyIdx_[i].build(); }

    loaded_ = true;
    return true;
}

void GeneModel::overlapExon(int32_t refid, int64_t s, int64_t e,
                            std::vector<uint32_t>& out) const {
    if(refid < 0 || (size_t)refid >= exonIdx_.size()) return;
    exonIdx_[refid].query(s, e, out);
}

void GeneModel::overlapBody(int32_t refid, int64_t s, int64_t e,
                            std::vector<uint32_t>& out) const {
    if(refid < 0 || (size_t)refid >= bodyIdx_.size()) return;
    bodyIdx_[refid].query(s, e, out);
}

bool GeneModel::blockInExonsOf(uint32_t g, int64_t s, int64_t e) const {
    for(uint32_t i = exonOff_[g]; i < exonOff_[g + 1]; i++) {
        if(exons_[i].start <= s && exons_[i].end >= e) return true;
        if(exons_[i].start >= e) break;   // exons are sorted by start
    }
    return false;
}

bool GeneModel::junctionInGene(uint32_t g, int64_t donor, int64_t acceptor) const {
    for(uint32_t i = juncOff_[g]; i < juncOff_[g + 1]; i++) {
        if(juncs_[i].donor == donor && juncs_[i].acceptor == acceptor) return true;
    }
    return false;
}

bool GeneModel::strandOk(uint32_t g, bool readFw, GeneStrand strand) const {
    if(strand == GENE_STRAND_UNSTRANDED) return true;
    bool geneFw = genes_[g].fw != 0;
    return (strand == GENE_STRAND_FORWARD) ? (readFw == geneFw) : (readFw != geneFw);
}

size_t GeneModel::assignGenes(int32_t refid,
                              const std::vector<std::pair<int64_t, int64_t> >& blocks,
                              bool readFw,
                              GeneFeature feat,
                              GeneStrand strand,
                              std::vector<uint32_t>& out,
                              std::vector<uint32_t>& scratch) const {
    size_t before = out.size();
    if(blocks.empty()) return 0;

    // Candidates from the first block, then intersect with each later block.
    // Intersecting rather than unioning is what makes this containment: a read
    // straddling two genes ends up assigned to neither.
    scratch.clear();
    if(feat == GENE_FEATURE_EXONIC) overlapExon(refid, blocks[0].first, blocks[0].second, scratch);
    else                            overlapBody(refid, blocks[0].first, blocks[0].second, scratch);
    std::sort(scratch.begin(), scratch.end());
    scratch.erase(std::unique(scratch.begin(), scratch.end()), scratch.end());
    if(scratch.empty()) return 0;

    std::vector<uint32_t> cand;
    cand.swap(scratch);

    for(size_t b = 1; b < blocks.size() && !cand.empty(); b++) {
        scratch.clear();
        if(feat == GENE_FEATURE_EXONIC) overlapExon(refid, blocks[b].first, blocks[b].second, scratch);
        else                            overlapBody(refid, blocks[b].first, blocks[b].second, scratch);
        std::sort(scratch.begin(), scratch.end());
        scratch.erase(std::unique(scratch.begin(), scratch.end()), scratch.end());

        std::vector<uint32_t> keep;
        keep.reserve(cand.size() < scratch.size() ? cand.size() : scratch.size());
        std::set_intersection(cand.begin(), cand.end(),
                              scratch.begin(), scratch.end(),
                              std::back_inserter(keep));
        cand.swap(keep);
    }

    for(size_t i = 0; i < cand.size(); i++) {
        uint32_t g = cand[i];
        if(!strandOk(g, readFw, strand)) continue;
        if(feat == GENE_FEATURE_EXONIC) {
            // Overlap is not enough: every block must sit inside one exon.
            bool ok = true;
            for(size_t b = 0; b < blocks.size(); b++) {
                if(!blockInExonsOf(g, blocks[b].first, blocks[b].second)) { ok = false; break; }
            }
            if(!ok) continue;
        } else {
            if(blocks.front().first < genes_[g].bodyStart ||
               blocks.back().second > genes_[g].bodyEnd) continue;
        }
        out.push_back(g);
    }
    return out.size() - before;
}
