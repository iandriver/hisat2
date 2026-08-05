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

#include "solo_barcode.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

namespace {

inline int baseCode(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

const char kBases[4] = {'A', 'C', 'G', 'T'};

} // namespace

bool soloPack(const char* s, int len, uint64_t& out) {
    uint64_t v = 0;
    for(int i = 0; i < len; i++) {
        int b = baseCode(s[i]);
        if(b < 0) return false;
        v = (v << 2) | (uint64_t)b;
    }
    out = v;
    return true;
}

void soloUnpack(uint64_t v, int len, char* out) {
    for(int i = len - 1; i >= 0; i--) {
        out[i] = kBases[v & 3];
        v >>= 2;
    }
    out[len] = '\0';
}

bool SoloWhitelist::load(const std::string& path, int cbLen, std::string& err) {
    if(cbLen <= 0 || cbLen > 16) {
        err = "barcode length must be between 1 and 16 for the packed whitelist";
        return false;
    }
    std::ifstream in(path.c_str());
    if(!in.good()) { err = "could not open barcode whitelist: " + path; return false; }

    cbLen_ = cbLen;
    codes_.clear();
    std::string line;
    size_t lineno = 0, badLen = 0, badBase = 0;
    while(std::getline(in, line)) {
        lineno++;
        // Tolerate CRLF and trailing whitespace; some whitelists ship with a
        // second column (e.g. a translated barcode), so stop at whitespace.
        size_t end = line.find_first_of(" \t\r\n");
        if(end != std::string::npos) line.resize(end);
        if(line.empty()) continue;
        if((int)line.size() != cbLen) { badLen++; continue; }
        uint64_t v;
        if(!soloPack(line.c_str(), cbLen, v)) { badBase++; continue; }
        codes_.push_back((uint32_t)v);
    }
    if(codes_.empty()) {
        std::ostringstream os;
        os << "no usable barcodes in " << path << " (" << lineno << " lines read, "
           << badLen << " wrong length, " << badBase << " with non-ACGT bases); "
           << "check that --solo-cb-len matches the whitelist";
        err = os.str();
        return false;
    }
    std::sort(codes_.begin(), codes_.end());
    codes_.erase(std::unique(codes_.begin(), codes_.end()), codes_.end());

    // Bucket by the top 16 bits so a lookup binary-searches only a few dozen
    // entries. For a 16bp barcode that is the first 8 bases.
    prefix_.assign(65537, 0);
    int shift = (cbLen_ * 2 > 16) ? (cbLen_ * 2 - 16) : 0;
    for(size_t i = 0; i < codes_.size(); i++) {
        uint32_t b = (uint32_t)(codes_[i] >> shift);
        if(b > 65535) b = 65535;
        prefix_[b + 1]++;
    }
    for(size_t i = 1; i < prefix_.size(); i++) prefix_[i] += prefix_[i - 1];

    loaded_ = true;
    return true;
}

uint32_t SoloWhitelist::find(uint32_t code) const {
    if(!loaded_) return SoloRead::kNoIdx;
    int shift = (cbLen_ * 2 > 16) ? (cbLen_ * 2 - 16) : 0;
    uint32_t b = code >> shift;
    if(b > 65535) b = 65535;
    std::vector<uint32_t>::const_iterator lo = codes_.begin() + prefix_[b];
    std::vector<uint32_t>::const_iterator hi = codes_.begin() + prefix_[b + 1];
    std::vector<uint32_t>::const_iterator it = std::lower_bound(lo, hi, code);
    if(it != hi && *it == code) return (uint32_t)(it - codes_.begin());
    return SoloRead::kNoIdx;
}

void SoloWhitelist::resolve(SoloRead& sr, bool correct1MM) const {
    if(!sr.hasBarcode()) { sr.status = SOLO_CB_ABSENT; return; }
    if(!loaded_) {
        // With no whitelist every barcode is taken at face value; the packed
        // code itself is the identity.
        sr.cbIdx = sr.cbPacked;
        sr.status = SOLO_CB_EXACT;
        return;
    }
    uint32_t hit = find(sr.cbPacked);
    if(hit != SoloRead::kNoIdx) {
        sr.cbIdx = hit;
        sr.status = SOLO_CB_EXACT;
        return;
    }
    if(!correct1MM) { sr.cbIdx = SoloRead::kNoIdx; sr.status = SOLO_CB_NOMATCH; return; }

    uint32_t found = SoloRead::kNoIdx;
    int nfound = 0;
    for(int pos = 0; pos < sr.cbLen; pos++) {
        int shift = 2 * (sr.cbLen - 1 - pos);
        uint32_t cur = (sr.cbPacked >> shift) & 3u;
        for(uint32_t alt = 0; alt < 4; alt++) {
            if(alt == cur) continue;
            uint32_t cand = (sr.cbPacked & ~(3u << shift)) | (alt << shift);
            uint32_t h = find(cand);
            if(h != SoloRead::kNoIdx) {
                nfound++;
                if(nfound == 1) found = h;
                else if(nfound > 1) {
                    // More than one neighbour: CellRanger picks by posterior
                    // over observed counts, which needs a completed pass.
                    sr.cbIdx = SoloRead::kNoIdx;
                    sr.status = SOLO_CB_AMBIG;
                    return;
                }
            }
        }
    }
    if(nfound == 1) { sr.cbIdx = found; sr.status = SOLO_CB_1MM; }
    else            { sr.cbIdx = SoloRead::kNoIdx; sr.status = SOLO_CB_NOMATCH; }
}

bool soloExtractFromSeq(const char* seq, const char* qual, size_t len,
                        const SoloParams& p, SoloRead& sr) {
    sr.reset();
    size_t cbEnd  = (size_t)p.cbStart  - 1 + (size_t)p.cbLen;
    size_t umiEnd = (size_t)p.umiStart - 1 + (size_t)p.umiLen;
    if(len < cbEnd || len < umiEnd) return false;

    uint64_t cb = 0, umi = 0;
    bool cbOk  = soloPack(seq + p.cbStart  - 1, p.cbLen,  cb);
    bool umiOk = soloPack(seq + p.umiStart - 1, p.umiLen, umi);

    sr.cbLen  = (uint8_t)p.cbLen;
    sr.umiLen = (uint8_t)p.umiLen;
    if(!cbOk || !umiOk) {
        // An N anywhere in the barcode or UMI makes it unusable; recording the
        // reason separately keeps it out of the "no whitelist match" bucket,
        // which would otherwise look like a chemistry mismatch.
        sr.status = SOLO_CB_HAS_N;
        sr.cbLen = 0;
        return false;
    }
    sr.cbPacked  = (uint32_t)cb;
    sr.umiPacked = umi;

    uint8_t minq = 255;
    if(qual != NULL) {
        for(int i = 0; i < p.cbLen; i++) {
            uint8_t q = (uint8_t)qual[p.cbStart - 1 + i];
            if(q < minq) minq = q;
        }
    }
    sr.cbMinQual = (minq == 255) ? 0 : minq;
    return true;
}

bool soloExtractFromName(const char* name, size_t len,
                         const SoloParams& p, SoloRead& sr) {
    sr.reset();
    // umi_tools writes "<original name>_<CB>_<UMI>"; the barcode and UMI are
    // the last two underscore-separated fields.
    size_t last = std::string::npos, prev = std::string::npos;
    for(size_t i = len; i-- > 0; ) {
        if(name[i] == '_') {
            if(last == std::string::npos) last = i;
            else { prev = i; break; }
        }
    }
    if(last == std::string::npos || prev == std::string::npos) return false;

    size_t cbLen  = last - prev - 1;
    size_t umiLen = len - last - 1;
    if(cbLen == 0 || umiLen == 0 || cbLen > 16 || umiLen > 32) return false;

    uint64_t cb = 0, umi = 0;
    if(!soloPack(name + prev + 1, (int)cbLen, cb) ||
       !soloPack(name + last + 1, (int)umiLen, umi)) {
        sr.status = SOLO_CB_HAS_N;
        return false;
    }
    sr.cbPacked  = (uint32_t)cb;
    sr.umiPacked = umi;
    sr.cbLen     = (uint8_t)cbLen;
    sr.umiLen    = (uint8_t)umiLen;
    sr.cbMinQual = 0;   // quality is not preserved by the name convention
    return true;
}
