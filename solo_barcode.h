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

#ifndef SOLO_BARCODE_H_
#define SOLO_BARCODE_H_

#include <stdint.h>
#include <string>
#include <vector>

/**
 * Cell-barcode and UMI handling for single-cell quantification.
 *
 * HISAT2 has no notion of a barcode, so this is entirely new. Barcodes are
 * 2-bit packed: a 16bp 10x barcode fits exactly in a uint32_t, and UMIs up to
 * 32bp in a uint64_t. Packing keeps the per-read state POD and small, which
 * matters because it is carried on every Read and will later be written to
 * fixed-width count records.
 */

/** Outcome of resolving one read's barcode against the whitelist. */
enum SoloCBStatus {
    SOLO_CB_ABSENT = 0,   // no barcode extracted (solo off, or read too short)
    SOLO_CB_EXACT,        // matched the whitelist exactly
    SOLO_CB_1MM,          // corrected to a unique 1-mismatch neighbour
    SOLO_CB_AMBIG,        // several 1-mismatch neighbours; needs counts to resolve
    SOLO_CB_NOMATCH,      // no whitelist match within one mismatch
    SOLO_CB_HAS_N         // contained a non-ACGT base, so cannot be packed
};

/** How the barcode reaches the aligner. */
enum SoloInputMode {
    SOLO_INPUT_NONE = 0,
    SOLO_INPUT_READNAME,  // "...:_CB_UMI" suffix, the UMI-tools convention
    SOLO_INPUT_MATE       // a mate carries CB+UMI (10x: R1)
};

/** Barcode/UMI geometry and behaviour. */
struct SoloParams {
    SoloInputMode mode;
    int  barcodeMate;   // which mate carries CB+UMI (1 or 2); 2 matches STARsolo
    int  cbStart;       // 1-based offset of the barcode within the barcode read
    int  cbLen;
    int  umiStart;      // 1-based
    int  umiLen;
    bool correct1MM;    // attempt 1-mismatch whitelist correction
    bool emitRaw;       // also emit CR:Z/UR:Z (uncorrected barcode/UMI)

    SoloParams()
        : mode(SOLO_INPUT_NONE), barcodeMate(2),
          cbStart(1), cbLen(16), umiStart(17), umiLen(12),
          correct1MM(true), emitRaw(false) {}

    bool enabled() const { return mode != SOLO_INPUT_NONE; }
};

/** Per-read barcode state. POD; lives on Read and is cleared in reset(). */
struct SoloRead {
    uint32_t cbPacked;   // 2-bit packed barcode as sequenced
    uint32_t cbIdx;      // whitelist index after correction, or kNoIdx
    uint64_t umiPacked;  // 2-bit packed UMI
    uint8_t  cbLen;      // 0 when no barcode was extracted
    uint8_t  umiLen;
    uint8_t  status;     // SoloCBStatus
    uint8_t  cbMinQual;  // lowest base quality across the barcode

    static const uint32_t kNoIdx = 0xffffffffu;

    void reset() {
        cbPacked = 0; cbIdx = kNoIdx; umiPacked = 0;
        cbLen = 0; umiLen = 0; status = SOLO_CB_ABSENT; cbMinQual = 0;
    }
    bool hasBarcode() const { return cbLen > 0; }
    bool corrected()  const { return status == SOLO_CB_EXACT || status == SOLO_CB_1MM; }
};

/** Packs ACGT to 2 bits. Returns false if any base is not ACGT. */
bool soloPack(const char* s, int len, uint64_t& out);
/** Unpacks into a caller-supplied buffer of at least len+1 bytes. */
void soloUnpack(uint64_t v, int len, char* out);

/**
 * The 10x barcode whitelist -- ~3.7M entries for 3'/5' v3.
 *
 * Stored as a sorted vector of packed barcodes plus a 65,536-entry table of
 * bucket offsets keyed on the top 16 bits, so a lookup is one table read plus
 * a binary search over a few dozen entries. ~15MB for 3.7M barcodes, and
 * immutable after load, so lookups need no locking.
 */
class SoloWhitelist {
public:
    SoloWhitelist() : loaded_(false), cbLen_(0) {}

    bool load(const std::string& path, int cbLen, std::string& err);
    bool loaded() const { return loaded_; }
    size_t size() const { return codes_.size(); }
    int cbLen() const { return cbLen_; }

    /** Packed barcode at a whitelist index; for printing a corrected CB. */
    uint32_t codeAt(uint32_t idx) const { return codes_[idx]; }

    /** Index of an exact match, or SoloRead::kNoIdx. */
    uint32_t find(uint32_t code) const;

    /**
     * Resolves a barcode, filling status and cbIdx.
     *
     * Exact matches win. Otherwise, if correction is enabled, every
     * single-substitution neighbour is probed; a unique hit is accepted as
     * SOLO_CB_1MM and several as SOLO_CB_AMBIG. CellRanger disambiguates the
     * ambiguous case with a posterior weighted by each candidate's observed
     * count, which is not knowable in a streaming pass -- that resolution is
     * deferred to the counting stage, so the raw barcode is retained here.
     */
    void resolve(SoloRead& sr, bool correct1MM) const;

private:
    std::vector<uint32_t> codes_;    // sorted, unique
    std::vector<uint32_t> prefix_;   // 65537 entries: bucket starts by top 16 bits
    bool loaded_;
    int  cbLen_;
};

/**
 * Extracts CB+UMI from a barcode read's sequence.
 * Returns false if the read is too short for the configured geometry.
 */
bool soloExtractFromSeq(const char* seq, const char* qual, size_t len,
                        const SoloParams& p, SoloRead& sr);

/**
 * Extracts CB+UMI from a read name of the form "<name>_<CB>_<UMI>"
 * (the umi_tools convention). Returns false if the name has no such suffix.
 */
bool soloExtractFromName(const char* name, size_t len,
                         const SoloParams& p, SoloRead& sr);

#endif /* SOLO_BARCODE_H_ */
