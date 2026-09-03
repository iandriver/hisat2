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

#ifndef SOLO_PRIMING_H_
#define SOLO_PRIMING_H_

#include <stdint.h>
#include <stddef.h>

/*
 * Internal priming.
 *
 * An oligo(dT) primer is supposed to find the poly(A) tail. It will also anneal
 * to any A-rich stretch inside a transcript, and the resulting molecule looks
 * like a genuine 3' end: it carries a barcode and a UMI, it maps, and it is
 * counted. Nothing downstream of the count matrix can tell it apart from a real
 * one -- the only evidence is the genome immediately past where the read ends,
 * which is A-rich for a mispriming event and not for a real polyadenylation
 * site, whose A's live on the transcript rather than in the reference.
 *
 * The thresholds are the ones the field converged on (CellRanger, scAPA and
 * Sierra all use this shape): a 20 nt window read outward from the alignment's
 * 3' end, called primed at 12 or more A's, or a run of 6.
 */

const int kPrimeWindow   = 20;   // nt of genome to read past the 3' end
const int kPrimeMinCount = 12;   // A's in the window that call it primed
const int kPrimeMinRun   = 6;    // or this many consecutive A's

/**
 * True if a decoded reference window looks like an internal priming site.
 *
 * `w` holds 2-bit base codes (A=0, C=1, G=2, T=3, N=4) as BitPairReference
 * returns them. `base` is the code to count: A for an alignment on the forward
 * strand, T for one on the reverse, since the transcript's A's are the
 * reference's T's there.
 *
 * A window truncated by the end of a contig is judged on the same density
 * rather than being waved through, which would make every contig edge look
 * clean.
 */
inline bool soloIsPrimingWindow(const char* w, size_t n, int base) {
    // A window shorter than the run threshold cannot be evidence of anything;
    // scaling the count rule down that far would flag a single A at a contig
    // edge.
    if(n < (size_t)kPrimeMinRun) return false;
    int count = 0, run = 0, best = 0;
    for(size_t i = 0; i < n; i++) {
        if(w[i] == (char)base) { count++; if(++run > best) best = run; }
        else run = 0;
    }
    if(best >= kPrimeMinRun) return true;
    // Scale the count threshold to a short window so density, not length,
    // decides. Rounds up: 12/20 of a 10 nt window is 6. Windows below
    // kPrimeMinRun already returned above.
    const int need = (int)((kPrimeMinCount * (int)n + kPrimeWindow - 1) / kPrimeWindow);
    return count >= need;
}

#endif /* SOLO_PRIMING_H_ */
