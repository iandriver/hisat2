/*
 * Unit test for the SAM->BAM converter.
 *
 * Covers the cases HISAT2's own output does not reach, and which would
 * otherwise never be exercised: float and array tags, every CIGAR operation,
 * odd-length SEQ (which lands mid-byte in the 4-bit packing), '*' for SEQ and
 * QUAL, unmapped records, and integers at the boundaries where the tag encoder
 * switches width.
 *
 * Writes a BAM to the path given as argv[1]; the caller checks it with
 * samtools.  Run via tests/run_tests.sh.
 *
 *   c++ -std=c++17 -DWITH_ZLIB -iquote . -o unit_bam tests/unit_bam.cpp bam.cpp -lz
 */

#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>
#include "bam.h"

using namespace std;

static const char *kSam =
	"@HD\tVN:1.6\tSO:unsorted\n"
	"@SQ\tSN:chr1\tLN:100000\n"
	"@SQ\tSN:chr2\tLN:5000\n"
	"@RG\tID:rg1\tSM:sample1\n"
	// plain mapped read, tags at each integer width boundary
	"r1\t0\tchr1\t100\t42\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\t"
		"AS:i:-1\tXN:i:0\tNM:i:127\tZa:i:-128\tZb:i:255\tZc:i:256\tZd:i:-32768\t"
		"Ze:i:65535\tZf:i:65536\tZg:i:-2147483648\n"
	// odd-length SEQ: the last base occupies the high nibble of a final byte
	"r2\t0\tchr1\t200\t30\t7M\t*\t0\t0\tACGTACG\tIIIIIII\tXS:A:+\n"
	// every CIGAR operation.  H must be terminal and the query-consuming ops
	// (M I S = X) must sum to the SEQ length, or samtools rejects the record.
	"r3\t0\tchr1\t300\t20\t1H2S3M1I1P1D3N2=1X1S1H\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:2A5\n"
	// unmapped: refID -1, POS -1, and SEQ/QUAL still present
	"r4\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tYT:Z:UU\n"
	// SEQ and QUAL both absent
	"r5\t0\tchr2\t50\t10\t10M\t*\t0\t0\t*\t*\tNH:i:1\n"
	// QUAL absent but SEQ present
	"r6\t0\tchr2\t60\t10\t4M\t*\t0\t0\tACGT\t*\tNH:i:2\n"
	// float and array tags -- HISAT2 emits neither, so nothing else covers them
	"r7\t0\tchr2\t70\t10\t4M\t*\t0\t0\tACGT\tIIII\t"
		"Zh:f:3.5\tZi:B:c,-1,2,-3\tZj:B:i,100000,-100000\tZk:B:f,1.5,2.5\n"
	// paired, mate on another reference, negative TLEN
	"r8\t99\tchr1\t400\t60\t4M\t=\t500\t-104\tACGT\tIIII\tYS:i:-5\n"
	"r9\t147\tchr1\t500\tchr2_placeholder_removed\n";

int main(int argc, char **argv) {
	if(argc < 2) { cerr << "usage: unit_bam <out.bam>" << endl; return 2; }

	// Drop the deliberately malformed trailing line; it is only here as a
	// reminder that emitRecord() exits on a short record, which cannot be
	// tested in-process.
	string sam(kSam);
	size_t lastNl = sam.rfind("r9\t");
	sam.resize(lastNl);

	FILE *f = fopen(argv[1], "wb");
	if(f == NULL) { cerr << "could not open " << argv[1] << endl; return 2; }

	BamWriter bw;
	bw.init(f, true);
	// Feed the text in awkward chunk sizes so record boundaries land inside
	// chunks: writeSamText() must reassemble lines across calls.
	const size_t chunk = 7;
	for(size_t i = 0; i < sam.size(); i += chunk) {
		size_t n = sam.size() - i < chunk ? sam.size() - i : chunk;
		bw.writeSamText(sam.data() + i, n);
	}
	bw.close();
	cout << "wrote " << argv[1] << endl;
	return 0;
}
