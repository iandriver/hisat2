/*
 * BAM output for HISAT2.  See bam.h for why this converts SAM text rather
 * than encoding records directly.
 */

#ifdef WITH_ZLIB

#include "bam.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

// -----------------------------------------------------------------------
// little-endian putters.  BAM is little-endian regardless of host.
// -----------------------------------------------------------------------

static inline void put32(vector<char>& v, uint32_t x) {
	v.push_back((char)(x & 0xff));
	v.push_back((char)((x >> 8) & 0xff));
	v.push_back((char)((x >> 16) & 0xff));
	v.push_back((char)((x >> 24) & 0xff));
}


// -----------------------------------------------------------------------
// Hot-path helpers.  emitRecord() runs once per alignment record -- tens of
// millions of times on a real run -- so it parses out of the line in place and
// appends through a raw pointer.  The first version built a std::string per
// numeric field and push_back'd the record a byte at a time, which cost about
// 2s per million reads on its own.
// -----------------------------------------------------------------------

static inline int64_t parseI64(const char *p, size_t len) {
	int64_t v = 0;
	size_t i = 0;
	bool neg = false;
	if(i < len && (p[i] == '-' || p[i] == '+')) { neg = (p[i] == '-'); i++; }
	for(; i < len; i++) {
		if(p[i] < '0' || p[i] > '9') break;
		v = v * 10 + (p[i] - '0');
	}
	return neg ? -v : v;
}

/** Appends into a caller-sized buffer; no bounds checks in the inner loop. */
struct Buf {
	char *p;
	char *base;
	explicit Buf(char *b) : p(b), base(b) {}
	size_t size() const { return (size_t)(p - base); }
	void u8(uint8_t x)  { *p++ = (char)x; }
	void u16(uint16_t x) { *p++ = (char)(x & 0xff); *p++ = (char)((x >> 8) & 0xff); }
	void u32(uint32_t x) {
		*p++ = (char)(x & 0xff);        *p++ = (char)((x >> 8) & 0xff);
		*p++ = (char)((x >> 16) & 0xff); *p++ = (char)((x >> 24) & 0xff);
	}
	void mem(const void *s, size_t n) { memcpy(p, s, n); p += n; }
};

// -----------------------------------------------------------------------
// BgzfWriter
// -----------------------------------------------------------------------

void BgzfWriter::init(FILE *out, bool ownsFile, int level, int nthreads) {
	out_ = out;
	ownsFile_ = ownsFile;
	level_ = level;
	in_.reserve(BLOCK_SZ);
	if(nthreads > 1 && level != 0) {
		// Enough slots that workers stay fed while the producer waits on the
		// oldest one: blocks must be *written* in order even though they are
		// compressed out of order.  Each slot is 64 KB, so this is ~1-3 MB.
		size_t nslots = (size_t)nthreads * 2 + 2;
		ring_.resize(nslots);
		for(int i = 0; i < nthreads; i++) {
			workers_.push_back(std::thread(&BgzfWriter::workerLoop, this));
		}
	}
}

void BgzfWriter::write(const char *buf, size_t len) {
	if(out_ == NULL) return;
	while(len > 0) {
		size_t room = BLOCK_SZ - in_.size();
		size_t take = len < room ? len : room;
		in_.insert(in_.end(), buf, buf + take);
		buf += take;
		len -= take;
		if(in_.size() == BLOCK_SZ) {
			if(parallel()) submitBlock(); else flushBlockSerial();
		}
	}
}

/**
 * Compress b.in into b.out as a complete BGZF member.  Touches no shared
 * state, which is what lets several of these run at once.
 */
void BgzfWriter::deflateBlock(Block& b) {
	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	// Negative windowBits selects a raw deflate stream: BGZF supplies its own
	// gzip header and trailer.
	if(deflateInit2(&zs, level_, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
		cerr << "Error: could not initialize BAM compression" << endl;
		exit(1);
	}
	uLong bound = deflateBound(&zs, (uLong)b.in.size());
	if(b.out.size() < (size_t)bound + 64) b.out.resize((size_t)bound + 64);

	zs.next_in   = (Bytef*)&b.in[0];
	zs.avail_in  = (uInt)b.in.size();
	zs.next_out  = (Bytef*)&b.out[0] + 18;
	zs.avail_out = (uInt)(b.out.size() - 26);
	int r = deflate(&zs, Z_FINISH);
	if(r != Z_STREAM_END) {
		cerr << "Error: BAM compression failed (zlib code " << r << ")" << endl;
		exit(1);
	}
	size_t clen = (size_t)zs.total_out;
	deflateEnd(&zs);

	size_t bsize = 18 + clen + 8;
	if(bsize > 65536) {
		// Cannot happen with BLOCK_SZ == 0xff00, but a silent overlong block
		// would produce a file no BAM reader could parse.
		cerr << "Error: BGZF block overflow (" << bsize << " bytes)" << endl;
		exit(1);
	}

	unsigned char *h = (unsigned char*)&b.out[0];
	static const unsigned char kHdr[18] = {
		0x1f, 0x8b, 0x08, 0x04,           // magic, deflate, FEXTRA
		0, 0, 0, 0,                        // MTIME
		0, 0xff,                           // XFL, OS (unknown)
		6, 0,                              // XLEN
		'B', 'C', 2, 0, 0, 0               // SI1 SI2 SLEN BSIZE(lo,hi)
	};
	memcpy(h, kHdr, 18);
	uint16_t bs = (uint16_t)(bsize - 1);
	h[16] = (unsigned char)(bs & 0xff);
	h[17] = (unsigned char)((bs >> 8) & 0xff);

	uLong crc = crc32(0L, Z_NULL, 0);
	crc = crc32(crc, (const Bytef*)&b.in[0], (uInt)b.in.size());
	unsigned char *t = h + 18 + clen;
	uint32_t isize = (uint32_t)b.in.size();
	t[0] = (unsigned char)(crc & 0xff);
	t[1] = (unsigned char)((crc >> 8) & 0xff);
	t[2] = (unsigned char)((crc >> 16) & 0xff);
	t[3] = (unsigned char)((crc >> 24) & 0xff);
	t[4] = (unsigned char)(isize & 0xff);
	t[5] = (unsigned char)((isize >> 8) & 0xff);
	t[6] = (unsigned char)((isize >> 16) & 0xff);
	t[7] = (unsigned char)((isize >> 24) & 0xff);
	b.outLen = bsize;
}

void BgzfWriter::flushBlockSerial() {
	if(in_.empty()) return;
	Block b;
	b.in.swap(in_);
	b.out.swap(out_buf_);
	deflateBlock(b);
	if(fwrite(&b.out[0], 1, b.outLen, out_) != b.outLen) {
		cerr << "Error: could not write BAM output" << endl;
		exit(1);
	}
	b.in.clear();
	in_.swap(b.in);
	out_buf_.swap(b.out);
}

void BgzfWriter::workerLoop() {
	for(;;) {
		size_t idx;
		{
			std::unique_lock<std::mutex> lk(mu_);
			while(queued_.empty() && !quit_) cvWork_.wait(lk);
			if(queued_.empty() && quit_) return;
			idx = queued_.front();
			queued_.pop_front();
		}
		deflateBlock(ring_[idx]);
		{
			std::lock_guard<std::mutex> lk(mu_);
			ring_[idx].state = 2;
		}
		cvDone_.notify_all();
	}
}

/** Hand the filled block to the ring, waiting for a free slot if need be. */
void BgzfWriter::submitBlock() {
	if(in_.empty()) return;
	// Wait for the slot we are about to reuse, writing out finished blocks in
	// order until it frees up.  This is also the backpressure that stops the
	// aligner running arbitrarily far ahead of compression.
	while(ring_[head_].state != 0) drainOne();
	ring_[head_].in.swap(in_);
	in_.clear();
	in_.reserve(BLOCK_SZ);
	{
		std::lock_guard<std::mutex> lk(mu_);
		ring_[head_].state = 1;
		queued_.push_back(head_);
	}
	cvWork_.notify_one();
	head_ = (head_ + 1) % ring_.size();
}

/** Write the oldest outstanding block, blocking until it has been compressed. */
void BgzfWriter::drainOne() {
	Block& b = ring_[tail_];
	{
		std::unique_lock<std::mutex> lk(mu_);
		if(b.state == 0) return;             // nothing outstanding here
		while(b.state != 2) cvDone_.wait(lk);
	}
	if(fwrite(&b.out[0], 1, b.outLen, out_) != b.outLen) {
		cerr << "Error: could not write BAM output" << endl;
		exit(1);
	}
	b.in.clear();
	{
		std::lock_guard<std::mutex> lk(mu_);
		b.state = 0;
	}
	tail_ = (tail_ + 1) % ring_.size();
}

void BgzfWriter::drainAll() {
	for(size_t i = 0; i < ring_.size(); i++) {
		if(ring_[tail_].state != 0) drainOne();
		else tail_ = (tail_ + 1) % ring_.size();
	}
}

void BgzfWriter::close() {
	if(out_ == NULL) return;
	if(parallel()) {
		if(!in_.empty()) submitBlock();
		drainAll();
		{
			std::lock_guard<std::mutex> lk(mu_);
			quit_ = true;
		}
		cvWork_.notify_all();
		for(size_t i = 0; i < workers_.size(); i++) workers_[i].join();
		workers_.clear();
	} else {
		flushBlockSerial();
	}
	// The canonical empty BGZF block.  samtools treats a BAM without it as
	// truncated, which is exactly the check worth passing.
	static const unsigned char eof[28] = {
		0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
	};
	fwrite(eof, 1, 28, out_);
	fflush(out_);
	if(ownsFile_) fclose(out_);
	out_ = NULL;
}

// -----------------------------------------------------------------------
// BamWriter
// -----------------------------------------------------------------------

void BamWriter::writeSamText(const char *buf, size_t len) {
	size_t start = 0;
	for(size_t i = 0; i < len; i++) {
		if(buf[i] != '\n') continue;
		if(pending_.empty()) {
			handleLine(buf + start, i - start);
		} else {
			pending_.append(buf + start, i - start);
			handleLine(pending_.data(), pending_.size());
			pending_.clear();
		}
		start = i + 1;
	}
	if(start < len) pending_.append(buf + start, len - start);
}

void BamWriter::handleLine(const char *line, size_t len) {
	// Strip a trailing CR so CRLF input does not end up inside the last field.
	if(len > 0 && line[len-1] == '\r') len--;
	if(len == 0) return;
	if(!headerDone_ && line[0] == '@') {
		header_.append(line, len);
		header_.push_back('\n');
		return;
	}
	if(!headerDone_) emitHeader();
	emitRecord(line, len);
}

void BamWriter::emitHeader() {
	headerDone_ = true;

	// Reference list comes from the @SQ lines.  A BAM whose header lacks them
	// cannot express refID, so refuse rather than write a file that looks fine
	// until something tries to resolve a reference name.
	size_t pos = 0;
	while(pos < header_.size()) {
		size_t eol = header_.find('\n', pos);
		if(eol == string::npos) eol = header_.size();
		if(header_.compare(pos, 3, "@SQ") == 0) {
			string name;
			int32_t len = -1;
			size_t f = pos;
			while(f < eol) {
				size_t fe = header_.find('\t', f);
				if(fe == string::npos || fe > eol) fe = eol;
				if(header_.compare(f, 3, "SN:") == 0) {
					name = header_.substr(f + 3, fe - f - 3);
				} else if(header_.compare(f, 3, "LN:") == 0) {
					len = (int32_t)strtol(header_.substr(f + 3, fe - f - 3).c_str(), NULL, 10);
				}
				f = fe + 1;
			}
			if(!name.empty() && len >= 0) {
				refNames_.push_back(name);
				refLens_.push_back(len);
			}
		}
		pos = eol + 1;
	}
	// Index the names only once refNames_ has stopped growing: the map keys are
	// views into those strings, and a push_back that reallocated the vector
	// would leave every key dangling.
	refId_.reserve(refNames_.size());
	for(size_t i = 0; i < refNames_.size(); i++) {
		refId_[std::string_view(refNames_[i])] = (int32_t)i;
	}
	if(refNames_.empty()) {
		cerr << "Error: BAM output needs @SQ header lines, but none were written." << endl
		     << "Remove --no-hd/--no-sq, or write SAM instead." << endl;
		exit(1);
	}

	vector<char> h;
	h.push_back('B'); h.push_back('A'); h.push_back('M'); h.push_back('\1');
	put32(h, (uint32_t)header_.size());
	h.insert(h.end(), header_.begin(), header_.end());
	put32(h, (uint32_t)refNames_.size());
	for(size_t i = 0; i < refNames_.size(); i++) {
		put32(h, (uint32_t)(refNames_[i].size() + 1));
		h.insert(h.end(), refNames_[i].begin(), refNames_[i].end());
		h.push_back('\0');
		put32(h, (uint32_t)refLens_[i]);
	}
	bgzf_.write(&h[0], h.size());
	header_.clear();
}

/** BAI binning scheme, from the SAM specification. */
static int reg2bin(int64_t beg, int64_t end) {
	end--;
	if(beg >> 14 == end >> 14) return (int)(((1 << 15) - 1) / 7 + (beg >> 14));
	if(beg >> 17 == end >> 17) return (int)(((1 << 12) - 1) / 7 + (beg >> 17));
	if(beg >> 20 == end >> 20) return (int)(((1 << 9)  - 1) / 7 + (beg >> 20));
	if(beg >> 23 == end >> 23) return (int)(((1 << 6)  - 1) / 7 + (beg >> 23));
	if(beg >> 26 == end >> 26) return (int)(((1 << 3)  - 1) / 7 + (beg >> 26));
	return 0;
}

static inline uint8_t seqNibble(char c) {
	switch(c) {
		case '=': return 0;  case 'A': case 'a': return 1;
		case 'C': case 'c': return 2;  case 'M': case 'm': return 3;
		case 'G': case 'g': return 4;  case 'R': case 'r': return 5;
		case 'S': case 's': return 6;  case 'V': case 'v': return 7;
		case 'T': case 't': return 8;  case 'W': case 'w': return 9;
		case 'Y': case 'y': return 10; case 'H': case 'h': return 11;
		case 'K': case 'k': return 12; case 'D': case 'd': return 13;
		case 'B': case 'b': return 14; default:  return 15; // N and anything else
	}
}

static const char *kCigarOps = "MIDNSHP=X";

/** Append an optional field, given the SAM "TAG:TYPE:VALUE" text. */
static void putTag(Buf& r, const char *f, size_t flen) {
	if(flen < 5 || f[2] != ':' || f[4] != ':') return; // not a well-formed tag
	const char type = f[3];
	const char *val = f + 5;
	const size_t vlen = flen - 5;
	char *mark = r.p;
	r.u8((uint8_t)f[0]);
	r.u8((uint8_t)f[1]);
	switch(type) {
		case 'A':
			r.u8('A');
			r.u8(vlen > 0 ? (uint8_t)val[0] : 0);
			break;
		case 'i': {
			// Narrow to the smallest type that fits, as htslib does.  Purely a
			// size win; every reader widens these back to the same integer.
			int64_t v = parseI64(val, vlen);
			if(v >= -128 && v <= 127)          { r.u8('c'); r.u8((uint8_t)(int8_t)v); }
			else if(v >= 0 && v <= 255)        { r.u8('C'); r.u8((uint8_t)v); }
			else if(v >= -32768 && v <= 32767) { r.u8('s'); r.u16((uint16_t)(int16_t)v); }
			else if(v >= 0 && v <= 65535)      { r.u8('S'); r.u16((uint16_t)v); }
			else                               { r.u8('i'); r.u32((uint32_t)(int32_t)v); }
			break;
		}
		case 'f': {
			char tmp[64];
			size_t n = vlen < sizeof(tmp) - 1 ? vlen : sizeof(tmp) - 1;
			memcpy(tmp, val, n); tmp[n] = '\0';
			float fv = (float)atof(tmp);
			uint32_t bits;
			memcpy(&bits, &fv, 4);
			r.u8('f');
			r.u32(bits);
			break;
		}
		case 'Z':
		case 'H':
			r.u8((uint8_t)type);
			r.mem(val, vlen);
			r.u8(0);
			break;
		case 'B': {
			// B:subtype,v1,v2,...
			if(vlen < 1) { r.p = mark; return; }
			const char sub = val[0];
			r.u8('B');
			r.u8((uint8_t)sub);
			char *countAt = r.p;
			r.u32(0);
			uint32_t count = 0;
			size_t i = 1;
			while(i < vlen) {
				if(val[i] != ',') { i++; continue; }
				size_t b = i + 1, e = b;
				while(e < vlen && val[e] != ',') e++;
				if(sub == 'f') {
					char tmp[64];
					size_t n = (e - b) < sizeof(tmp) - 1 ? (e - b) : sizeof(tmp) - 1;
					memcpy(tmp, val + b, n); tmp[n] = '\0';
					float fv = (float)atof(tmp);
					uint32_t bits; memcpy(&bits, &fv, 4);
					r.u32(bits);
				} else {
					int64_t v = parseI64(val + b, e - b);
					switch(sub) {
						case 'c': r.u8((uint8_t)(int8_t)v); break;
						case 'C': r.u8((uint8_t)v); break;
						case 's': r.u16((uint16_t)(int16_t)v); break;
						case 'S': r.u16((uint16_t)v); break;
						case 'i': r.u32((uint32_t)(int32_t)v); break;
						case 'I': r.u32((uint32_t)v); break;
						default: r.p = mark; return; // unknown subtype: drop the tag
					}
				}
				count++;
				i = e;
			}
			countAt[0] = (char)(count & 0xff);
			countAt[1] = (char)((count >> 8) & 0xff);
			countAt[2] = (char)((count >> 16) & 0xff);
			countAt[3] = (char)((count >> 24) & 0xff);
			break;
		}
		default:
			// Unknown type: drop the tag rather than write something a reader
			// would misparse and run off the end of the record with.
			r.p = mark;
			break;
	}
}

void BamWriter::emitRecord(const char *line, size_t len) {
	// Split into fields.  Field 11 onward are optional tags.
	const char *f[11];
	size_t fl[11];
	int nf = 0;
	const char *tagStart = NULL;
	size_t start = 0;
	for(size_t i = 0; i <= len && nf < 11; i++) {
		if(i == len || line[i] == '\t') {
			f[nf] = line + start;
			fl[nf] = i - start;
			nf++;
			start = i + 1;
			if(nf == 11) tagStart = (start <= len) ? line + start : NULL;
		}
	}
	if(nf < 11) {
		cerr << "Error: malformed SAM record while writing BAM (" << nf
		     << " fields)" << endl;
		exit(1);
	}

	const int     flag  = (int)parseI64(f[1], fl[1]);
	const int64_t pos   = parseI64(f[3], fl[3]) - 1;
	const int     mapq  = (int)parseI64(f[4], fl[4]);
	const int64_t pnext = parseI64(f[7], fl[7]) - 1;
	const int64_t tlen  = parseI64(f[8], fl[8]);

	int32_t refID = -1;
	if(!(fl[2] == 1 && f[2][0] == '*')) {
		unordered_map<string_view,int32_t>::const_iterator it =
			refId_.find(string_view(f[2], fl[2]));
		if(it == refId_.end()) {
			cerr << "Error: reference \"" << string(f[2], fl[2])
			     << "\" is not in the BAM header" << endl;
			exit(1);
		}
		refID = it->second;
	}
	int32_t nextRefID = -1;
	if(fl[6] == 1 && f[6][0] == '=') {
		nextRefID = refID;
	} else if(!(fl[6] == 1 && f[6][0] == '*')) {
		unordered_map<string_view,int32_t>::const_iterator it =
			refId_.find(string_view(f[6], fl[6]));
		if(it != refId_.end()) nextRefID = it->second;
	}

	const bool seqStar  = (fl[9] == 1 && f[9][0] == '*');
	const int32_t lseq  = seqStar ? 0 : (int32_t)fl[9];
	const bool qualStar = (fl[10] == 1 && f[10][0] == '*');

	// One conservative sizing up front so the appender needs no bounds checks.
	// Every part of a BAM record is smaller than its SAM text except the fixed
	// 36-byte core and the NUL terminators, so the line length plus a fixed
	// margin is always enough.
	if(rec_.size() < len + 256) rec_.resize(len + 256);
	Buf r(&rec_[0]);

	r.u32(0);                    // block_size, patched below
	r.u32((uint32_t)refID);
	r.u32((uint32_t)pos);
	r.u8((uint8_t)(fl[0] + 1));
	r.u8((uint8_t)mapq);
	char *binAt = r.p;
	r.u16(0);                    // bin, patched once the CIGAR span is known
	char *ncigAt = r.p;
	r.u16(0);                    // n_cigar_op, patched below
	r.u16((uint16_t)flag);
	r.u32((uint32_t)lseq);
	r.u32((uint32_t)nextRefID);
	r.u32((uint32_t)pnext);
	r.u32((uint32_t)tlen);
	r.mem(f[0], fl[0]);
	r.u8(0);

	// CIGAR
	uint32_t nops = 0;
	int64_t refSpan = 0;
	if(!(fl[5] == 1 && f[5][0] == '*')) {
		uint32_t n = 0;
		for(size_t i = 0; i < fl[5]; i++) {
			const char c = f[5][i];
			if(c >= '0' && c <= '9') {
				n = n * 10 + (uint32_t)(c - '0');
			} else {
				const char *q = strchr(kCigarOps, c);
				if(q == NULL || c == '\0') {
					cerr << "Error: unrecognized CIGAR operation '" << c << "'" << endl;
					exit(1);
				}
				const uint32_t opi = (uint32_t)(q - kCigarOps);
				r.u32((n << 4) | opi);
				nops++;
				// M, D, N, = and X consume reference
				if(opi == 0 || opi == 2 || opi == 3 || opi == 7 || opi == 8) refSpan += n;
				n = 0;
			}
		}
	}
	ncigAt[0] = (char)(nops & 0xff);
	ncigAt[1] = (char)((nops >> 8) & 0xff);
	const uint16_t bin = (uint16_t)(refID < 0 || pos < 0 ? 4680
	                                : reg2bin(pos, pos + (refSpan > 0 ? refSpan : 1)));
	binAt[0] = (char)(bin & 0xff);
	binAt[1] = (char)((bin >> 8) & 0xff);

	// SEQ, two bases per byte, high nibble first
	if(lseq > 0) {
		const char *s = f[9];
		for(int32_t i = 0; i < lseq; i += 2) {
			const uint8_t hi = seqNibble(s[i]);
			const uint8_t lo = (i + 1 < lseq) ? seqNibble(s[i+1]) : 0;
			r.u8((uint8_t)((hi << 4) | lo));
		}
		// QUAL: 0xff throughout means "unavailable", which is how '*' is stored
		if(qualStar) {
			memset(r.p, 0xff, (size_t)lseq);
			r.p += lseq;
		} else {
			const char *q = f[10];
			const int32_t n = (int32_t)fl[10] < lseq ? (int32_t)fl[10] : lseq;
			for(int32_t i = 0; i < n; i++) r.u8((uint8_t)(q[i] - 33));
			if(n < lseq) { memset(r.p, 0xff, (size_t)(lseq - n)); r.p += lseq - n; }
		}
	}

	// Optional fields
	if(tagStart != NULL) {
		const char *end = line + len;
		const char *p = tagStart;
		while(p < end) {
			const char *q = p;
			while(q < end && *q != '\t') q++;
			putTag(r, p, (size_t)(q - p));
			p = q + 1;
		}
	}

	const uint32_t blockSize = (uint32_t)(r.size() - 4);
	rec_[0] = (char)(blockSize & 0xff);
	rec_[1] = (char)((blockSize >> 8) & 0xff);
	rec_[2] = (char)((blockSize >> 16) & 0xff);
	rec_[3] = (char)((blockSize >> 24) & 0xff);
	bgzf_.write(&rec_[0], r.size());
}

void BamWriter::close() {
	if(!bgzf_.isOpen()) return;
	if(!pending_.empty()) {
		handleLine(pending_.data(), pending_.size());
		pending_.clear();
	}
	// A run that aligned nothing still owes a valid header.
	if(!headerDone_ && !header_.empty()) emitHeader();
	bgzf_.close();
}

#endif /* WITH_ZLIB */
