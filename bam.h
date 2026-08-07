/*
 * BAM output for HISAT2.
 *
 * Two pieces: a BGZF block writer, and a converter that turns the SAM text
 * HISAT2 already produces into BAM records.
 *
 * Converting the text rather than encoding binary directly at the point of
 * alignment is deliberate. Records are formatted in exactly one place --
 * AlnSinkSam::appendMate plus SamConfig::printAlignedOptFlags, together some
 * 700 lines of field and optional-tag logic. A second binary encoder would
 * have to mirror all of it, and would silently fall out of step the first
 * time a tag was added to one and not the other. Going through the text keeps
 * a single source of truth, so BAM output is correct by construction and any
 * future tag works without further effort. The cost is re-parsing records we
 * just wrote, which is still far cheaper than the `| samtools view -b` this
 * replaces -- that pays the same parse plus a process boundary.
 */

#ifndef BAM_H_
#define BAM_H_

#ifdef WITH_ZLIB

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>
#include <zlib.h>

/**
 * Writes a BGZF stream: a series of gzip members, each carrying a BC extra
 * field giving the compressed block length, and each holding at most 64 KB of
 * uncompressed data.  Closing appends the 28-byte empty block that marks a
 * well-formed end of file; without it samtools reports the file as truncated.
 */
class BgzfWriter {
public:
	BgzfWriter() : out_(NULL), ownsFile_(false), level_(Z_DEFAULT_COMPRESSION),
	               head_(0), tail_(0), quit_(false) {}
	~BgzfWriter() { close(); }

	/**
	 * Take ownership of an already-open binary stream, or of stdout when
	 * ownsFile is false.
	 *
	 * nthreads > 1 compresses blocks on a pool of background threads.  This is
	 * the difference between BAM output being free and being the bottleneck:
	 * zlib manages about 38 MB/s at level 6, so a run producing 400 MB of SAM
	 * text spends ~11s compressing, against ~5s aligning.  BGZF blocks are
	 * independent by construction, so they parallelise exactly.
	 */
	void init(FILE *out, bool ownsFile, int level = Z_DEFAULT_COMPRESSION,
	          int nthreads = 1);

	/** Append uncompressed bytes; blocks are emitted as the buffer fills. */
	void write(const char *buf, size_t len);

	/** Flush the pending block, write the EOF marker, and close. */
	void close();

	bool isOpen() const { return out_ != NULL; }

private:
	// htslib uses this size too: it leaves room for the worst-case deflate
	// expansion plus the 18-byte header and 8-byte footer inside the 64 KB
	// a BGZF block is allowed to occupy.
	static const size_t BLOCK_SZ = 0xff00;

	/** Compress in_ and write it, on the calling thread. */
	void flushBlockSerial();

	// --- parallel path ---
	struct Block {
		std::vector<char> in;
		std::vector<char> out;   // full BGZF member, header and trailer included
		size_t            outLen;
		int               state; // 0 empty, 1 queued, 2 compressed
		Block() : outLen(0), state(0) {}
	};
	void   deflateBlock(Block& b);   // in -> out, no I/O; safe to run concurrently
	void   workerLoop();
	void   submitBlock();            // hand in_ to the ring
	void   drainOne();               // write the oldest finished block
	void   drainAll();
	bool   parallel() const { return !workers_.empty(); }

	FILE   *out_;
	bool    ownsFile_;
	int     level_;
	std::vector<char> in_;    // block being filled
	std::vector<char> out_buf_;

	std::vector<Block>       ring_;
	std::vector<std::thread> workers_;
	std::deque<size_t>       queued_;
	size_t                   head_;   // next slot the producer will fill
	size_t                   tail_;   // next slot to be written out
	bool                     quit_;
	std::mutex               mu_;
	std::condition_variable  cvWork_;   // workers wait for queued blocks
	std::condition_variable  cvDone_;   // producer waits for tail_ to finish
};

/**
 * Accepts SAM text in arbitrary chunks and writes BAM.
 *
 * Header lines are accumulated until the first alignment record arrives, at
 * which point the BAM header is emitted.  The reference list comes from the
 * @SQ lines in that text, so nothing extra has to be plumbed through from the
 * index.
 */
class BamWriter {
public:
	BamWriter() : headerDone_(false), pending_(), header_() {}
	~BamWriter() { close(); }

	void init(FILE *out, bool ownsFile, int level = Z_DEFAULT_COMPRESSION,
	          int nthreads = 1) {
		bgzf_.init(out, ownsFile, level, nthreads);
	}

	/** Feed SAM text.  Chunks need not align to record boundaries. */
	void writeSamText(const char *buf, size_t len);

	/** Emit any buffered record, finish the BGZF stream and close. */
	void close();

	bool isOpen() const { return bgzf_.isOpen(); }

private:
	void handleLine(const char *line, size_t len);
	void emitHeader();
	void emitRecord(const char *line, size_t len);

	BgzfWriter  bgzf_;
	bool        headerDone_;
	std::string pending_;     // partial line carried between chunks
	std::string header_;      // SAM header text seen so far
	std::vector<std::string> refNames_;
	std::vector<int32_t>     refLens_;
	// Keyed by views into refNames_, so a lookup costs no allocation.  The
	// converter runs once per alignment record, which on a real run is tens of
	// millions of times; a std::string temporary per lookup showed up clearly.
	std::unordered_map<std::string_view,int32_t> refId_;
	std::vector<char> rec_;   // scratch for the record being encoded
};

#endif /* WITH_ZLIB */
#endif /* BAM_H_ */
