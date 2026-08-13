/*
 * hifiasm_overlaps.h
 *
 * File-based bridge API for running hifiasm's initial candidate overlap
 * detection (+ base-level alignment/filter) from another C/C++ program such
 * as dinara.
 *
 * The function takes read files (FASTA/FASTQ, optionally gzipped) and an
 * output prefix, runs the same pipeline as the CLI, and writes the overlaps
 * to "<output_prefix>.ovlp.paf" (or "<output_prefix>.candidates.paf" when the
 * raw pre-alignment candidate set is requested). It returns the path of the
 * PAF that was written, which the caller then parses.
 *
 * PAF CIGARs use the extended =/X convention (op0='=' match, op1='X' mismatch,
 * op2='I' insertion, op3='D' deletion), exposed via the cg:Z: tag.
 */

#ifndef HIFIASM_OVERLAPS_H
#define HIFIASM_OVERLAPS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Options controlling overlap detection. Zero-initialize this struct
 * (e.g. `hifiasm_ovlp_opt_t o = {0};`) and then override what you need;
 * fields left at 0 fall back to the values documented below. */
typedef struct {
    int      threads;        /* number of worker threads (0 -> 1) */
    int      is_ont;         /* 1 = Nanopore/ONT preset, 0 = HiFi (default) */
    int      k_mer_length;   /* k-mer length (0 -> preset default) */
    int      mz_win;         /* minimizer window (0 -> preset default) */
    double   max_ov_diff;    /* max overlap error rate (0 -> preset default) */
    int      raw_candidates; /* 1 = emit pre-alignment candidates (.candidates.paf);
                                0 = emit alignment-filtered overlaps (.ovlp.paf) */
} hifiasm_ovlp_opt_t;

/*
 * Run overlap detection.
 *
 *   read_files    : array of input read file paths (FASTA/FASTQ[.gz])
 *   n_read_files  : number of entries in read_files (>= 1)
 *   output_prefix : output prefix; the PAF is written to
 *                   "<output_prefix>.ovlp.paf" (or ".candidates.paf")
 *   opt           : options (may be NULL for all-defaults, HiFi, 1 thread)
 *   out_paf_path  : if non-NULL, receives a malloc()'d string with the path of
 *                   the written PAF. Caller must free() it.
 *
 * Returns 0 on success, non-zero on error.
 *
 * NOTE: This routine uses process-global state (hifiasm's option and read
 * stores) and is NOT thread-safe; do not call it concurrently from multiple
 * threads. It may be called sequentially multiple times in one process.
 */
int hifiasm_detect_overlaps(const char *const *read_files,
                            int n_read_files,
                            const char *output_prefix,
                            const hifiasm_ovlp_opt_t *opt,
                            char **out_paf_path);

/*
 * One overlap returned by hifiasm_detect_overlaps_mem(). Fields mirror the PAF
 * columns the file-based path writes, so a caller can build the same records it
 * would have parsed from "<prefix>.ovlp.paf" without the text round-trip:
 *
 *   q_id, t_id       : read indices into the returned name table (see below).
 *                      q is the PAF query, t the PAF target.
 *   q_start, q_end   : half-open query interval [q_start, q_end) on the
 *                      forward strand of read q_id (PAF cols 3-4).
 *   t_start, t_end   : half-open target interval [t_start, t_end) on the
 *                      forward strand of read t_id (PAF cols 8-9). These are
 *                      forward-strand coordinates for both orientations, exactly
 *                      as the PAF path emits them.
 *   n_match          : residue matches (PAF col 10).
 *   block_len        : alignment block length (PAF col 11); the same value the
 *                      file path uses as the chain score / length filter key.
 *   is_same_strand   : 1 if the overlap is forward-forward (PAF strand '+'),
 *                      0 if reverse-complement (PAF strand '-').
 *   cigar_offset     : start index of this overlap's CIGAR in the token arena
 *                      returned via out_cigar (in uint16_t token units).
 *                      Meaningful only when cigar_len > 0.
 *   cigar_len        : number of CIGAR tokens for this overlap; 0 if none
 *                      (e.g. the raw pre-alignment candidate set carries no
 *                      CIGAR). The tokens are hifiasm's packed uint16_t form
 *                      (op << 14 | length), op 0='=' (match), 1='X' (mismatch),
 *                      2='I' (insertion, consumes query), 3='D' (deletion,
 *                      consumes target) -- the standard SAM/PAF convention.
 *   cigar_t_start    : the CIGAR's target anchor, in the SAME alignment frame
 *                      the tokens run in. The tokens align the forward query
 *                      (anchored at q_start) against the target in ALIGNMENT
 *                      orientation: for a same-strand overlap this equals
 *                      t_start; for a reverse-complement overlap the target is
 *                      reverse-complemented, so cigar_t_start is the offset from
 *                      the start of the reverse-complemented target (i.e.
 *                      t_len - t_end), NOT the forward t_start. A consumer that
 *                      reconstructs the target as the reverse complement for
 *                      '-' overlaps (as aligners typically do) can walk the
 *                      tokens directly from (q_start, cigar_t_start) with no
 *                      token reversal. Meaningful only when cigar_len > 0.
 */
typedef struct {
    uint32_t q_id;
    uint32_t t_id;
    uint32_t q_start;
    uint32_t q_end;
    uint32_t t_start;
    uint32_t t_end;
    uint32_t n_match;
    uint32_t block_len;
    uint8_t  is_same_strand;
    uint64_t cigar_offset;
    uint32_t cigar_len;
    uint32_t cigar_t_start;
} hifiasm_overlap_t;

/*
 * Run overlap detection and return the overlaps IN MEMORY instead of writing a
 * PAF file. This is the deep-integration path: it runs the same pipeline as
 * hifiasm_detect_overlaps() (candidate detection + base-level alignment/filter,
 * or the raw pre-alignment set when opt->raw_candidates is set), but hands the
 * surviving overlaps back as an array plus the read-name table, so the caller
 * resolves q_id/t_id to its own read ids by NAME (robust to load order) exactly
 * as it would when parsing the PAF.
 *
 *   read_files    : input read file paths (FASTA/FASTQ[.gz])
 *   n_read_files  : number of entries (>= 1)
 *   opt           : options (may be NULL for all-defaults, HiFi, 1 thread)
 *   out_ov        : receives a malloc()'d array of hifiasm_overlap_t.
 *                   *out_ov is NULL when there are zero overlaps.
 *   out_n_ov      : receives the number of overlaps in *out_ov.
 *   out_names     : receives a malloc()'d char buffer holding every read name
 *                   concatenated (NOT individually NUL-terminated).
 *   out_name_off  : receives a malloc()'d array of (*out_n_reads + 1) offsets
 *                   into *out_names; read i's name is the byte range
 *                   [out_name_off[i], out_name_off[i+1]). This mirrors
 *                   hifiasm's internal name_index layout.
 *   out_n_reads   : receives the number of reads (== hifiasm's total_reads).
 *   out_cigar     : receives a malloc()'d uint16_t arena holding every overlap's
 *                   CIGAR tokens concatenated. Each hifiasm_overlap_t references
 *                   its slice via (cigar_offset, cigar_len). May be NULL if the
 *                   caller does not want CIGARs, in which case cigar_len is 0 on
 *                   every overlap and *out_cigar is NULL. *out_cigar is also
 *                   NULL when there are zero CIGAR tokens overall.
 *   out_cigar_len : receives the total number of uint16_t tokens in *out_cigar.
 *                   May be NULL only if out_cigar is NULL.
 *
 * Ownership: on success the caller owns *out_ov, *out_names, *out_name_off and
 * *out_cigar and must release them with hifiasm_overlaps_mem_free(). Returns 0
 * on success, non-zero on error (in which case all out-params are set to NULL/0).
 *
 * NOTE: like hifiasm_detect_overlaps(), this uses hifiasm's process-global
 * option and read stores and is NOT thread-safe; do not call it concurrently
 * with other bridge entry points that touch those globals.
 */
int hifiasm_detect_overlaps_mem(const char *const *read_files,
                                int n_read_files,
                                const hifiasm_ovlp_opt_t *opt,
                                hifiasm_overlap_t **out_ov,
                                uint64_t *out_n_ov,
                                char **out_names,
                                uint64_t **out_name_off,
                                uint64_t *out_n_reads,
                                uint16_t **out_cigar,
                                uint64_t *out_cigar_len);

/*
 * Release the buffers returned by hifiasm_detect_overlaps_mem() /
 * hifiasm_detect_overlaps_from_store(). Any argument may be NULL (no-op for that
 * one). Safe to call with all-NULL.
 */
void hifiasm_overlaps_mem_free(hifiasm_overlap_t *ov,
                               char *names,
                               uint64_t *name_off,
                               uint16_t *cigar);


/* ---------------------------------------------------------------------------
 * Shared-read-store path (deepest integration): feed reads from memory.
 * ---------------------------------------------------------------------------
 * The functions above each read the input FASTA/FASTQ files themselves. When a
 * caller (dinara) has ALREADY loaded the reads, that re-reads the same bytes
 * from disk. The functions below instead take the reads in memory and load them
 * once into hifiasm's read store, so BOTH the marker filter build and overlap
 * detection can run without touching the input files again.
 *
 * The raw read bases are independent of k and homopolymer-compression, so the
 * SAME loaded store serves the no-HPC k=50 filter/sketch AND the HPC k=51
 * overlap detection; only the per-call k/w/HPC options differ.
 */

/* One read handed to the shared-store loader. seq must be raw ASCII bases
 * (A/C/G/T/N/...); non-ACGT bases are recorded as ambiguous, matching hifiasm's
 * file loader. name need not be NUL-terminated (name_len gives its length). */
typedef struct {
    const char *seq;
    uint64_t    seq_len;
    const char *name;
    uint32_t    name_len;
} hifiasm_read_t;

/*
 * Load `reads` into hifiasm's process-global read store (R_INF). Must be called
 * before hifiasm_build_filter_from_store() / hifiasm_detect_overlaps_from_store()
 * and released afterwards with hifiasm_reads_store_release(). Reads keep the
 * given order, so read index i corresponds to reads[i]. Returns 0 on success.
 *
 * NOTE: process-global and NOT thread-safe; one loaded store at a time.
 */
int hifiasm_reads_store_load(const hifiasm_read_t *reads, uint64_t n_reads);

/* Free the read store loaded by hifiasm_reads_store_load(). No-op if nothing is
 * loaded. */
void hifiasm_reads_store_release(void);

/*
 * Run overlap detection over the ALREADY-LOADED store, returning overlaps in
 * memory. Same outputs/ownership as hifiasm_detect_overlaps_mem(), but reads
 * nothing from disk. The name table is taken from the loaded store.
 */
int hifiasm_detect_overlaps_from_store(const hifiasm_ovlp_opt_t *opt,
                                       hifiasm_overlap_t **out_ov,
                                       uint64_t *out_n_ov,
                                       char **out_names,
                                       uint64_t **out_name_off,
                                       uint64_t *out_n_reads,
                                       uint16_t **out_cigar,
                                       uint64_t *out_cigar_len);

/* The store-fed filter builder (hifiasm_build_filter_from_store) is declared
 * further down, next to hifiasm_build_filter(), because it shares that
 * function's option and result types (hifiasm_filter_opt_t / hifiasm_filter_t).
 */

/*
 * One minimizer produced by hifiasm_sketch_minimizers().
 *
 *   pos : START position of the k-mer on the forward strand of the input
 *         sequence (0-based). hifiasm's internal sketch stores the END
 *         position; this bridge converts to START (end + 1 - span) so callers
 *         that key markers on the k-mer start (e.g. dinara) can use it
 *         directly. Guaranteed 0 <= pos and pos + span <= len.
 *   span: length of the k-mer in bases. With HPC disabled this always equals
 *         k; it is still returned so HPC-mode callers (span may exceed k) work.
 *   rev : 0 if the forward k-mer is canonical, 1 if the reverse complement is.
 *   hash: hifiasm's canonical yak hash of the k-mer (order-independent seed
 *         identity). Provided for callers that want hifiasm's identity; callers
 *         that need the actual k-mer bases should re-extract them from the
 *         sequence at [pos, pos+span) since the hash is not invertible.
 */
typedef struct {
    uint32_t pos;
    uint32_t span;
    uint32_t rev;
    uint64_t hash;
} hifiasm_minimizer_t;

/*
 * Sketch a single sequence into its minimizers using hifiasm's sketcher,
 * WITHOUT frequency filtering, position-table refinement, or subsampling
 * (hf/pt/sample_dist are all disabled), so the output depends only on
 * (seq, len, w, k, is_hpc).
 *
 *   seq     : sequence bases (ACGT/acgt; other chars are treated as ambiguous
 *             and break k-mers, matching hifiasm). Not required to be
 *             NUL-terminated; only the first `len` bytes are read.
 *   len     : number of bases in `seq` (> 0).
 *   w       : minimizer window (e.g. 50). Must satisfy 0 < w < 256.
 *   k       : k-mer length (e.g. 50). Must satisfy 0 < k <= 63.
 *   is_hpc  : 0 = no homopolymer compression (raw bases, span == k);
 *             1 = HPC (span may exceed k).
 *   out_mz  : receives a malloc()'d array of hifiasm_minimizer_t sorted by
 *             ascending START position. Caller must free() it. On zero
 *             minimizers, *out_mz is set to NULL.
 *   out_n   : receives the number of minimizers written to *out_mz.
 *
 * Returns 0 on success, non-zero on error (invalid arguments or allocation
 * failure). This function uses only local/heap state and does not touch the
 * process-global option or read stores, so unlike hifiasm_detect_overlaps it
 * is safe to call concurrently from multiple threads with distinct arguments.
 */
int hifiasm_sketch_minimizers(const char *seq,
                              int len,
                              int w,
                              int k,
                              int is_hpc,
                              hifiasm_minimizer_t **out_mz,
                              int *out_n);

/*
 * ---------------------------------------------------------------------------
 * Frequency-filter handle (for marker generation with overlap-path parity)
 * ---------------------------------------------------------------------------
 *
 * hifiasm's overlap path does not sketch reads raw: it first builds a
 * high-occurrence k-mer table ("hf") over all reads and, while sketching,
 * drops any minimizer whose k-mer occurs too frequently (repeats, low-signal
 * seeds). It then subsamples the survivors by genomic distance. A caller that
 * wants its markers to match hifiasm's overlap seeds must apply the same two
 * filters.
 *
 * This handle wraps that hf table. Build it once over the read set with
 * hifiasm_build_filter(), then pass it (together with a sample distance) to
 * hifiasm_sketch_minimizers_ctx_filtered() for each read. The table is a
 * standalone hash of k-mer values with no dependency on hifiasm's internal
 * read store, so it can outlive the build call and be shared read-only across
 * threads.
 *
 * IMPORTANT: the filter is k/w/HPC-specific. Build it with the SAME k, w and
 * is_hpc you will sketch with, or the frequency counts won't correspond to the
 * minimizers you filter. For dinara's no-HPC marker path that means
 * k=w=<marker k> and is_hpc=0.
 */
typedef struct hifiasm_filter_s hifiasm_filter_t;

/* Options for hifiasm_build_filter(). Zero-initialize and override as needed. */
typedef struct {
    int      threads;      /* worker threads for counting (0 -> 1)            */
    int      k_mer_length; /* k-mer length; MUST match the sketch k (0 -> 50) */
    int      mz_win;       /* minimizer window; MUST match sketch w (0 -> 50) */
    int      is_hpc;       /* 1 = HPC counting, 0 = no-HPC (dinara markers)   */
    int64_t  min_read_len; /* drop reads shorter than this; <0 keeps all.
                              Pass <0 so the count set matches every read you
                              intend to sketch (default: keep all).           */
} hifiasm_filter_opt_t;

/*
 * Build a high-occurrence k-mer filter table over the given read files.
 *
 *   read_files   : input read file paths (FASTA/FASTQ[.gz])
 *   n_read_files : number of entries (>= 1)
 *   opt          : options (NULL -> defaults: 1 thread, k=w=50, no-HPC,
 *                  keep all reads)
 *
 * Returns a filter handle on success, NULL on failure. Free it with
 * hifiasm_filter_destroy().
 *
 * NOTE: like hifiasm_detect_overlaps(), this uses hifiasm's process-global
 * option and read stores while it runs and is NOT thread-safe; do not call it
 * concurrently with other bridge entry points that touch those globals. The
 * RETURNED handle, however, is independent of those globals and is safe to use
 * concurrently from many sketch threads once built.
 */
hifiasm_filter_t *hifiasm_build_filter(const char *const *read_files,
                                       int n_read_files,
                                       const hifiasm_filter_opt_t *opt);

/* Destroy a filter handle. NULL is accepted (no-op). */
void hifiasm_filter_destroy(hifiasm_filter_t *hf);

/*
 * Build the filter over the ALREADY-LOADED store (see hifiasm_reads_store_load)
 * instead of reading files. Same options, return value and ownership as
 * hifiasm_build_filter(); reads nothing from disk. Free with
 * hifiasm_filter_destroy().
 */
hifiasm_filter_t *hifiasm_build_filter_from_store(const hifiasm_filter_opt_t *opt);

/*
 * A k-mer whose occurrence count is >= HIFIASM_FILTER_HIGH_OCC is considered
 * high-occurrence and is dropped as a minimizer seed by the filtered sketch.
 * This mirrors the internal threshold hifiasm's sketcher uses.
 */
#define HIFIASM_FILTER_HIGH_OCC (1 << 28)

/*
 * Look up the occurrence count recorded for a k-mer hash in the filter table.
 *
 *   hf   : filter handle (NULL -> returns 0)
 *   hash : the `hash` field of a hifiasm_minimizer_t (hifiasm's canonical yak
 *          hash, which is exactly the key this table is indexed by)
 *
 * Returns the recorded count, or a value >= HIFIASM_FILTER_HIGH_OCC for k-mers
 * marked high-occurrence (these are the seeds the filtered sketch drops), or 0
 * for k-mers not present in the table (rare k-mers below the counting floor).
 *
 * Mainly useful for inspection/testing: a minimizer emitted by
 * hifiasm_sketch_minimizers_ctx_filtered() with this same hf must always have
 * a count < HIFIASM_FILTER_HIGH_OCC.
 */
int32_t hifiasm_filter_count(const hifiasm_filter_t *hf, uint64_t hash);

/*
 * Test-only: return the raw underlying yak filter table (const void*) so tests
 * can call hifiasm's mz1_ha_sketch directly with the same hf the bridge uses,
 * and assert bit-for-bit parity. Not intended for production callers; NULL-safe.
 */
const void *hifiasm_filter_raw_for_test(const hifiasm_filter_t *hf);

/*
 * Reusable sketch context.
 *
 * hifiasm_sketch_minimizers() allocates and frees its work buffers on every
 * call, which is fine for one-shot use but wasteful when sketching many reads
 * (e.g. as an assembler's marker source). The context variant keeps the sketch
 * scratch buffers alive across calls, so steady-state sketching performs no
 * per-read heap allocation once the buffers have grown to fit.
 *
 * A context is NOT thread-safe: give each worker thread its own context. It
 * carries no read/option global state, so distinct contexts can be used
 * concurrently.
 */
typedef struct hifiasm_sketch_ctx_s hifiasm_sketch_ctx_t;

/* Create a context. Returns NULL on allocation failure. */
hifiasm_sketch_ctx_t *hifiasm_sketch_ctx_init(void);

/* Destroy a context and free its buffers. NULL is accepted (no-op). */
void hifiasm_sketch_ctx_destroy(hifiasm_sketch_ctx_t *ctx);

/*
 * Sketch a single sequence using context-owned storage. Semantics of seq/len/
 * w/k/is_hpc match hifiasm_sketch_minimizers().
 *
 * On success, *out_mz points into the context's internal buffer and *out_n is
 * the number of minimizers. The pointer is valid only until the next call on
 * the same context (or hifiasm_sketch_ctx_destroy); the CALLER MUST NOT free
 * it. On zero minimizers, *out_mz is NULL and *out_n is 0.
 *
 * The returned minimizers are sorted by ascending START position and
 * deduplicated (each START position appears once). Returns 0 on success,
 * non-zero on error (invalid arguments or allocation failure).
 */
int hifiasm_sketch_minimizers_ctx(hifiasm_sketch_ctx_t *ctx,
                                  const char *seq,
                                  int len,
                                  int w,
                                  int k,
                                  int is_hpc,
                                  const hifiasm_minimizer_t **out_mz,
                                  int *out_n);

/*
 * Like hifiasm_sketch_minimizers_ctx(), but applies hifiasm's overlap-path
 * minimizer filters so the output matches the seeds hifiasm uses for overlap
 * detection:
 *
 *   1. Frequency filter: drops any minimizer whose k-mer is marked
 *      high-occurrence in `hf`. Pass the handle from hifiasm_build_filter(),
 *      built with the SAME k/w/is_hpc used here. hf may be NULL to skip this
 *      filter (then this behaves like the unfiltered variant plus subsampling).
 *   2. Distance subsampling: when sample_dist > w, thins the surviving
 *      minimizers so roughly one high-quality seed is kept per sample_dist
 *      bases (hifiasm's select_mz_h). Pass sample_dist <= 0 (or <= w) to skip
 *      subsampling. hifiasm's overlap default is 500.
 *
 * Position-table refinement (hifiasm's optional third stage) is intentionally
 * NOT applied here: the overlap sketch itself runs with pt=NULL, so this
 * matches it. Semantics of seq/len/w/k/is_hpc and the returned pointer/lifetime
 * are identical to hifiasm_sketch_minimizers_ctx().
 *
 * Returns 0 on success, non-zero on error.
 */
int hifiasm_sketch_minimizers_ctx_filtered(hifiasm_sketch_ctx_t *ctx,
                                           const char *seq,
                                           int len,
                                           int w,
                                           int k,
                                           int is_hpc,
                                           const hifiasm_filter_t *hf,
                                           int sample_dist,
                                           const hifiasm_minimizer_t **out_mz,
                                           int *out_n);

/* -----------------------------------------------------------------------------
 * hifiasm's tuned overlap-path sketch parameters (single source of truth).
 *
 * These mirror the defaults hifiasm's own overlap detection uses
 * (CommandLines.cpp: k_mer_length=51, mz_win=51, mz_sample_dist=500,
 * mz_rewin=1000). A caller that wants its minimizer selection to match the
 * seeds hifiasm uses for overlaps should feed exactly these values to the
 * filtered sketch (k, w, sample_dist) and build the frequency filter at the
 * same k/w. Keep these in sync with hifiasm's defaults if they ever change.
 *
 * Note on subsampling: distance subsampling is applied by the sketch ONLY when
 * sample_dist > w. With these constants (500 > 51) it is active; a caller that
 * overrides sample_dist to <= w silently disables it.
 * -------------------------------------------------------------------------- */
#define HIFIASM_OVLP_K            51  /* overlap-path k-mer length            */
#define HIFIASM_OVLP_W            51  /* overlap-path minimizer window        */
#define HIFIASM_OVLP_SAMPLE_DIST 500 /* overlap-path distance subsampling    */
#define HIFIASM_OVLP_REWIN      1000 /* overlap-path select_mz_h window (ws)  */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HIFIASM_OVERLAPS_H */
