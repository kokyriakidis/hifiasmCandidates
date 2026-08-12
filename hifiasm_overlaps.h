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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HIFIASM_OVERLAPS_H */
