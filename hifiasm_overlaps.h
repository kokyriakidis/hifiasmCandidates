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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HIFIASM_OVERLAPS_H */
