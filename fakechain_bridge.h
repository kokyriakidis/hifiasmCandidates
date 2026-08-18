#ifndef FAKECHAIN_BRIDGE_H
#define FAKECHAIN_BRIDGE_H

/*
 * End-to-end fake-chain bridge.
 *
 * Composes the whole two-pass pipeline behind one call:
 *   1. hifiasm raw candidate overlaps + intervals (HPC k=51), NO base-level
 *      alignment/filter (hifiasm_detect_overlaps_mem with raw_candidates=1).
 *   2. myloasm per-read markers: open syncmers + SNPmers at k=21
 *      (myloasm_index_reads).
 *   3. For each candidate pair, fakechain_pair() matches markers between the
 *      pair (m:n, canonical key equality) and CHAINS them with myloasm's DP,
 *      producing an ORDERED, TAGGED chained-anchor array per overlap
 *      (see fakechain.h). Dinara consumes these anchors and does all downstream
 *      scoring/phasing; this bridge does NOT score.
 *
 * Read identities are reconciled by NAME: hifiasm returns a name table for its
 * q_id/t_id; myloasm returns each read's base id. This is robust to load order
 * differences between the two engines.
 *
 * Output shape: one fakechain_overlap_t per surviving candidate (those that
 * produced a chain of >= min_anchors). Each references its chained anchors as a
 * slice [anchor_offset, anchor_offset + n_anchor) of a single flat anchor arena,
 * so the whole result is three malloc()'d buffers (overlaps, anchors, name
 * table) freed together with fakechain_bridge_free().
 */

#include <stdint.h>
#include <stddef.h>

#include "hifiasm_overlaps.h"  /* hifiasm_ovlp_opt_t */
#include "fakechain.h"         /* fakechain_opt_t, fakechain_anchor_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    hifiasm_ovlp_opt_t ovlp;      /* candidate-detection options. raw_candidates
                                     is FORCED to 1 by this bridge.            */
    fakechain_opt_t    fakechain; /* fake-chain / chaining options             */
    int                myloasm_k; /* SNPmer/syncmer k (0 -> 21)                */
    int                myloasm_c; /* syncmer compression (0 -> 11)             */
    int                min_overlap_len; /* drop a chained overlap unless BOTH its
                                     query and target spans are >= this many
                                     RAW read bases. 0 -> keep all.           */
} fakechain_bridge_opt_t;

/* Fill defaults: HiFi candidates, raw_candidates=1, fakechain defaults,
 * myloasm k=21/c=11. threads is left 0 (caller may set opt->ovlp.threads and
 * it is also used for the myloasm pass). */
void fakechain_bridge_opt_init(fakechain_bridge_opt_t *opt);

/*
 * One refined overlap: a read pair, its per-pair relative strand, the candidate
 * interval it came from, and a slice of the shared anchor arena holding its
 * ordered chained anchors. No score is emitted here -- dinara scores from the
 * anchors.
 *
 *   q_id, t_id       : read indices into the returned name table.
 *   q_start/q_end    : bounding forward query interval of the chain.
 *   t_start/t_end    : bounding forward target interval of the chain.
 *   is_same_strand   : 1 forward/forward, 0 reverse-complement (from hifiasm).
 *   chain_score      : myloasm DP chain score (informational; dinara may ignore).
 *   n_anchor         : number of chained anchors for this overlap.
 *   n_snpmer_anchor  : of those, SNPmer anchors.
 *   anchor_offset    : start index of this overlap's anchors in the flat
 *                      anchor arena returned via out_anchors.
 */
typedef struct {
    uint32_t q_id;
    uint32_t t_id;
    uint32_t q_start, q_end;
    uint32_t t_start, t_end;
    uint8_t  is_same_strand;
    int32_t  chain_score;
    uint32_t n_anchor;
    uint32_t n_snpmer_anchor;
    uint64_t anchor_offset;
} fakechain_overlap_t;

/*
 * Run the full pipeline over `read_files` and return refined overlaps + their
 * chained anchors in memory.
 *
 *   read_files/n_read_files : input FASTA/FASTQ[.gz] paths (>= 1).
 *   opt                     : options (NULL -> defaults).
 *   out_ov / out_n_ov       : malloc()'d array of fakechain_overlap_t and its
 *                             length. *out_ov is NULL when zero overlaps.
 *   out_anchors / out_n_anchors : malloc()'d flat arena of fakechain_anchor_t
 *                             (all overlaps' anchors concatenated, each overlap's
 *                             slice given by anchor_offset/n_anchor) and its
 *                             total length. *out_anchors is NULL when zero
 *                             anchors. Either pointer may be NULL if unwanted.
 *   out_names/out_name_off/out_n_reads : hifiasm's read-name table (same layout
 *                             as hifiasm_detect_overlaps_mem), so q_id/t_id
 *                             resolve to names. Any may be NULL if not needed.
 *
 * Ownership: caller frees *out_ov, *out_anchors, *out_names, *out_name_off with
 * fakechain_bridge_free(). Returns 0 on success, non-zero on error (all
 * out-params set to NULL/0).
 *
 * NOTE: uses hifiasm's process-global stores (via the overlap bridge) and the
 * myloasm FFI; NOT thread-safe. One call at a time per process.
 */
int fakechain_bridge_overlaps(const char *const *read_files,
                              int n_read_files,
                              const fakechain_bridge_opt_t *opt,
                              fakechain_overlap_t **out_ov,
                              uint64_t *out_n_ov,
                              fakechain_anchor_t **out_anchors,
                              uint64_t *out_n_anchors,
                              char **out_names,
                              uint64_t **out_name_off,
                              uint64_t *out_n_reads);

/* Release buffers from fakechain_bridge_overlaps(). Any argument may be NULL. */
void fakechain_bridge_free(fakechain_overlap_t *ov,
                           fakechain_anchor_t *anchors,
                           char *names,
                           uint64_t *name_off);

#ifdef __cplusplus
}
#endif

#endif /* FAKECHAIN_BRIDGE_H */
