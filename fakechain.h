#ifndef FAKECHAIN_H
#define FAKECHAIN_H

/*
 * "Fake-chain" refinement of candidate overlaps using myloasm markers.
 *
 * Pipeline context:
 *   1. hifiasm produces candidate overlaps + their query/target INTERVALS with
 *      its default HPC k=51/w=51 minimizers (good for finding pairs).
 *   2. The base-level alignment/filter is SKIPPED (raw candidate set).
 *   3. myloasm indexes every read (myloasm_index_reads): open syncmers + SNPmers
 *      at k=21, each carrying the read's CANONICAL k-mer key at that position.
 *   4. THIS module refines each candidate by MATCHING those markers between the
 *      pair and CHAINING them with myloasm's own DP (myloasm_chain, exposed via
 *      the FFI). Markers are CLIPPED to the hifiasm-derived interval on BOTH
 *      reads before matching (query pos in [q_s,q_e), target pos in [t_s,t_e)),
 *      so the chain runs only over the in-interval anchor set, not the full
 *      reads. This is safe because hifiasm candidate coordinates and myloasm
 *      marker positions are both RAW read coordinates (same space), and the
 *      interval already brackets the colinear region hifiasm found. The output
 *      is the ORDERED chained anchor list, tagged by marker type, handed to
 *      dinara for downstream scoring/phasing.
 *
 * Why myloasm markers and not hifiasm's here: myloasm uses OPEN SYNCMERS, which
 * (unlike minimizers) are context-independent and give unbiased sequence-
 * identity estimates (myloasm paper section 4.3). SNPmers add phase-informative
 * anchors at heterozygous sites. Syncmers and SNPmers are chained TOGETHER in a
 * single DP pass.
 *
 * Why myloasm's chainer and not hifiasm's: hifiasm's chain_DP is calibrated for
 * its long HPC k=51 anchors (min_score = k) and uses a repeat-defence term that
 * needs a per-anchor multiplicity count the marker FFI does not populate. For
 * the dense, uniform, short (k=20/21) raw markers used here, myloasm's
 * dp_anchors_v2 (match_score=1, gap_cost=c, band=max_mult*20) is the correct
 * model and keeps repeats (m:n matches) instead of forcing 1:1.
 *
 * Key handling: myloasm keys are already CANONICAL (min of forward /
 * reverse-complement), hence strand-invariant -- the same locus yields the same
 * key on both reads for same-strand AND reverse-complement overlaps. Matching is
 * therefore pure key equality; no bit-order or RC handling is needed here. The
 * per-pair relative strand comes from hifiasm (is_same_strand) and is encoded
 * into bit 31 of each anchor's query position for the DP.
 *
 * Even-length clip for dinara: keys are clipped from k=21 to EVEN marker_k
 * (default 20) by masking the low 2*marker_k bits. With k=21 the middle (het)
 * base sits at bits [20,21]; marker_k=20 keeps it, so a SNPmer's two alleles
 * stay distinct after the clip. Syncmer and SNPmer keys are matched in SEPARATE
 * key spaces (a syncmer never matches a SNPmer), so no cross-type collision.
 *
 * Repeats: a clipped key may occur multiple times on either read. We keep ALL
 * m:n matches (repeat multiplicity preserved); myloasm's DP resolves them via
 * band = max_mult * 20. Nothing is forced to 1:1.
 *
 * fakechain_pair is PURE apart from calling the (stateless) myloasm_chain FFI:
 * it takes two reads' marker lists + a candidate interval + relative strand and
 * returns the chained, ordered, tagged anchors. No global state, so it is
 * unit-testable and safe to call concurrently on disjoint outputs.
 */

#include <stdint.h>
#include <stddef.h>

#include "myloasm_ffi.h"   /* MyloMarker, MyloAnchor, MyloChainAnchor */

#ifdef __cplusplus
extern "C" {
#endif

/* Anchor type tags. Carried out of the chain so dinara can weight/phase by
 * type. Packed into the low bits of each output anchor's `tag`. */
enum {
    FAKECHAIN_SYNCMER = 0,
    FAKECHAIN_SNPMER  = 1
};

typedef struct {
    int marker_k;      /* even clip length, <= 21 (0 -> 20)                    */
    int gap_cost;      /* DP gap cost c (0 -> myloasm default 11)             */
    int match_score;   /* DP per-anchor score (0 -> myloasm default 1)        */
    int min_chain_len; /* DP minimum chain length (0 -> myloasm default 3)    */
    int band_mult;     /* band = band_mult * max_mult; if 0 -> 20 (myloasm)   */
    int min_anchors;   /* reject chain if fewer kept anchors (0 -> 1)         */
} fakechain_opt_t;

/* Fill `opt` with defaults (marker_k=20, gap_cost=11, match_score=1,
 * min_chain_len=3, band_mult=20, min_anchors=1). */
void fakechain_opt_init(fakechain_opt_t *opt);

/* One chained anchor handed to dinara. Positions are RAW read positions
 * (forward, 0-based) at the marker start; strand is uniform per overlap and
 * reported once in the result. `key` is the clipped (marker_k) canonical key,
 * shared by both reads at this locus. `tag` is FAKECHAIN_SYNCMER/SNPMER. */
typedef struct {
    uint32_t q_pos;
    uint32_t t_pos;
    uint64_t key;
    uint8_t  tag;
    uint8_t  _pad[7];
} fakechain_anchor_t;

typedef struct {
    /* Bounding forward, half-open intervals derived from the first/last chained
     * anchor (target min/max, so it covers both strands). Equal to the input
     * candidate interval when no chain is found. */
    uint32_t q_start, q_end;
    uint32_t t_start, t_end;
    uint8_t  is_same_strand;

    int32_t  score;          /* myloasm chain score (DP), 0 if no chain       */
    uint32_t n_anchor;       /* chained anchors in `anchors` (ordered by q)   */
    uint32_t n_snpmer_anchor;/* of those, SNPmer anchors                       */
    uint32_t n_match;        /* total m:n matches fed to the DP (pre-chain)   */
    uint32_t max_mult;       /* max key multiplicity across the two reads      */
    uint8_t  ok;             /* 1 if n_anchor >= min_anchors                   */

    /* malloc()'d array of `n_anchor` chained anchors, ordered by query pos.
     * NULL when n_anchor == 0. Free with fakechain_result_free(). */
    fakechain_anchor_t *anchors;
} fakechain_result_t;

/* Release the anchor array held by `r` (idempotent; leaves r->anchors NULL). */
void fakechain_result_free(fakechain_result_t *r);

/*
 * Refine ONE candidate overlap: match myloasm markers between the pair (m:n,
 * canonical key equality, syncmers + SNPmers in separate key spaces), merge the
 * matches into one tagged anchor set, and chain them with myloasm's DP.
 *
 *   q_sync,q_n_sync : query read's syncmer markers (from MyloReadMarkers)
 *   q_snp, q_n_snp  : query read's SNPmer markers
 *   t_sync,t_n_sync : target read's syncmer markers
 *   t_snp, t_n_snp  : target read's SNPmer markers
 *                     Each MyloMarker is {pos (forward, 0-based), key
 *                     (canonical k=21)}. Lists need not be sorted.
 *   q_s,q_e         : candidate query interval [q_s,q_e) (forward). Used only
 *                     to restrict query-side markers fed to the matcher.
 *   t_s,t_e         : candidate target interval [t_s,t_e) (forward). Used only
 *                     to restrict target-side markers fed to the matcher.
 *   is_same_strand  : 1 forward/forward, 0 reverse-complement (from hifiasm).
 *   opt             : options (NULL -> defaults).
 *   out             : receives the refined overlap + chained anchors. The
 *                     caller must fakechain_result_free(out).
 *
 * Returns 0 on success (including "no chain", out->ok==0 and out->anchors NULL),
 * non-zero on error (invalid arguments or a chainer/allocation failure).
 */
int fakechain_pair(const MyloMarker *q_sync, size_t q_n_sync,
                   const MyloMarker *q_snp,  size_t q_n_snp,
                   const MyloMarker *t_sync, size_t t_n_sync,
                   const MyloMarker *t_snp,  size_t t_n_snp,
                   uint32_t q_s, uint32_t q_e,
                   uint32_t t_s, uint32_t t_e,
                   uint8_t is_same_strand,
                   const fakechain_opt_t *opt,
                   fakechain_result_t *out);

#ifdef __cplusplus
}
#endif

#endif /* FAKECHAIN_H */
