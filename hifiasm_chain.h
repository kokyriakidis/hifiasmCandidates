#ifndef HIFIASM_CHAIN_H
#define HIFIASM_CHAIN_H

/*
 * hifiasm_chain -- thin per-pair minimizer chainer.
 *
 * Exposes hifiasm's colinear-chaining DP (lchain_dp in Hash_Table.cpp) as a
 * single stateless call over a caller-supplied set of seed anchors. This is the
 * chaining core hifiasm runs after it collects minimizer hits between a read
 * pair; here it is decoupled from hifiasm's read stores and minimizer sketch so
 * a caller that already has its own markers can chain them with the SAME DP and
 * the SAME tuned parameters hifiasm uses for HiFi overlap detection.
 *
 * Intended use (dinara): dinara's default markers are hifiasm's non-HPC 51/51
 * minimizers, so for a candidate read pair dinara can (1) collect the markers
 * that fall inside the candidate overlap box on each read, (2) match them by
 * their shared k-mer identity into (query_pos, target_pos) anchors, (3) call
 * hifiasm_chain_pair() to get the colinear subset back, and (4) map the kept
 * anchor positions back to its own marker ordinals. No read sequences, no
 * re-sketching, and no name reconciliation cross this ABI -- only positions.
 *
 * Coordinate frame: all positions are RAW, forward, 0-based read positions at
 * the seed (marker) START, in the SAME strand frame for both reads. The caller
 * is responsible for orienting the pair before calling: for a same-strand
 * overlap feed forward positions on both reads; for a reverse-complement
 * overlap feed the target read's positions in its reverse-complement frame
 * (i.e. use the oriented read the caller will map ordinals against). The DP
 * itself is strand-agnostic -- it only requires that a colinear chain has
 * strictly increasing coordinates on both axes.
 */

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * One seed anchor fed to the chainer. q_pos/t_pos are the shared seed's start
 * position on the query and target read (raw, forward, oriented frame). id is
 * an opaque caller tag (e.g. a packed pair of marker ordinals) that the chainer
 * carries through untouched, so the caller can recover which of its markers a
 * kept anchor corresponds to without a second lookup.
 */
typedef struct {
    uint32_t q_pos;
    uint32_t t_pos;
    uint64_t id;
} hifiasm_chain_anchor_t;

/*
 * Tunables for the DP. Zero-initialize (or pass NULL to hifiasm_chain_pair) for
 * hifiasm's HiFi overlap defaults, which mirror set_lchain_dp_op(is_accurate=0,
 * mz_k=51) plus the HiFi band-width rate used in candidates.cpp:
 *
 *   max_skip  = 25      max_iter = 5000     max_dis = 5000
 *   bw_rate   = 0.02    (0.05 for ONT)
 *   chn_pen_gap/chn_pen_skip derived from mz_k via expf(-0.1*mz_k).
 *
 * seed_len is the per-anchor seed length used for the DP match score (the low
 * byte of hifiasm's k_mer_hit.cnt); it should equal the sketch k (51). A value
 * of 0 defaults to 51.
 */
typedef struct {
    int    seed_len;    /* per-anchor seed length (0 -> 51)                    */
    int    mz_k;        /* k used to derive gap/skip penalties (0 -> 51)       */
    int    is_ont;      /* 0 -> HiFi bw_rate 0.02; 1 -> ONT bw_rate 0.05       */
    int    max_skip;    /* 0 -> 25                                             */
    int    max_iter;    /* 0 -> 5000                                           */
    int    max_dis;     /* 0 -> 5000                                           */
    double bw_rate;     /* band-width rate; <= 0 -> is_ont default             */
} hifiasm_chain_opt_t;

/* Fill opt with the HiFi defaults described above. */
void hifiasm_chain_opt_init(hifiasm_chain_opt_t *opt);

/*
 * Chain one read pair.
 *
 *   anchors / n_anchor : input seed anchors. Need NOT be sorted; the chainer
 *                        sorts a private copy by query position (the input
 *                        contract of lchain_dp) and does not modify the caller
 *                        array. Anchors sharing the SAME (q_pos, t_pos) are
 *                        redundant seeds and should be de-duplicated by the
 *                        caller for best results, but are tolerated.
 *   q_len / t_len      : read lengths (bases) of query and target, used by the
 *                        band-width and gap penalties. Must be > 0.
 *   opt                : tunables (NULL -> hifiasm HiFi defaults).
 *   out_kept           : receives a malloc()'d array of the CHAINED anchors, a
 *                        colinear subset of the input, ordered by ascending
 *                        query position. Each element is the input anchor
 *                        (q_pos, t_pos, id preserved) that the DP kept. NULL
 *                        when the chain is empty. Caller frees with
 *                        hifiasm_chain_free().
 *   out_n_kept         : number of kept anchors (chain length).
 *   out_score          : DP chain score (hifiasm's shared_seed), or 0.
 *
 * Returns 0 on success (including the empty-chain case: out_n_kept == 0),
 * non-zero on invalid arguments or allocation failure (all out-params cleared).
 *
 * Thread-safety: uses only local/heap state (its own Chain_Data scratch per
 * call); safe to call concurrently with distinct arguments.
 */
int hifiasm_chain_pair(const hifiasm_chain_anchor_t *anchors,
                       uint64_t n_anchor,
                       uint32_t q_len,
                       uint32_t t_len,
                       const hifiasm_chain_opt_t *opt,
                       hifiasm_chain_anchor_t **out_kept,
                       uint64_t *out_n_kept,
                       int32_t *out_score);

/* Release an array returned via out_kept. NULL is tolerated. */
void hifiasm_chain_free(hifiasm_chain_anchor_t *kept);

#ifdef __cplusplus
}
#endif

#endif /* HIFIASM_CHAIN_H */
