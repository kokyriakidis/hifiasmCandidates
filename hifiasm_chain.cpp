/*
 * hifiasm_chain.cpp -- see hifiasm_chain.h.
 *
 * Wraps hifiasm's colinear-chaining DP (lchain_dp, Hash_Table.cpp) so a caller
 * with its own seed anchors can chain them with the exact DP and tuned
 * parameters hifiasm uses for HiFi overlap detection, without touching
 * hifiasm's read stores or minimizer sketch.
 *
 * Mapping onto lchain_dp's k_mer_hit contract:
 *   self_offset <- q_pos   (x / query axis; lchain "self")
 *   offset      <- t_pos   (y / target axis)
 *   cnt          = (1<<8) | seed_len
 *        low byte  -> per-anchor seed span (q_span in comput_sc_ch)
 *        high bits -> weight divisor in normal_w(sc, cnt>>8); 1 => sc/1 == sc
 *   readID       = input anchor index (0..n-1), used ONLY to recover the
 *                  caller's opaque `id` after chaining. lchain_dp never reads
 *                  readID for the DP itself (it is called per read-pair, all one
 *                  id) and copies whole k_mer_hit structs into des[], so the
 *                  index rides through untouched. 31 bits is ample for the seed
 *                  count of a single read pair.
 *
 * DP parameters mirror set_lchain_dp_op(is_accurate=0, mz_k) plus the HiFi
 * band-width rate 0.02 (0.05 ONT) used at the candidates.cpp overlap call.
 */

#include "hifiasm_chain.h"
#include "Hash_Table.h"   /* k_mer_hit, Chain_Data, overlap_region, lchain_dp */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>

/* Chain_Data lifecycle: defined in Hash_Table.cpp (C++ linkage) but not declared
 * in the public header. init == memset-zero, destory == free the six scratch
 * arrays. Declared here with matching C++ linkage. */
void init_Chain_Data(Chain_Data *x);
void destory_Chain_Data(Chain_Data *x);

/* set_lchain_dp_op(is_accurate=0, mz_k) reproduced locally: it is file-static in
 * anchor.cpp. div = 0.1 for the inaccurate (overlap) path; penalties scale by
 * expf(-div*mz_k). Values verified against anchor.cpp:2272. */
static void chain_dp_penalties(int mz_k, double *chn_pen_gap, double *chn_pen_skip) {
    const double div = 0.1;
    const double pen_gap = 0.5;
    const double pen_skip = 0.0005;
    const double tmp = expf(-div * (double)mz_k);
    *chn_pen_gap  = pen_gap * tmp;
    *chn_pen_skip = pen_skip * tmp;
}

void hifiasm_chain_opt_init(hifiasm_chain_opt_t *opt) {
    if (!opt) return;
    memset(opt, 0, sizeof(*opt));
    opt->seed_len = 51;
    opt->mz_k     = 51;
    opt->is_ont   = 0;
    opt->max_skip = 25;
    opt->max_iter = 5000;
    opt->max_dis  = 5000;
    opt->bw_rate  = 0.02;   /* HiFi (candidates.cpp:213); 0.05 for ONT */
}

int hifiasm_chain_pair(const hifiasm_chain_anchor_t *anchors,
                       uint64_t n_anchor,
                       uint32_t q_len,
                       uint32_t t_len,
                       const hifiasm_chain_opt_t *opt,
                       hifiasm_chain_anchor_t **out_kept,
                       uint64_t *out_n_kept,
                       int32_t *out_score) {
    if (out_kept)   *out_kept = NULL;
    if (out_n_kept) *out_n_kept = 0;
    if (out_score)  *out_score = 0;

    if (!anchors || !out_kept || !out_n_kept) return 1;
    if (q_len == 0 || t_len == 0) return 1;
    if (n_anchor == 0) return 0;                 /* empty chain, not an error */
    if (n_anchor > 0x7fffffffu) return 1;        /* index must fit readID:31  */

    hifiasm_chain_opt_t o;
    hifiasm_chain_opt_init(&o);
    if (opt) {
        if (opt->seed_len) o.seed_len = opt->seed_len;
        if (opt->mz_k)     o.mz_k     = opt->mz_k;
        o.is_ont = opt->is_ont;
        if (opt->max_skip) o.max_skip = opt->max_skip;
        if (opt->max_iter) o.max_iter = opt->max_iter;
        if (opt->max_dis)  o.max_dis  = opt->max_dis;
        o.bw_rate = (opt->bw_rate > 0.0) ? opt->bw_rate
                                         : (opt->is_ont ? 0.05 : 0.02);
    }
    const uint32_t seed_span = (o.seed_len > 0 && o.seed_len < 256)
                                   ? (uint32_t)o.seed_len : 51u;
    const uint32_t cnt = (1u << 8) | (seed_span & 0xffu);

    double chn_pen_gap, chn_pen_skip;
    chain_dp_penalties(o.mz_k, &chn_pen_gap, &chn_pen_skip);

    /* Build the k_mer_hit input array, sorted by (offset, self_offset) =
     * (t_pos, q_pos) as lchain_dp expects (sorted by offset). readID carries the
     * original anchor index so we can recover the caller's opaque id. */
    std::vector<k_mer_hit> a(n_anchor);
    for (uint64_t i = 0; i < n_anchor; ++i) {
        k_mer_hit &h = a[i];
        h.readID = (uint32_t)i;
        h.strand = 0;
        h.self_offset = anchors[i].q_pos;
        h.offset      = anchors[i].t_pos;
        h.cnt         = cnt;
    }
    std::sort(a.begin(), a.end(), [](const k_mer_hit &x, const k_mer_hit &y) {
        if (x.offset != y.offset) return x.offset < y.offset;
        return x.self_offset < y.self_offset;
    });

    /* Output/scratch for lchain_dp: des[] holds the kept chained anchors. */
    std::vector<k_mer_hit> des(n_anchor);
    Chain_Data dp;
    init_Chain_Data(&dp);
    overlap_region res;
    memset(&res, 0, sizeof(res));

    const uint64_t cL = lchain_dp(
        a.data(), (int64_t)n_anchor, des.data(), &dp, &res,
        (int64_t)o.max_skip, (int64_t)o.max_iter, (int64_t)o.max_dis,
        chn_pen_gap, chn_pen_skip, o.bw_rate,
        (int64_t)q_len, (int64_t)t_len, /*quick_check=*/0);

    destory_Chain_Data(&dp);

    if (cL == 0) {
        return 0;                                /* no chain */
    }

    hifiasm_chain_anchor_t *kept =
        (hifiasm_chain_anchor_t *)malloc((size_t)cL * sizeof(*kept));
    if (!kept) return 2;

    /* des[] is ordered chain-start -> chain-end (ascending on both axes). Each
     * element's readID indexes the ORIGINAL caller anchor array, from which we
     * copy the opaque id back out (des carries only positions). */
    for (uint64_t i = 0; i < cL; ++i) {
        const uint32_t src = des[i].readID;
        kept[i].q_pos = des[i].self_offset;
        kept[i].t_pos = des[i].offset;
        kept[i].id    = (src < n_anchor) ? anchors[src].id : 0;
    }

    *out_kept = kept;
    *out_n_kept = cL;
    if (out_score) *out_score = res.shared_seed;
    return 0;
}

void hifiasm_chain_free(hifiasm_chain_anchor_t *kept) {
    free(kept);
}
