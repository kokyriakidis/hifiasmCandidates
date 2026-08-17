/*
 * fakechain.cpp -- see fakechain.h for the design rationale.
 *
 * Given two reads' myloasm marker lists (open syncmers + SNPmers, each with a
 * canonical k=21 key) and a candidate overlap interval, this:
 *   1. clips keys to an even marker_k (default 20),
 *   2. matches markers between the pair by clipped-key equality, keeping ALL
 *      m:n matches (repeats preserved) in SEPARATE key spaces per type,
 *   3. merges syncmer + SNPmer matches into one tagged anchor set,
 *   4. chains them with myloasm's own DP (myloasm_chain FFI), and
 *   5. emits the ordered chained anchors, tagged by type, for dinara.
 *
 * There is no base-level alignment and no hand-rolled chaining here: myloasm's
 * dp_anchors_v2 does the colinear chaining, keeping behaviour identical to
 * upstream myloasm.
 */

#include "fakechain.h"

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <unordered_map>

void fakechain_opt_init(fakechain_opt_t *opt) {
    if (!opt) return;
    opt->marker_k = 20;
    opt->gap_cost = 11;
    opt->match_score = 1;
    opt->min_chain_len = 3;
    opt->band_mult = 20;
    opt->min_anchors = 1;
}

void fakechain_result_free(fakechain_result_t *r) {
    if (!r) return;
    free(r->anchors);
    r->anchors = NULL;
    r->n_anchor = 0;
}

namespace {

/* Clip a canonical k=21 key to the low 2*marker_k bits (keeps the middle base
 * at bits [20,21] for marker_k>=11, so SNPmer alleles stay distinct). Because
 * the input key is canonical, this deterministic function yields identical
 * clipped keys for the same locus on both reads regardless of strand. */
static inline uint64_t clip_key(uint64_t key, int marker_k) {
    if (marker_k >= 32) return key;               /* no clip */
    const uint64_t mask = (marker_k <= 0)
        ? 0ull
        : ((1ull << (2 * marker_k)) - 1ull);
    return key & mask;
}

/* What each merged anchor's `slot` recovers after chaining: its clipped key and
 * type tag. Query/target positions come back from the chain output directly. */
struct SlotInfo {
    uint64_t key;
    uint8_t  tag;
};

/* Match one marker type between the pair (m:n, clipped-key equality), appending
 * a MyloAnchor per match and a parallel SlotInfo. Only markers whose pos lies in
 * the respective interval are considered. Updates max_mult with the largest
 * target-side key multiplicity seen for a matched key (myloasm's convention). */
static void match_type(const MyloMarker *q_mk, size_t q_n,
                       const MyloMarker *t_mk, size_t t_n,
                       uint32_t q_s, uint32_t q_e,
                       uint32_t t_s, uint32_t t_e,
                       int marker_k, uint8_t tag, uint8_t is_same_strand,
                       std::vector<MyloAnchor> &anchors,
                       std::vector<SlotInfo> &slots,
                       size_t &max_mult) {
    /* Index the target side: clipped key -> list of target positions in [t_s,t_e). */
    std::unordered_map<uint64_t, std::vector<uint32_t> > tindex;
    tindex.reserve(t_n * 2 + 1);
    for (size_t j = 0; j < t_n; ++j) {
        uint32_t p = t_mk[j].pos;
        if (p < t_s || p >= t_e) continue;
        tindex[clip_key(t_mk[j].key, marker_k)].push_back(p);
    }
    if (tindex.empty()) return;

    const uint32_t rel_strand = is_same_strand ? 0u : 1u;

    for (size_t i = 0; i < q_n; ++i) {
        uint32_t qp = q_mk[i].pos;
        if (qp < q_s || qp >= q_e) continue;
        uint64_t k = clip_key(q_mk[i].key, marker_k);
        auto it = tindex.find(k);
        if (it == tindex.end()) continue;
        const std::vector<uint32_t> &tps = it->second;
        if (tps.size() > max_mult) max_mult = tps.size();
        for (size_t m = 0; m < tps.size(); ++m) {
            MyloAnchor a;
            a.pos1 = (rel_strand << 31) | qp;   /* strand in bit 31 */
            a.pos2 = tps[m];                    /* raw target pos   */
            a.slot = (uint32_t)slots.size();
            a._pad = 0;
            anchors.push_back(a);
            SlotInfo si;
            si.key = k;
            si.tag = tag;
            slots.push_back(si);
        }
    }
}

} // namespace

int fakechain_pair(const MyloMarker *q_sync, size_t q_n_sync,
                   const MyloMarker *q_snp,  size_t q_n_snp,
                   const MyloMarker *t_sync, size_t t_n_sync,
                   const MyloMarker *t_snp,  size_t t_n_snp,
                   uint32_t q_s, uint32_t q_e,
                   uint32_t t_s, uint32_t t_e,
                   uint8_t is_same_strand,
                   const fakechain_opt_t *opt,
                   fakechain_result_t *out) {
    if (!out) return 1;
    memset(out, 0, sizeof(*out));
    out->is_same_strand = is_same_strand;
    out->q_start = q_s; out->q_end = q_e;
    out->t_start = t_s; out->t_end = t_e;

    if (q_e <= q_s || t_e <= t_s) return 2;

    fakechain_opt_t o;
    fakechain_opt_init(&o);
    if (opt) {
        if (opt->marker_k)      o.marker_k = opt->marker_k;
        if (opt->gap_cost)      o.gap_cost = opt->gap_cost;
        if (opt->match_score)   o.match_score = opt->match_score;
        if (opt->min_chain_len) o.min_chain_len = opt->min_chain_len;
        if (opt->band_mult)     o.band_mult = opt->band_mult;
        if (opt->min_anchors)   o.min_anchors = opt->min_anchors;
    }
    if (o.marker_k <= 0 || (o.marker_k & 1)) return 3;   /* must be positive, even */

    /* --- Match markers between the pair (m:n), merged + tagged. --- */
    std::vector<MyloAnchor> anchors;
    std::vector<SlotInfo>   slots;
    size_t max_mult = 0;

    match_type(q_sync, q_n_sync, t_sync, t_n_sync, q_s, q_e, t_s, t_e,
               o.marker_k, (uint8_t)FAKECHAIN_SYNCMER, is_same_strand,
               anchors, slots, max_mult);
    match_type(q_snp, q_n_snp, t_snp, t_n_snp, q_s, q_e, t_s, t_e,
               o.marker_k, (uint8_t)FAKECHAIN_SNPMER, is_same_strand,
               anchors, slots, max_mult);

    out->n_match  = (uint32_t)anchors.size();
    out->max_mult = (uint32_t)max_mult;

    if (anchors.empty()) return 0;                 /* ok stays 0, no chain */

    /* --- Chain with myloasm's DP. band = band_mult * max_mult. --- */
    size_t band = (size_t)o.band_mult * max_mult;  /* max_mult >= 1 here */
    std::vector<MyloChainAnchor> chained(anchors.size());
    size_t out_n = 0;
    int chain_score = 0;
    int is_reverse = 0;

    int crc = myloasm_chain(anchors.data(), anchors.size(),
                            o.gap_cost, o.match_score, band, o.min_chain_len,
                            chained.data(), chained.size(),
                            &out_n, &chain_score, &is_reverse);
    if (crc != 0) return 100 + crc;

    out->score = chain_score;

    if (out_n == 0) return 0;                       /* no chain of >= min length */

    /* --- Emit ordered, tagged chained anchors. --- */
    fakechain_anchor_t *arr =
        (fakechain_anchor_t *)malloc(out_n * sizeof(fakechain_anchor_t));
    if (!arr) return 4;

    uint32_t n_snp = 0;
    uint32_t q_first = chained[0].qpos, q_last = chained[0].qpos;
    uint32_t t_min = chained[0].tpos, t_max = chained[0].tpos;
    for (size_t i = 0; i < out_n; ++i) {
        const MyloChainAnchor &ca = chained[i];
        const SlotInfo &si = slots[ca.slot];       /* slot is our own index */
        arr[i].q_pos = ca.qpos;                     /* strand already stripped */
        arr[i].t_pos = ca.tpos;
        arr[i].key   = si.key;
        arr[i].tag   = si.tag;
        memset(arr[i]._pad, 0, sizeof(arr[i]._pad));
        if (si.tag == FAKECHAIN_SNPMER) n_snp++;
        if (ca.qpos < q_first) q_first = ca.qpos;
        if (ca.qpos > q_last)  q_last  = ca.qpos;
        if (ca.tpos < t_min)   t_min   = ca.tpos;
        if (ca.tpos > t_max)   t_max   = ca.tpos;
    }

    out->anchors = arr;
    out->n_anchor = (uint32_t)out_n;
    out->n_snpmer_anchor = n_snp;

    if ((int)out_n < o.min_anchors) {
        out->ok = 0;
        return 0;
    }

    /* Bounding forward intervals from the chain endpoints, extended by one
     * marker_k to cover the last marker, clamped to the candidate interval. */
    uint32_t qa = q_first;
    uint32_t qb = q_last + (uint32_t)o.marker_k;
    if (qb > q_e) qb = q_e;
    if (qa < q_s) qa = q_s;
    uint32_t ta = t_min;
    uint32_t tb = t_max + (uint32_t)o.marker_k;
    if (tb > t_e) tb = t_e;
    if (ta < t_s) ta = t_s;

    out->q_start = qa; out->q_end = qb;
    out->t_start = ta; out->t_end = tb;
    out->ok = 1;
    return 0;
}
