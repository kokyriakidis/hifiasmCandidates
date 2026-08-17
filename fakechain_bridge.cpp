/*
 * fakechain_bridge.cpp -- see fakechain_bridge.h.
 *
 * Glue only: run the raw candidate overlap bridge, run the myloasm marker
 * index, reconcile reads by name, fake-chain each candidate, and gather the
 * chained anchors into one flat arena. All the marker/match/chain logic lives
 * in fakechain.cpp; all overlap detection lives in hifiasm_overlaps.cpp; this
 * file just wires them together and packages the results.
 */

#include "fakechain_bridge.h"

#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <unordered_map>

#include "myloasm_ffi.h"   /* myloasm_index_reads, MyloReadIndex */

void fakechain_bridge_opt_init(fakechain_bridge_opt_t *opt) {
    if (!opt) return;
    memset(opt, 0, sizeof(*opt));
    /* ovlp: zero-init means HiFi defaults; force raw candidates. */
    opt->ovlp.raw_candidates = 1;
    fakechain_opt_init(&opt->fakechain);
    opt->myloasm_k = 21;
    opt->myloasm_c = 11;
}

int fakechain_bridge_overlaps(const char *const *read_files,
                              int n_read_files,
                              const fakechain_bridge_opt_t *opt,
                              fakechain_overlap_t **out_ov,
                              uint64_t *out_n_ov,
                              fakechain_anchor_t **out_anchors,
                              uint64_t *out_n_anchors,
                              char **out_names,
                              uint64_t **out_name_off,
                              uint64_t *out_n_reads) {
    if (out_ov)        *out_ov = NULL;
    if (out_n_ov)      *out_n_ov = 0;
    if (out_anchors)   *out_anchors = NULL;
    if (out_n_anchors) *out_n_anchors = 0;
    if (out_names)     *out_names = NULL;
    if (out_name_off)  *out_name_off = NULL;
    if (out_n_reads)   *out_n_reads = 0;

    if (!read_files || n_read_files < 1 || !out_ov || !out_n_ov) return 1;

    fakechain_bridge_opt_t o;
    fakechain_bridge_opt_init(&o);
    if (opt) {
        o.ovlp = opt->ovlp;
        o.ovlp.raw_candidates = 1;              /* always raw for this path */
        o.fakechain = opt->fakechain;
        if (opt->myloasm_k) o.myloasm_k = opt->myloasm_k;
        if (opt->myloasm_c) o.myloasm_c = opt->myloasm_c;
    }

    /* --- Step 1: raw candidate overlaps + intervals + name table. --- */
    hifiasm_overlap_t *cand = NULL;
    uint64_t n_cand = 0;
    char *names = NULL;
    uint64_t *name_off = NULL;
    uint64_t n_reads = 0;
    int rc = hifiasm_detect_overlaps_mem(read_files, n_read_files, &o.ovlp,
                                         &cand, &n_cand,
                                         &names, &name_off, &n_reads,
                                         /*out_cigar=*/NULL, /*out_cigar_len=*/NULL);
    if (rc != 0) return 10 + rc;

    /* --- Step 2: myloasm per-read markers. --- */
    MyloReadIndex idx;
    memset(&idx, 0, sizeof(idx));
    int mrc = myloasm_index_reads(read_files, (size_t)n_read_files,
                                  o.myloasm_k, o.myloasm_c, o.ovlp.threads,
                                  &idx);
    if (mrc != 0) {
        hifiasm_overlaps_mem_free(cand, names, name_off, NULL);
        return 30 + mrc;
    }

    /* --- Map myloasm read base-id -> index in idx.reads. --- */
    std::unordered_map<std::string, size_t> mylo_by_name;
    mylo_by_name.reserve(idx.n_reads * 2);
    for (size_t i = 0; i < idx.n_reads; ++i) {
        const MyloReadMarkers &r = idx.reads[i];
        mylo_by_name.emplace(std::string(r.name, r.name_len), i);
    }

    /* Map hifiasm read id -> myloasm index (via the hifiasm name table). Base id
     * = first whitespace-delimited token of the hifiasm name, to match
     * myloasm's base_id. */
    std::vector<size_t> h2m(n_reads, (size_t)-1);
    for (uint64_t i = 0; i < n_reads; ++i) {
        const char *nm = names + name_off[i];
        size_t nlen = (size_t)(name_off[i + 1] - name_off[i]);
        size_t base = nlen;
        for (size_t j = 0; j < nlen; ++j) {
            char ch = nm[j];
            if (ch == ' ' || ch == '\t') { base = j; break; }
        }
        auto it = mylo_by_name.find(std::string(nm, base));
        if (it != mylo_by_name.end()) h2m[i] = it->second;
    }

    /* --- Step 3: fake-chain each candidate; gather anchors into one arena. --- */
    std::vector<fakechain_overlap_t> kept;
    kept.reserve(n_cand);
    std::vector<fakechain_anchor_t> anchor_arena;

    for (uint64_t i = 0; i < n_cand; ++i) {
        const hifiasm_overlap_t &c = cand[i];
        if (c.q_id >= n_reads || c.t_id >= n_reads) continue;
        size_t mq = h2m[c.q_id], mt = h2m[c.t_id];
        if (mq == (size_t)-1 || mt == (size_t)-1) continue; /* no markers */

        const MyloReadMarkers &Q = idx.reads[mq];
        const MyloReadMarkers &T = idx.reads[mt];

        fakechain_result_t fr;
        int frc = fakechain_pair(
            Q.syncmers, Q.n_syncmers, Q.snpmers, Q.n_snpmers,
            T.syncmers, T.n_syncmers, T.snpmers, T.n_snpmers,
            c.q_start, c.q_end, c.t_start, c.t_end,
            c.is_same_strand, &o.fakechain, &fr);
        if (frc != 0 || !fr.ok || fr.n_anchor == 0) {
            fakechain_result_free(&fr);
            continue;
        }

        fakechain_overlap_t out;
        memset(&out, 0, sizeof(out));
        out.q_id = c.q_id;
        out.t_id = c.t_id;
        out.q_start = fr.q_start;
        out.q_end   = fr.q_end;
        out.t_start = fr.t_start;
        out.t_end   = fr.t_end;
        out.is_same_strand = c.is_same_strand;
        out.chain_score = fr.score;
        out.n_anchor = fr.n_anchor;
        out.n_snpmer_anchor = fr.n_snpmer_anchor;
        out.anchor_offset = (uint64_t)anchor_arena.size();

        anchor_arena.insert(anchor_arena.end(),
                            fr.anchors, fr.anchors + fr.n_anchor);
        kept.push_back(out);

        fakechain_result_free(&fr);
    }

    myloasm_read_index_free(&idx);

    /* --- Package results. --- */
    fakechain_overlap_t *ov = NULL;
    if (!kept.empty()) {
        ov = (fakechain_overlap_t *)malloc(kept.size() * sizeof(fakechain_overlap_t));
        if (!ov) {
            hifiasm_overlaps_mem_free(cand, names, name_off, NULL);
            return 2;
        }
        memcpy(ov, kept.data(), kept.size() * sizeof(fakechain_overlap_t));
    }

    fakechain_anchor_t *anchors = NULL;
    if (!anchor_arena.empty()) {
        anchors = (fakechain_anchor_t *)malloc(
            anchor_arena.size() * sizeof(fakechain_anchor_t));
        if (!anchors) {
            free(ov);
            hifiasm_overlaps_mem_free(cand, names, name_off, NULL);
            return 3;
        }
        memcpy(anchors, anchor_arena.data(),
               anchor_arena.size() * sizeof(fakechain_anchor_t));
    }

    *out_ov = ov;
    *out_n_ov = kept.size();
    if (out_anchors)   *out_anchors = anchors;    else free(anchors);
    if (out_n_anchors) *out_n_anchors = anchor_arena.size();

    /* The candidate array is no longer needed; free ONLY it (names/name_off are
     * transferred to the caller unless unwanted). */
    free(cand);

    if (out_names)    *out_names = names;       else free(names);
    if (out_name_off) *out_name_off = name_off; else free(name_off);
    if (out_n_reads)  *out_n_reads = n_reads;

    return 0;
}

void fakechain_bridge_free(fakechain_overlap_t *ov,
                           fakechain_anchor_t *anchors,
                           char *names,
                           uint64_t *name_off) {
    /* ov and anchors are plain malloc; names/name_off came from the overlap
     * bridge but are plain malloc buffers too. Free each directly. */
    free(ov);
    free(anchors);
    hifiasm_overlaps_mem_free(NULL, names, name_off, NULL);
}
