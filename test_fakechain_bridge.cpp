/*
 * End-to-end smoke test for fakechain_bridge_overlaps().
 *
 * Runs the full two-pass pipeline (hifiasm raw candidates -> myloasm markers ->
 * per-pair match + myloasm chain) over a real HiFi read file and asserts the
 * OUTPUT INVARIANTS that dinara relies on, without hard-coding exact counts
 * (which depend on the data/engine internals):
 *
 *   - the call succeeds and returns a name table;
 *   - every overlap's anchor slice [anchor_offset, +n_anchor) lies inside the
 *     returned anchor arena and slices are contiguous / non-overlapping in emit
 *     order;
 *   - q_id/t_id are valid read indices and differ;
 *   - each overlap has n_anchor >= 1 and n_snpmer_anchor <= n_anchor;
 *   - anchors are ordered by query position and carry a valid type tag;
 *   - target positions are monotone in the strand-appropriate direction
 *     (non-decreasing for same-strand, non-increasing for RC), matching the
 *     colinear chain the DP produced;
 *   - bounding coords bracket the chained anchors.
 *
 * If the pipeline finds zero overlaps on this input the test still passes (it
 * only asserts the empty-result contract), so it never gives a false failure on
 * data changes; it exists to catch crashes, leaks-of-contract, and broken
 * slice/ordering invariants across the FFI boundary.
 *
 * Build/run:  make test_fakechain_bridge && ./test_fakechain_bridge <reads>
 */
#include <cstdio>
#include <cstdint>
#include <cstring>

#include "fakechain_bridge.h"

static int g_fail = 0;
#define CHECK(c,m) do{ if(!(c)){ std::fprintf(stderr,"FAIL: %s\n",(m)); ++g_fail; } }while(0)

int main(int argc, char **argv) {
    const char *reads = (argc > 1)
        ? argv[1]
        : "GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq";
    const char *files[1] = { reads };

    fakechain_bridge_opt_t opt;
    fakechain_bridge_opt_init(&opt);
    opt.ovlp.threads = 4;

    fakechain_overlap_t *ov = NULL;
    fakechain_anchor_t  *anchors = NULL;
    uint64_t n_ov = 0, n_anchors = 0, n_reads = 0;
    char *names = NULL;
    uint64_t *name_off = NULL;

    int rc = fakechain_bridge_overlaps(files, 1, &opt,
                                       &ov, &n_ov,
                                       &anchors, &n_anchors,
                                       &names, &name_off, &n_reads);
    CHECK(rc == 0, "bridge returns success");
    if (rc != 0) {
        std::fprintf(stderr, "bridge failed rc=%d\n", rc);
        return 1;
    }
    CHECK(n_reads > 0, "read name table populated");
    CHECK(names != NULL && name_off != NULL, "name buffers returned");

    std::fprintf(stderr,
        "bridge: reads=%llu overlaps=%llu anchors=%llu\n",
        (unsigned long long)n_reads, (unsigned long long)n_ov,
        (unsigned long long)n_anchors);

    /* Empty-result contract. */
    if (n_ov == 0) {
        CHECK(ov == NULL || n_ov == 0, "zero overlaps -> consistent");
        CHECK(n_anchors == 0, "zero overlaps -> zero anchors");
    }

    uint64_t running = 0;                 /* contiguous slice cursor */
    uint64_t total_anchor_seen = 0;
    for (uint64_t i = 0; i < n_ov; ++i) {
        const fakechain_overlap_t &o = ov[i];

        CHECK(o.q_id < n_reads && o.t_id < n_reads, "overlap read ids in range");
        CHECK(o.q_id != o.t_id, "overlap is between two distinct reads");
        CHECK(o.n_anchor >= 1, "overlap has >= 1 anchor");
        CHECK(o.n_snpmer_anchor <= o.n_anchor, "snpmer anchors <= total");
        CHECK(o.is_same_strand == 0 || o.is_same_strand == 1, "strand flag 0/1");

        /* Slice inside the arena, contiguous in emit order. */
        CHECK(o.anchor_offset == running, "anchor slices contiguous in emit order");
        CHECK(o.anchor_offset + o.n_anchor <= n_anchors, "anchor slice within arena");
        running += o.n_anchor;
        total_anchor_seen += o.n_anchor;

        const fakechain_anchor_t *a = anchors + o.anchor_offset;
        uint32_t qmin = a[0].q_pos, qmax = a[0].q_pos;
        uint32_t tmin = a[0].t_pos, tmax = a[0].t_pos;
        bool q_ordered = true, t_monotone = true, tags_ok = true;
        for (uint32_t k = 0; k < o.n_anchor; ++k) {
            if (a[k].tag != FAKECHAIN_SYNCMER && a[k].tag != FAKECHAIN_SNPMER)
                tags_ok = false;
            if (k) {
                if (a[k].q_pos < a[k-1].q_pos) q_ordered = false;   /* non-decreasing */
                if (o.is_same_strand) {
                    if (a[k].t_pos < a[k-1].t_pos) t_monotone = false;
                } else {
                    if (a[k].t_pos > a[k-1].t_pos) t_monotone = false;
                }
            }
            if (a[k].q_pos < qmin) qmin = a[k].q_pos;
            if (a[k].q_pos > qmax) qmax = a[k].q_pos;
            if (a[k].t_pos < tmin) tmin = a[k].t_pos;
            if (a[k].t_pos > tmax) tmax = a[k].t_pos;
        }
        CHECK(tags_ok, "all anchor tags valid");
        CHECK(q_ordered, "anchors ordered by query position");
        CHECK(t_monotone, "target positions monotone for the overlap's strand");
        /* Bounding coords bracket the chained anchors. */
        CHECK(o.q_start <= qmin && qmax < o.q_end, "q bounds bracket anchors");
        CHECK(o.t_start <= tmin && tmax < o.t_end, "t bounds bracket anchors");

        if (g_fail) break;   /* stop on first structural failure for a clean log */
    }
    CHECK(total_anchor_seen == n_anchors, "arena fully covered by overlap slices");

    fakechain_bridge_free(ov, anchors, names, name_off);

    if (g_fail) { std::fprintf(stderr, "\n%d CHECK(s) FAILED\n", g_fail); return 1; }
    std::fprintf(stderr, "\nfakechain bridge end-to-end smoke test passed.\n");
    return 0;
}
