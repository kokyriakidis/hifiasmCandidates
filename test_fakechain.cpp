/*
 * Unit tests for fakechain_pair (fakechain.cpp).
 *
 * fakechain_pair now: clips myloasm keys to an even marker_k (default 20),
 * MATCHES markers between a read pair (m:n, canonical key equality, syncmers +
 * SNPmers in separate key spaces), merges them into one tagged anchor set, and
 * CHAINS them with myloasm's DP (myloasm_chain FFI). The output is the ordered,
 * tagged chained-anchor array plus counts and bounding coords.
 *
 * These tests use SYNTHETIC markers so we control every key/pos. Because the DP
 * default min_chain_length is 3, most fixtures use >= 3 colinear anchors. For
 * perfectly diagonal anchors (dist_q == dist_t) each link adds match_score(1),
 * so the chain score equals the number of chained anchors.
 *
 * Properties covered:
 *   1. Same-strand chaining: syncmers + SNPmers chain together, ordered by
 *      query pos, tagged; score == n_anchor for a clean diagonal.
 *   2. RC overlap: identical canonical keys with target positions DECREASING
 *      chain (myloasm's strand-bit encoding); n_anchor == matches.
 *   3. k21->k20 clip keeps SNP alleles distinct: mismatched middle base does
 *      NOT match (no SNPmer anchor); matching allele DOES (one SNPmer anchor).
 *      Last-base differences ARE clipped away (they still match).
 *   4. m:n repeats kept pre-chain: a key repeated on the target yields multiple
 *      matches (n_match, max_mult), and the DP selects the colinear one.
 *   5. Interval restriction: markers outside [s,e) are not matched.
 *   6. Chain-length gating: 2 colinear anchors give no chain at default
 *      min_chain_len=3, but do chain when min_chain_len=2.
 *   7. min_anchors gating (host side): a real chain below min_anchors sets ok=0.
 *
 * Build/run:  make test_fakechain && ./test_fakechain
 */
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>

#include "fakechain.h"

static int g_fail = 0;
#define CHECK(c,m) do{ if(!(c)){ std::fprintf(stderr,"FAIL: %s\n",(m)); ++g_fail; } \
                       else { std::fprintf(stderr,"ok: %s\n",(m)); } }while(0)

/* Build a k=21 key from 21 base codes (0..3), base i at bits [2i,2i+1]. We only
 * need DETERMINISTIC keys and control over the middle base (index 10). */
static uint64_t mk_key(const int b[21]) {
    uint64_t k = 0;
    for (int i = 0; i < 21; ++i) k |= ((uint64_t)(b[i] & 3)) << (2 * i);
    return k;
}

static MyloMarker M(uint32_t pos, uint64_t key) {
    MyloMarker m; m.pos = pos; m._pad = 0; m.key = key; return m;
}

int main(void) {
    /* ---- Fixtures: a handful of distinct 21-mer keys. ---- */
    int base[21];
    for (int i = 0; i < 21; ++i) base[i] = (i * 7 + 1) & 3;   /* arbitrary */
    uint64_t kA = mk_key(base);

    base[3] ^= 1; uint64_t kB = mk_key(base); base[3] ^= 1;
    base[7] ^= 2; uint64_t kC = mk_key(base); base[7] ^= 2;
    base[15] ^= 3; uint64_t kD = mk_key(base); base[15] ^= 3;

    /* Two SNP alleles: differ ONLY in the middle base (index 10). */
    int allele[21]; memcpy(allele, base, sizeof(base));
    allele[10] = 0; uint64_t snp0 = mk_key(allele);
    allele[10] = 2; uint64_t snp1 = mk_key(allele);  /* diff middle base */

    /* Key differing only in the LAST base (index 20, bits [40,41]). */
    int last[21]; memcpy(last, base, sizeof(base));
    last[20] ^= 1; uint64_t kLast = mk_key(last);

    /* --- Property 3a: the clip itself keeps middle, drops last (bit math). --- */
    {
        const uint64_t mask20 = (1ull << 40) - 1ull;
        CHECK((snp0 & mask20) != (snp1 & mask20),
              "k21->k20 clip keeps SNP alleles distinct (middle base retained)");
        CHECK((kA & mask20) == (kLast & mask20),
              "k21->k20 clip removes the last base (bits [40,41])");
    }

    /* --- Property 1: same-strand, syncmers + SNPmer chain together. --- */
    {
        /* Diagonal, target offset +1000. 3 syncmer + 1 snpmer = 4 anchors. */
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB), M(300,kC) };
        std::vector<MyloMarker> qp = { M(150,snp0) };
        std::vector<MyloMarker> ts = { M(1100,kA), M(1200,kB), M(1300,kC) };
        std::vector<MyloMarker> tp = { M(1150,snp0) };

        fakechain_result_t r;
        int rc = fakechain_pair(qs.data(),qs.size(), qp.data(),qp.size(),
                                ts.data(),ts.size(), tp.data(),tp.size(),
                                90,400, 1090,1400, /*same_strand=*/1, NULL, &r);
        CHECK(rc==0, "same-strand: returns success");
        CHECK(r.ok==1, "same-strand: ok");
        CHECK(r.n_anchor==4, "same-strand: 4 chained anchors (3 syncmer + 1 snpmer)");
        CHECK(r.n_snpmer_anchor==1, "same-strand: 1 snpmer anchor");
        CHECK(r.n_match==4, "same-strand: 4 matches fed to DP");
        CHECK(r.max_mult==1, "same-strand: max_mult == 1 (unique keys)");
        CHECK(r.score==4, "same-strand: diagonal chain score == n_anchor");
        CHECK(r.anchors!=NULL, "same-strand: anchors array allocated");
        /* Ordered by query pos ascending. */
        bool ordered = true, has_snp = false;
        for (uint32_t i = 0; i < r.n_anchor; ++i) {
            if (i && r.anchors[i].q_pos <= r.anchors[i-1].q_pos) ordered = false;
            if (r.anchors[i].tag == FAKECHAIN_SNPMER) has_snp = true;
        }
        CHECK(ordered, "same-strand: anchors ordered by query pos");
        CHECK(has_snp, "same-strand: SNPmer anchor tagged in output");
        /* SNPmer anchor carries the clipped snp0 key. */
        const uint64_t mask20 = (1ull << 40) - 1ull;
        bool snp_key_ok = false;
        for (uint32_t i = 0; i < r.n_anchor; ++i)
            if (r.anchors[i].tag == FAKECHAIN_SNPMER)
                snp_key_ok = (r.anchors[i].key == (snp0 & mask20));
        CHECK(snp_key_ok, "same-strand: SNPmer anchor key is clipped snp0");
        fakechain_result_free(&r);
    }

    /* --- Property 2: RC overlap, target positions DECREASING. --- */
    {
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB), M(300,kC) };
        std::vector<MyloMarker> ts = { M(1300,kA), M(1200,kB), M(1100,kC) };
        fakechain_result_t r;
        int rc = fakechain_pair(qs.data(),qs.size(), NULL,0,
                                ts.data(),ts.size(), NULL,0,
                                90,400, 1090,1400, /*same_strand=*/0, NULL, &r);
        CHECK(rc==0 && r.ok==1, "RC: success + ok");
        CHECK(r.n_anchor==3, "RC: all 3 canonical keys chained across strands");
        CHECK(r.is_same_strand==0, "RC: strand flag carried through");
        CHECK(r.score==3, "RC: diagonal chain score == 3");
        fakechain_result_free(&r);
    }

    /* --- Property 3b: mismatched SNP allele does NOT anchor; matching does. --- */
    {
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB), M(300,kC) };
        std::vector<MyloMarker> ts = { M(1100,kA), M(1200,kB), M(1300,kC) };
        std::vector<MyloMarker> qp = { M(150,snp0) };

        /* Wrong allele on target -> no SNPmer anchor. */
        std::vector<MyloMarker> tp_wrong = { M(1150,snp1) };
        fakechain_result_t rw;
        fakechain_pair(qs.data(),qs.size(), qp.data(),qp.size(),
                       ts.data(),ts.size(), tp_wrong.data(),tp_wrong.size(),
                       90,400, 1090,1400, 1, NULL, &rw);
        CHECK(rw.ok==1 && rw.n_anchor==3 && rw.n_snpmer_anchor==0,
              "clip: mismatched middle base -> no SNPmer anchor (3 syncmer only)");
        fakechain_result_free(&rw);

        /* Matching allele on target -> one SNPmer anchor. */
        std::vector<MyloMarker> tp_ok = { M(1150,snp0) };
        fakechain_result_t rok;
        fakechain_pair(qs.data(),qs.size(), qp.data(),qp.size(),
                       ts.data(),ts.size(), tp_ok.data(),tp_ok.size(),
                       90,400, 1090,1400, 1, NULL, &rok);
        CHECK(rok.n_anchor==4 && rok.n_snpmer_anchor==1,
              "clip: matching middle base -> one SNPmer anchor");
        fakechain_result_free(&rok);

        /* Last-base difference is clipped away -> still matches as a syncmer. */
        std::vector<MyloMarker> qs2 = { M(100,kLast), M(200,kB), M(300,kC) };
        std::vector<MyloMarker> ts2 = { M(1100,kA),   M(1200,kB), M(1300,kC) };
        fakechain_result_t rl;
        fakechain_pair(qs2.data(),qs2.size(), NULL,0,
                       ts2.data(),ts2.size(), NULL,0,
                       90,400, 1090,1400, 1, NULL, &rl);
        CHECK(rl.ok==1 && rl.n_anchor==3,
              "clip: last-base difference clipped -> kLast matches kA (3 anchors)");
        fakechain_result_free(&rl);
    }

    /* --- Property 4: m:n repeats kept pre-chain; DP picks the colinear one. --- */
    {
        /* kB repeated on target: one on-diagonal (1200), one off (1250). */
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB), M(300,kC) };
        std::vector<MyloMarker> ts = { M(1100,kA), M(1200,kB), M(1250,kB), M(1300,kC) };
        fakechain_result_t r;
        int rc = fakechain_pair(qs.data(),qs.size(), NULL,0,
                                ts.data(),ts.size(), NULL,0,
                                90,400, 1090,1400, 1, NULL, &r);
        CHECK(rc==0, "m:n: success");
        CHECK(r.n_match==4, "m:n: repeated key kept -> 4 matches fed to DP");
        CHECK(r.max_mult==2, "m:n: max_mult == 2 (kB twice on target)");
        CHECK(r.n_anchor==3, "m:n: DP selects the colinear subset (3 anchors)");
        /* The kept kB anchor is the on-diagonal one (target 1200). */
        bool picked_diag = false;
        for (uint32_t i = 0; i < r.n_anchor; ++i)
            if (r.anchors[i].q_pos == 200) picked_diag = (r.anchors[i].t_pos == 1200);
        CHECK(picked_diag, "m:n: DP kept the on-diagonal kB (t==1200)");
        fakechain_result_free(&r);
    }

    /* --- Property 5: interval restriction. --- */
    {
        /* kD sits OUTSIDE the query interval [90,400) (pos 5000); ignored. */
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB), M(300,kC), M(5000,kD) };
        std::vector<MyloMarker> ts = { M(1100,kA), M(1200,kB), M(1300,kC), M(6000,kD) };
        fakechain_result_t r;
        int rc = fakechain_pair(qs.data(),qs.size(), NULL,0,
                                ts.data(),ts.size(), NULL,0,
                                90,400, 1090,1400, 1, NULL, &r);
        CHECK(rc==0, "interval: success");
        CHECK(r.n_match==3, "interval: out-of-interval marker not matched (3 matches)");
        CHECK(r.n_anchor==3, "interval: chain over in-interval anchors only");
        fakechain_result_free(&r);
    }

    /* --- Property 6: chain-length gating. --- */
    {
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB) };
        std::vector<MyloMarker> ts = { M(1100,kA), M(1200,kB) };

        /* Default min_chain_len=3 -> 2 anchors do not form a chain. */
        fakechain_result_t r;
        int rc = fakechain_pair(qs.data(),qs.size(), NULL,0,
                                ts.data(),ts.size(), NULL,0,
                                90,400, 1090,1400, 1, NULL, &r);
        CHECK(rc==0 && r.ok==0 && r.n_anchor==0,
              "chain-len: 2 anchors < default min_chain_len(3) -> no chain");
        CHECK(r.n_match==2, "chain-len: both matches were fed to the DP");
        fakechain_result_free(&r);

        /* Lower min_chain_len to 2 -> the 2-anchor chain is returned. */
        fakechain_opt_t opt; fakechain_opt_init(&opt); opt.min_chain_len = 2;
        fakechain_result_t r2;
        fakechain_pair(qs.data(),qs.size(), NULL,0, ts.data(),ts.size(), NULL,0,
                       90,400, 1090,1400, 1, &opt, &r2);
        CHECK(r2.ok==1 && r2.n_anchor==2,
              "chain-len: min_chain_len=2 -> 2-anchor chain returned");
        fakechain_result_free(&r2);
    }

    /* --- Property 7: host-side min_anchors gating. --- */
    {
        std::vector<MyloMarker> qs = { M(100,kA), M(200,kB), M(300,kC) };
        std::vector<MyloMarker> ts = { M(1100,kA), M(1200,kB), M(1300,kC) };
        fakechain_opt_t opt; fakechain_opt_init(&opt); opt.min_anchors = 5;
        fakechain_result_t r;
        fakechain_pair(qs.data(),qs.size(), NULL,0, ts.data(),ts.size(), NULL,0,
                       90,400, 1090,1400, 1, &opt, &r);
        CHECK(r.n_anchor==3 && r.ok==0,
              "min_anchors: chain below threshold -> ok==0 (anchors still emitted)");
        fakechain_result_free(&r);
    }

    /* --- Robustness: free on a zeroed / no-anchor result is safe. --- */
    {
        std::vector<MyloMarker> qs = { M(100,kA) };
        std::vector<MyloMarker> ts = { M(9000,kB) };   /* no shared key */
        fakechain_result_t r;
        int rc = fakechain_pair(qs.data(),qs.size(), NULL,0,
                                ts.data(),ts.size(), NULL,0,
                                90,400, 8090,9400, 1, NULL, &r);
        CHECK(rc==0 && r.ok==0 && r.n_anchor==0 && r.anchors==NULL,
              "no-match: ok==0, no anchors, NULL array");
        fakechain_result_free(&r);   /* must be safe */
        fakechain_result_free(&r);   /* idempotent */
    }

    if (g_fail) { std::fprintf(stderr, "\n%d CHECK(s) FAILED\n", g_fail); return 1; }
    std::fprintf(stderr, "\nAll fakechain tests passed.\n");
    return 0;
}
