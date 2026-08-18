/*
 * Unit test for hifiasm_chain_pair(): a pure, data-free check of the DP wrapper.
 *
 * Feeds hand-built anchor sets and asserts the chaining contract:
 *   - a clean colinear run is returned intact, ordered by query position;
 *   - an off-diagonal outlier is dropped from the kept chain;
 *   - the opaque `id` rides through onto the kept anchors;
 *   - empty input is the empty-chain case (rc 0, 0 kept), not an error;
 *   - argument validation rejects NULL / zero-length reads.
 *
 * No read files, no global stores -- runs anywhere in milliseconds.
 *
 * Build/run:  make test_hifiasm_chain && ./test_hifiasm_chain
 */
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>

#include "hifiasm_chain.h"

static int g_fail = 0;
#define CHECK(c,m) do{ if(!(c)){ std::fprintf(stderr,"FAIL: %s\n",(m)); ++g_fail; } }while(0)

int main(void) {
    const uint32_t QLEN = 20000, TLEN = 20000;

    /* --- 1. Clean colinear diagonal: q and t advance together by ~500. --- */
    {
        std::vector<hifiasm_chain_anchor_t> in;
        for (uint32_t i = 0; i < 12; ++i) {
            hifiasm_chain_anchor_t a;
            a.q_pos = 1000 + i * 500;
            a.t_pos = 3000 + i * 500;   /* constant +2000 offset (colinear) */
            a.id    = 0xA000ull + i;
            in.push_back(a);
        }
        hifiasm_chain_anchor_t *kept = NULL;
        uint64_t n_kept = 0; int32_t score = 0;
        int rc = hifiasm_chain_pair(in.data(), in.size(), QLEN, TLEN, NULL,
                                    &kept, &n_kept, &score);
        CHECK(rc == 0, "colinear: rc 0");
        CHECK(n_kept == in.size(), "colinear: all anchors kept");
        CHECK(kept != NULL, "colinear: kept array non-null");
        if (kept) {
            for (uint64_t i = 1; i < n_kept; ++i) {
                CHECK(kept[i].q_pos > kept[i-1].q_pos, "colinear: q strictly increasing");
                CHECK(kept[i].t_pos > kept[i-1].t_pos, "colinear: t strictly increasing");
            }
            /* id preserved: first kept anchor should be the first input id. */
            CHECK(kept[0].id == 0xA000ull, "colinear: opaque id preserved");
            CHECK(score > 0, "colinear: positive chain score");
        }
        hifiasm_chain_free(kept);
    }

    /* --- 2. Colinear run + one off-diagonal outlier that breaks the band. --- */
    {
        std::vector<hifiasm_chain_anchor_t> in;
        for (uint32_t i = 0; i < 10; ++i) {
            hifiasm_chain_anchor_t a;
            a.q_pos = 1000 + i * 500;
            a.t_pos = 1000 + i * 500;   /* main diagonal */
            a.id    = i;
            in.push_back(a);
        }
        /* Outlier: query mid-range but target far off the diagonal. Its huge
         * gap should exclude it from the best chain. */
        hifiasm_chain_anchor_t bad;
        bad.q_pos = 3200;               /* between anchors 4 and 5 on q      */
        bad.t_pos = 15000;              /* ~11.5k off the diagonal           */
        bad.id    = 0xBAD;
        in.push_back(bad);

        hifiasm_chain_anchor_t *kept = NULL;
        uint64_t n_kept = 0; int32_t score = 0;
        int rc = hifiasm_chain_pair(in.data(), in.size(), QLEN, TLEN, NULL,
                                    &kept, &n_kept, &score);
        CHECK(rc == 0, "outlier: rc 0");
        CHECK(n_kept >= 10, "outlier: main diagonal chain kept");
        bool kept_bad = false;
        for (uint64_t i = 0; i < n_kept; ++i) {
            if (kept && kept[i].id == 0xBAD) kept_bad = true;
        }
        CHECK(!kept_bad, "outlier: off-diagonal anchor dropped");
        if (kept) {
            for (uint64_t i = 1; i < n_kept; ++i) {
                CHECK(kept[i].q_pos > kept[i-1].q_pos, "outlier: q strictly increasing");
                CHECK(kept[i].t_pos > kept[i-1].t_pos, "outlier: t strictly increasing");
            }
        }
        hifiasm_chain_free(kept);
    }

    /* --- 3. Empty input: empty-chain case, not an error. --- */
    {
        hifiasm_chain_anchor_t *kept = (hifiasm_chain_anchor_t *)0x1;
        uint64_t n_kept = 99; int32_t score = 7;
        hifiasm_chain_anchor_t dummy = {0,0,0};
        int rc = hifiasm_chain_pair(&dummy, 0, QLEN, TLEN, NULL,
                                    &kept, &n_kept, &score);
        CHECK(rc == 0, "empty: rc 0");
        CHECK(n_kept == 0, "empty: 0 kept");
        CHECK(kept == NULL, "empty: kept cleared to NULL");
        CHECK(score == 0, "empty: score cleared");
    }

    /* --- 4. Argument validation. --- */
    {
        hifiasm_chain_anchor_t a = {0,0,0};
        hifiasm_chain_anchor_t *kept = NULL; uint64_t n_kept = 0;
        CHECK(hifiasm_chain_pair(NULL, 1, 100, 100, NULL, &kept, &n_kept, NULL) != 0,
              "validation: NULL anchors rejected");
        CHECK(hifiasm_chain_pair(&a, 1, 0, 100, NULL, &kept, &n_kept, NULL) != 0,
              "validation: zero q_len rejected");
        CHECK(hifiasm_chain_pair(&a, 1, 100, 0, NULL, &kept, &n_kept, NULL) != 0,
              "validation: zero t_len rejected");
        CHECK(hifiasm_chain_pair(&a, 1, 100, 100, NULL, NULL, &n_kept, NULL) != 0,
              "validation: NULL out_kept rejected");
    }

    if (g_fail) {
        std::fprintf(stderr, "hifiasm_chain: %d check(s) FAILED\n", g_fail);
        return 1;
    }
    std::fprintf(stderr, "hifiasm_chain: all checks passed\n");
    return 0;
}
