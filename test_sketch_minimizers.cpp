/*
 * Standalone test for hifiasm_sketch_minimizers().
 *
 * Validates the no-HPC minimizer bridge that dinara consumes:
 *   - every position is in range and START-anchored (pos + span <= len)
 *   - no-HPC span always equals k; HPC span >= k
 *   - output is sorted by START position and deduplicated of exact repeats
 *   - the rev/hash pair is self-consistent: re-hashing the canonical k-mer at
 *     [pos, pos+k) with hifiasm's own yak hash reproduces the reported hash
 *   - the call is deterministic and handles edge cases (short seq, N runs)
 *
 * Build (from the submodule root, after `make lib`):
 *   g++ -std=c++11 -O2 -I. test_sketch_minimizers.cpp \
 *       libhifiasm_overlaps.a -lz -lpthread -lm -o test_sketch_minimizers
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "hifiasm_overlaps.h"

/* hifiasm's base->2bit table and yak hash, declared here so the test can
 * independently recompute the canonical hash and cross-check the bridge. */
extern "C" const unsigned char seq_nt4_table[256];
static inline uint64_t yak_hash64_64(uint64_t key)
{
    key = ~key + (key << 21);
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key = key ^ key >> 28;
    key = key + (key << 31);
    return key;
}

static int g_failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); \
        ++g_failures; \
    } \
} while (0)

/* Recompute hifiasm's canonical hash for the k-mer at forward [pos, pos+k).
 * Mirrors sketch.cpp: rolling fwd/rev 2-bit encodings, canonical strand is the
 * lexicographically smaller of the two encodings, hash = h64(fwd)+h64(rev).
 * Returns false if the window contains an ambiguous base or is symmetric. */
static bool canonical_hash_at(const std::string& s, uint32_t pos, int k,
                              uint64_t* out_hash, uint32_t* out_rev)
{
    const uint64_t mask = (k >= 64) ? ~0ULL : ((1ULL << k) - 1);
    const uint64_t shift1 = (uint64_t)(k - 1);
    uint64_t fwd = 0, rev = 0;
    for (int i = 0; i < k; ++i) {
        int c = seq_nt4_table[(uint8_t)s[pos + i]];
        if (c >= 4) return false;
        fwd = ((fwd << 1) | (c & 1)) & mask;
        rev = (rev >> 1) | ((uint64_t)(1 - (c & 1)) << shift1);
    }
    /* The 2-bit-per-base scheme in sketch.cpp splits each base into two 1-bit
     * lanes; the canonical strand test there is on the second lane (kmer[1] vs
     * kmer[3]). Rather than replicate both lanes, we recompute exactly as the
     * sketch does using two independent bit lanes. */
    uint64_t f0 = 0, f1 = 0, r0 = 0, r1 = 0;
    for (int i = 0; i < k; ++i) {
        int c = seq_nt4_table[(uint8_t)s[pos + i]];
        f0 = ((f0 << 1) | (c & 1)) & mask;
        f1 = ((f1 << 1) | (c >> 1)) & mask;
        r0 = (r0 >> 1) | ((uint64_t)(1 - (c & 1)) << shift1);
        r1 = (r1 >> 1) | ((uint64_t)(1 - (c >> 1)) << shift1);
    }
    if (f1 == r1) return false; /* symmetric: sketch skips these */
    uint32_t z = (f1 < r1) ? 0 : 1;
    uint64_t hf = yak_hash64_64(z == 0 ? f0 : r0);
    uint64_t hr = yak_hash64_64(z == 0 ? f1 : r1);
    /* sketch.cpp: y = h64(kmer[z<<1|0]) + h64(kmer[z<<1|1]) */
    *out_hash = hf + hr;
    *out_rev = z;
    return true;
}

/* Build a reproducible pseudo-random ACGT sequence. */
static std::string random_seq(size_t n, uint32_t seed)
{
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.resize(n);
    uint32_t x = seed ? seed : 1;
    for (size_t i = 0; i < n; ++i) {
        x ^= x << 13; x ^= x >> 17; x ^= x << 5; /* xorshift32 */
        s[i] = bases[x & 3];
    }
    return s;
}

static void test_basic_invariants()
{
    const int k = 50, w = 50;
    std::string s = random_seq(5000, 12345);

    hifiasm_minimizer_t* mz = nullptr;
    int n = 0;
    int rc = hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k,
                                       /*is_hpc*/ 0, &mz, &n);
    CHECK(rc == 0, "no-HPC sketch returns 0");
    CHECK(n > 0, "no-HPC sketch produces minimizers");

    uint32_t prev = 0;
    bool first = true;
    for (int i = 0; i < n; ++i) {
        CHECK(mz[i].span == (uint32_t)k, "no-HPC span equals k");
        CHECK(mz[i].pos + mz[i].span <= s.size(), "k-mer fits in sequence");
        CHECK(mz[i].rev <= 1, "rev is 0 or 1");
        if (!first) CHECK(mz[i].pos >= prev, "positions are sorted ascending");
        prev = mz[i].pos;
        first = false;

        /* Cross-check hash/rev against an independent recomputation. */
        uint64_t exp_hash = 0; uint32_t exp_rev = 0;
        if (canonical_hash_at(s, mz[i].pos, k, &exp_hash, &exp_rev)) {
            CHECK(mz[i].hash == exp_hash, "hash matches recomputed canonical hash");
            CHECK(mz[i].rev == exp_rev, "rev matches recomputed canonical strand");
        }
    }
    free(mz);
}

static void test_determinism()
{
    const int k = 31, w = 31;
    std::string s = random_seq(3000, 99);

    hifiasm_minimizer_t* a = nullptr; int na = 0;
    hifiasm_minimizer_t* b = nullptr; int nb = 0;
    hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k, 0, &a, &na);
    hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k, 0, &b, &nb);
    CHECK(na == nb, "repeat calls produce same count");
    if (na == nb) {
        bool same = true;
        for (int i = 0; i < na; ++i) {
            if (a[i].pos != b[i].pos || a[i].hash != b[i].hash ||
                a[i].rev != b[i].rev || a[i].span != b[i].span) {
                same = false; break;
            }
        }
        CHECK(same, "repeat calls produce identical output");
    }
    free(a); free(b);
}

static void test_hpc_span()
{
    /* A sequence with long homopolymer runs: under HPC, some spans must exceed
     * k; under no-HPC, all spans equal k. */
    std::string s;
    for (int i = 0; i < 60; ++i) s += "AAAAACCCCCGGGGGTTTTT";
    const int k = 15, w = 15;

    hifiasm_minimizer_t* raw = nullptr; int nraw = 0;
    hifiasm_minimizer_t* hpc = nullptr; int nhpc = 0;
    hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k, 0, &raw, &nraw);
    hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k, 1, &hpc, &nhpc);

    for (int i = 0; i < nraw; ++i)
        CHECK(raw[i].span == (uint32_t)k, "no-HPC span == k on homopolymers");

    bool any_wide = false;
    for (int i = 0; i < nhpc; ++i) {
        CHECK(hpc[i].span >= (uint32_t)k, "HPC span >= k");
        CHECK(hpc[i].pos + hpc[i].span <= s.size(), "HPC k-mer fits in sequence");
        if (hpc[i].span > (uint32_t)k) any_wide = true;
    }
    CHECK(any_wide, "HPC produces at least one span > k on homopolymer input");
    free(raw); free(hpc);
}

static void test_edge_cases()
{
    hifiasm_minimizer_t* mz = (hifiasm_minimizer_t*)0x1; int n = 7;

    /* Sequence shorter than k -> zero minimizers, clean outputs. */
    std::string shortSeq = random_seq(10, 5);
    int rc = hifiasm_sketch_minimizers(shortSeq.c_str(), (int)shortSeq.size(),
                                       50, 50, 0, &mz, &n);
    CHECK(rc == 0, "seq shorter than k returns 0");
    CHECK(mz == nullptr, "seq shorter than k sets out array to NULL");
    CHECK(n == 0, "seq shorter than k sets count to 0");

    /* Invalid parameters are rejected, not asserted. */
    std::string s = random_seq(200, 7);
    rc = hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), 0, 50, 0, &mz, &n);
    CHECK(rc != 0, "w=0 is rejected");
    rc = hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), 50, 64, 0, &mz, &n);
    CHECK(rc != 0, "k>63 is rejected");
    rc = hifiasm_sketch_minimizers(nullptr, 100, 50, 50, 0, &mz, &n);
    CHECK(rc != 0, "null seq is rejected");

    /* Ambiguous bases (N) break k-mers but must not crash or emit bad spans. */
    std::string withN = random_seq(1000, 3);
    for (int i = 100; i < 130; ++i) withN[i] = 'N';
    mz = nullptr; n = 0;
    rc = hifiasm_sketch_minimizers(withN.c_str(), (int)withN.size(),
                                   31, 31, 0, &mz, &n);
    CHECK(rc == 0, "seq with N returns 0");
    for (int i = 0; i < n; ++i)
        CHECK(mz[i].pos + mz[i].span <= withN.size(), "N-adjacent k-mer in range");
    free(mz);
}

/* The context path must produce exactly the same output as the one-shot path,
 * be reusable across many reads, and deduplicate by START position. */
static void test_ctx_equivalence()
{
    const int k = 50, w = 50;
    hifiasm_sketch_ctx_t* ctx = hifiasm_sketch_ctx_init();
    CHECK(ctx != nullptr, "ctx init succeeds");
    if (!ctx) return;

    for (uint32_t seed = 1; seed <= 20; ++seed) {
        std::string s = random_seq(3000 + seed * 100, seed * 7919u + 1);

        hifiasm_minimizer_t* one = nullptr; int n1 = 0;
        hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k, 0, &one, &n1);

        const hifiasm_minimizer_t* cm = nullptr; int n2 = 0;
        int rc = hifiasm_sketch_minimizers_ctx(ctx, s.c_str(), (int)s.size(),
                                               w, k, 0, &cm, &n2);
        CHECK(rc == 0, "ctx sketch returns 0");
        CHECK(n1 == n2, "ctx and one-shot produce same count");
        if (n1 == n2) {
            bool same = true;
            for (int i = 0; i < n1; ++i) {
                if (one[i].pos != cm[i].pos || one[i].hash != cm[i].hash ||
                    one[i].rev != cm[i].rev || one[i].span != cm[i].span) {
                    same = false; break;
                }
            }
            CHECK(same, "ctx and one-shot produce identical output");
        }
        /* Output must be strictly increasing in START position (sorted+deduped). */
        for (int i = 1; i < n2; ++i)
            CHECK(cm[i].pos > cm[i - 1].pos, "ctx output strictly increasing (deduped)");

        free(one);
    }
    hifiasm_sketch_ctx_destroy(ctx);
    hifiasm_sketch_ctx_destroy(nullptr); /* NULL is a no-op */
}

/* A homopolymer-heavy sequence in no-HPC mode can produce several k-mers that
 * are identical in hash but the minimizer scheme selects distinct START
 * positions; ensure the deduped output never contains a repeated START. */
static void test_dedup_positions()
{
    std::string s;
    for (int i = 0; i < 40; ++i) s += "ACACACACGTGTGTGTAACCGGTT";
    const int k = 21, w = 21;
    hifiasm_minimizer_t* mz = nullptr; int n = 0;
    hifiasm_sketch_minimizers(s.c_str(), (int)s.size(), w, k, 0, &mz, &n);
    for (int i = 1; i < n; ++i)
        CHECK(mz[i].pos > mz[i - 1].pos, "one-shot output has no duplicate START");
    free(mz);
}

int main()
{
    test_basic_invariants();
    test_determinism();
    test_hpc_span();
    test_edge_cases();
    test_ctx_equivalence();
    test_dedup_positions();

    if (g_failures == 0) {
        std::printf("ALL TESTS PASSED\n");
        return 0;
    }
    std::printf("%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
