/*
 * Standalone test for the overlap-parity marker filter:
 *   hifiasm_build_filter() + hifiasm_sketch_minimizers_ctx_filtered()
 *
 * How the filter works (this drove the test design):
 *   hifiasm's minimizer comparator (sketch.cpp mz1_mzcmp) orders k-mers by
 *   their occurrence count FIRST, then by hash. During filtered sketching each
 *   candidate's count comes from the hf table (info.rid = cnt). So with hf the
 *   window minimizer is the RAREST k-mer in the window (lowest count, ties by
 *   hash); without hf every count is 0 and it degenerates to the ordinary
 *   lowest-hash minimizer. A separate hard filter drops k-mers whose count
 *   exceeds max_kmer_cnt outright. Distance subsampling (select_mz_h) then
 *   thins survivors when sample_dist > w.
 *
 * Because that comparator/selection logic lives entirely inside hifiasm's
 * mz1_ha_sketch, the strongest correctness statement is PARITY: the bridge must
 * reproduce, bit for bit, the seed set that mz1_ha_sketch produces when called
 * with the SAME arguments hifiasm's overlap path uses (anchor.cpp:110 ->
 * hf=ha_flt_tab, sample_dist=mz_sample_dist, ws=mz_rewin, pt=NULL). This test
 * calls mz1_ha_sketch directly as the oracle and checks the bridge matches it,
 * plus validity, determinism, and the hf=NULL passthrough.
 *
 * Build (from the submodule root, after `make lib`):
 *   g++ -std=c++11 -O2 -I. test_sketch_filter.cpp \
 *       libhifiasm_overlaps.a -lz -lpthread -lm -o test_sketch_filter
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unistd.h>

#include "hifiasm_overlaps.h"
#include "htab.h" /* ha_mz1_t, ha_mz1_v, st_mt_t, mz1_ha_sketch */

/* mz_rewin default hifiasm's overlap path uses for select_mz_h; the bridge
 * hardcodes the same value, so the oracle must use it too. */
#define ORACLE_WS 1000

static int g_failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); \
        ++g_failures; \
    } \
} while (0)

/* Reproducible pseudo-random ACGT sequence (xorshift32). */
static std::string random_seq(size_t n, uint32_t seed)
{
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.resize(n);
    uint32_t x = seed ? seed : 1;
    for (size_t i = 0; i < n; ++i) {
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        s[i] = bases[x & 3];
    }
    return s;
}

/*
 * Genome: large unique backbone (peak_hom tracks coverage) + a short repeat
 * unit at high copy number (its k-mers become high-count). Tiled at ~8x into a
 * temp FASTA. Fills *repeat_out with the repeat unit.
 */
static std::string write_genome_fasta(std::string* repeat_out)
{
    const int backbone_blocks = 40, block_len = 10000, repeat_copies = 60;
    std::string repeat = random_seq(200, 0xABCD);
    *repeat_out = repeat;

    std::string genome;
    for (int i = 0; i < backbone_blocks; ++i) {
        genome += random_seq(block_len, 7000u + (uint32_t)i);
        if (i < repeat_copies) genome += repeat;
    }
    for (int i = backbone_blocks; i < repeat_copies; ++i)
        genome += repeat + random_seq(500, 8000u + (uint32_t)i);

    char tmpl[] = "/tmp/hifiasm_flt_test_XXXXXX";
    int fd = mkstemp(tmpl);
    if (fd < 0) { perror("mkstemp"); std::exit(2); }
    FILE* fp = fdopen(fd, "w");
    if (!fp) { perror("fdopen"); std::exit(2); }

    const int win = 6000, step = 750;
    int rid = 0;
    for (size_t p = 0; p + win <= genome.size(); p += step)
        std::fprintf(fp, ">r%d\n%s\n", rid++, genome.substr(p, win).c_str());
    std::fclose(fp);
    return std::string(tmpl);
}

static void check_valid_sorted_dedup(const hifiasm_minimizer_t* mz, int n,
                                     int k, size_t len)
{
    bool first = true;
    uint32_t prev = 0;
    std::set<uint32_t> seen;
    for (int i = 0; i < n; ++i) {
        CHECK(mz[i].span == (uint32_t)k, "span==k (no-HPC)");
        CHECK((size_t)mz[i].pos + mz[i].span <= len, "k-mer fits in sequence");
        CHECK(mz[i].rev <= 1, "rev is 0 or 1");
        if (!first) CHECK(mz[i].pos >= prev, "positions sorted ascending");
        CHECK(seen.insert(mz[i].pos).second, "positions deduplicated");
        prev = mz[i].pos; first = false;
    }
}

/*
 * Oracle: call mz1_ha_sketch exactly as hifiasm's overlap path does, then
 * convert to START-anchored (pos,hash) pairs the way the bridge does, so the
 * two can be compared directly. Returns a pos->hash map (the bridge dedups by
 * START position, and exactly one k-mer starts at any position, so a map is the
 * right shape).
 */
static std::map<uint32_t, uint64_t>
oracle_seeds(const std::string& read, int w, int k,
             const void* raw_hf, int sample_dist)
{
    ha_mz1_v p; std::memset(&p, 0, sizeof(p));
    st_mt_t  mt; std::memset(&mt, 0, sizeof(mt));
    mz1_ha_sketch(read.c_str(), (int)read.size(), w, k, /*rid*/ 0,
                  /*is_hpc*/ 0, &p, raw_hf, sample_dist, /*k_flag*/ NULL,
                  /*dbg_ct*/ NULL, /*pt*/ NULL, /*min_freq*/ -1,
                  /*dp_min_len*/ -1, /*dp_e*/ -1.0f, &mt, /*ws*/ ORACLE_WS,
                  /*is_unique*/ 0, /*km*/ NULL);
    std::map<uint32_t, uint64_t> m;
    for (uint32_t i = 0; i < p.n; ++i) {
        uint32_t span = (uint32_t)p.a[i].span;
        uint32_t end  = (uint32_t)p.a[i].pos;
        m[end + 1u - span] = p.a[i].x; /* START -> hash */
    }
    free(p.a); free(mt.a);
    return m;
}

static bool maps_equal(const std::map<uint32_t,uint64_t>& a,
                       const hifiasm_minimizer_t* mz, int n)
{
    if ((int)a.size() != n) return false;
    for (int i = 0; i < n; ++i) {
        auto it = a.find(mz[i].pos);
        if (it == a.end() || it->second != mz[i].hash) return false;
    }
    return true;
}

int main(void)
{
    const int k = 50, w = 50, sample_dist = 500;

    std::string repeat;
    std::string fasta = write_genome_fasta(&repeat);
    const char* files[1] = { fasta.c_str() };

    hifiasm_filter_opt_t fopt;
    std::memset(&fopt, 0, sizeof(fopt));
    fopt.threads = 4; fopt.k_mer_length = k; fopt.mz_win = w;
    fopt.is_hpc = 0; fopt.min_read_len = -1;

    hifiasm_filter_t* hf = hifiasm_build_filter(files, 1, &fopt);
    CHECK(hf != nullptr, "filter builds successfully");
    if (!hf) { std::remove(fasta.c_str()); return 1; }
    const void* raw_hf = hifiasm_filter_raw_for_test(hf);
    CHECK(raw_hf != nullptr, "filter exposes raw table");

    hifiasm_sketch_ctx_t* ctx = hifiasm_sketch_ctx_init();
    CHECK(ctx != nullptr, "ctx init");

    std::string bigrep;
    for (int i = 0; i < 12; ++i) bigrep += repeat; /* ~2400bp pure repeat */
    std::string read = random_seq(3000, 424242) + bigrep +
                       random_seq(3000, 424243);

    /* ---- unfiltered baseline (for the passthrough check) ---- */
    const hifiasm_minimizer_t* umz = nullptr; int un = 0;
    int rc = hifiasm_sketch_minimizers_ctx(ctx, read.c_str(), (int)read.size(),
                                           w, k, 0, &umz, &un);
    CHECK(rc == 0 && un > 0, "unfiltered sketch ok");
    std::vector<hifiasm_minimizer_t> uvec(umz, umz + un);

    /* ---- filtered (hf + subsampling), the overlap-parity path ---- */
    const hifiasm_minimizer_t* fmz = nullptr; int fn = 0;
    rc = hifiasm_sketch_minimizers_ctx_filtered(ctx, read.c_str(),
                                                (int)read.size(), w, k, 0,
                                                hf, sample_dist, &fmz, &fn);
    CHECK(rc == 0, "filtered sketch ok");
    CHECK(fn > 0, "filtered produces minimizers");
    std::vector<hifiasm_minimizer_t> fvec(fmz, fmz + fn);
    check_valid_sorted_dedup(fvec.data(), fn, k, read.size());

    /* ---- (A) PARITY with the overlap sketch: bridge == mz1_ha_sketch ---- */
    std::map<uint32_t,uint64_t> oracle =
        oracle_seeds(read, w, k, raw_hf, sample_dist);
    std::fprintf(stderr, "[info] unfiltered=%d filtered=%d oracle=%zu\n",
                 un, fn, oracle.size());
    CHECK(maps_equal(oracle, fvec.data(), fn),
          "filtered bridge matches mz1_ha_sketch seed set exactly");

    /* ---- (B) filtered differs from unfiltered (hf changed selection) ---- */
    bool differs = ((int)oracle.size() != un);
    if (!differs) {
        std::map<uint32_t,uint64_t> um;
        for (auto& e : uvec) um[e.pos] = e.hash;
        differs = (um != oracle);
    }
    CHECK(differs, "hf changes minimizer selection vs unfiltered");

    /* ---- (C) determinism ---- */
    const hifiasm_minimizer_t* fmz2 = nullptr; int fn2 = 0;
    rc = hifiasm_sketch_minimizers_ctx_filtered(ctx, read.c_str(),
                                                (int)read.size(), w, k, 0,
                                                hf, sample_dist, &fmz2, &fn2);
    CHECK(rc == 0 && fn2 == fn, "filtered count deterministic");
    bool same = (fn2 == fn);
    for (int i = 0; same && i < fn; ++i)
        if (fmz2[i].pos != fvec[i].pos || fmz2[i].hash != fvec[i].hash)
            same = false;
    CHECK(same, "filtered output deterministic");

    /* ---- (D) hf=NULL, sample_dist<=0 reproduces the unfiltered sketch ---- */
    const hifiasm_minimizer_t* nmz = nullptr; int nn = 0;
    rc = hifiasm_sketch_minimizers_ctx_filtered(ctx, read.c_str(),
                                                (int)read.size(), w, k, 0,
                                                /*hf*/ nullptr,
                                                /*sample_dist*/ 0, &nmz, &nn);
    CHECK(rc == 0 && nn == un, "hf=NULL,d<=0 count equals unfiltered");
    bool eq = (nn == un);
    for (int i = 0; eq && i < nn; ++i)
        if (nmz[i].pos != uvec[i].pos || nmz[i].hash != uvec[i].hash) eq = false;
    CHECK(eq, "hf=NULL,d<=0 output equals unfiltered");

    hifiasm_sketch_ctx_destroy(ctx);
    hifiasm_filter_destroy(hf);
    std::remove(fasta.c_str());

    if (g_failures == 0) {
        std::printf("ALL FILTER TESTS PASSED\n");
        return 0;
    }
    std::fprintf(stderr, "%d FILTER TEST(S) FAILED\n", g_failures);
    return 1;
}
