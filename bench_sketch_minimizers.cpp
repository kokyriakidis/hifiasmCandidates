/*
 * Benchmark for hifiasm_sketch_minimizers() (no-HPC minimizer position source).
 *
 * Reports, for a synthetic read set at (k, w):
 *   - throughput   : bases/second and reads/second of sketching
 *   - density      : minimizers per base (expected ~2/(w+1) for a random
 *                    minimizer scheme; deviations indicate the selection's
 *                    actual sampling rate)
 *   - gap spacing  : mean / p50 / p95 / max distance between consecutive
 *                    minimizer START positions (seed spread => chaining reach)
 *
 * These are the hifiasm-side numbers for the head-to-head; the dinara side
 * (simd-minimizers on the same reads via Kmers.useHifiasmMinimizers=false vs
 * =true) is measured inside dinara's own build where both backends link.
 *
 * Build (after `make lib`):
 *   g++ -std=c++11 -O3 -msse4.2 -mpopcnt -I. bench_sketch_minimizers.cpp \
 *       libhifiasm_overlaps.a -lz -lpthread -lm -o bench_sketch_minimizers
 * Run:
 *   ./bench_sketch_minimizers [num_reads] [read_len] [k] [w]
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

#include "hifiasm_overlaps.h"

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

int main(int argc, char** argv)
{
    const int   num_reads = argc > 1 ? atoi(argv[1]) : 2000;
    const size_t read_len = argc > 2 ? (size_t)atol(argv[2]) : 15000;
    const int   k         = argc > 3 ? atoi(argv[3]) : 50;
    const int   w         = argc > 4 ? atoi(argv[4]) : 50;

    // Pre-generate reads so generation time is excluded from the measurement.
    std::vector<std::string> reads;
    reads.reserve(num_reads);
    for (int i = 0; i < num_reads; ++i)
        reads.push_back(random_seq(read_len, (uint32_t)(i * 2654435761u + 1)));

    uint64_t total_bases = 0, total_mz = 0;
    std::vector<uint32_t> gaps;
    gaps.reserve((size_t)num_reads * (2 * read_len / (w + 1) + 4));

    /* --- One-shot path: allocates + frees per read. --- */
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < num_reads; ++i) {
        hifiasm_minimizer_t* mz = nullptr;
        int n = 0;
        int rc = hifiasm_sketch_minimizers(reads[i].c_str(), (int)reads[i].size(),
                                           w, k, /*is_hpc*/ 0, &mz, &n);
        if (rc != 0) { std::fprintf(stderr, "sketch failed on read %d\n", i); return 1; }
        total_bases += reads[i].size();
        total_mz += (uint64_t)n;
        for (int j = 1; j < n; ++j) gaps.push_back(mz[j].pos - mz[j - 1].pos);
        free(mz);
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double secs =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();

    /* --- Context path: buffers reused, no per-read allocation at steady state. --- */
    hifiasm_sketch_ctx_t* ctx = hifiasm_sketch_ctx_init();
    uint64_t ctx_total_mz = 0;
    const auto c0 = std::chrono::steady_clock::now();
    for (int i = 0; i < num_reads; ++i) {
        const hifiasm_minimizer_t* mz = nullptr;
        int n = 0;
        int rc = hifiasm_sketch_minimizers_ctx(ctx, reads[i].c_str(),
                                               (int)reads[i].size(), w, k, 0, &mz, &n);
        if (rc != 0) { std::fprintf(stderr, "ctx sketch failed on read %d\n", i); return 1; }
        ctx_total_mz += (uint64_t)n;
    }
    const auto c1 = std::chrono::steady_clock::now();
    const double ctx_secs =
        std::chrono::duration_cast<std::chrono::duration<double>>(c1 - c0).count();
    hifiasm_sketch_ctx_destroy(ctx);

    std::sort(gaps.begin(), gaps.end());
    auto pct = [&](double p) -> uint32_t {
        if (gaps.empty()) return 0;
        size_t idx = (size_t)(p * (gaps.size() - 1));
        return gaps[idx];
    };
    double mean_gap = 0.0;
    for (uint32_t g : gaps) mean_gap += g;
    if (!gaps.empty()) mean_gap /= (double)gaps.size();

    std::printf("=== hifiasm no-HPC sketch benchmark ===\n");
    std::printf("reads            : %d x %zu bp  (%.1f Mbp total)\n",
                num_reads, read_len, total_bases / 1e6);
    std::printf("k, w             : %d, %d\n", k, w);
    std::printf("one-shot time    : %.3f s  (%.1f Mbp/s, %.0f reads/s)\n",
                secs, (total_bases / 1e6) / secs, num_reads / secs);
    std::printf("ctx time         : %.3f s  (%.1f Mbp/s, %.0f reads/s)\n",
                ctx_secs, (total_bases / 1e6) / ctx_secs, num_reads / ctx_secs);
    std::printf("ctx speedup      : %.2fx  (same minimizers: %s)\n",
                secs / ctx_secs, ctx_total_mz == total_mz ? "yes" : "NO");
    std::printf("minimizers       : %llu total\n",
                (unsigned long long)total_mz);
    std::printf("density (mz/base): %.4f  (ideal 2/(w+1) = %.4f)\n",
                (double)total_mz / (double)total_bases, 2.0 / (w + 1));
    std::printf("gap mean/p50/p95/max: %.1f / %u / %u / %u bp\n",
                mean_gap, pct(0.50), pct(0.95),
                gaps.empty() ? 0u : gaps.back());
    return 0;
}
