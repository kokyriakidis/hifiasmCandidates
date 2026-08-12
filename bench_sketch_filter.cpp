/*
 * Benchmark for the overlap-parity marker filter.
 *
 * Builds a repeat-containing genome, tiles it into reads, builds the no-HPC
 * hifiasm filter over those reads, then sketches every read three ways and
 * reports total minimizer counts and wall-clock time:
 *   - unfiltered            (hifiasm_sketch_minimizers_ctx)
 *   - filtered, no subsample (hf only)
 *   - filtered + subsample   (hf + sample_dist, the overlap-parity path)
 *
 * This quantifies how the filter reshapes the marker set (repeat seeds replaced
 * by unique ones; subsampling thins density) and the per-read cost of the extra
 * filtering, which is the relevant signal for dinara's marker phase.
 *
 * Build (from the submodule root, after `make lib`):
 *   g++ -std=c++11 -O2 -I. bench_sketch_filter.cpp \
 *       libhifiasm_overlaps.a -lz -lpthread -lm -o bench_sketch_filter
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <chrono>
#include <unistd.h>

#include "hifiasm_overlaps.h"

static std::string random_seq(size_t n, uint32_t seed)
{
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    std::string s; s.resize(n);
    uint32_t x = seed ? seed : 1;
    for (size_t i = 0; i < n; ++i) {
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        s[i] = bases[x & 3];
    }
    return s;
}

int main(void)
{
    const int k = 50, w = 50, sample_dist = 500;

    /* Genome: unique backbone + a high-copy short repeat. */
    std::string repeat = random_seq(200, 0xABCD);
    std::string genome;
    for (int i = 0; i < 60; ++i) {
        genome += random_seq(10000, 7000u + (uint32_t)i);
        genome += repeat;
    }

    /* Tile into reads; keep them in memory for sketching, and also write them
     * to a temp FASTA so the filter can be built over the same set. */
    const int win = 6000, step = 500;
    std::vector<std::string> reads;
    char tmpl[] = "/tmp/hifiasm_flt_bench_XXXXXX";
    int fd = mkstemp(tmpl);
    FILE* fp = fdopen(fd, "w");
    int rid = 0;
    for (size_t p = 0; p + win <= genome.size(); p += step) {
        std::string r = genome.substr(p, win);
        std::fprintf(fp, ">r%d\n%s\n", rid++, r.c_str());
        reads.push_back(std::move(r));
    }
    std::fclose(fp);
    std::fprintf(stderr, "genome=%zu bp, reads=%zu (win=%d step=%d)\n",
                 genome.size(), reads.size(), win, step);

    const char* files[1] = { tmpl };
    hifiasm_filter_opt_t fopt; std::memset(&fopt, 0, sizeof(fopt));
    fopt.threads = 4; fopt.k_mer_length = k; fopt.mz_win = w;
    fopt.is_hpc = 0; fopt.min_read_len = -1;
    hifiasm_filter_t* hf = hifiasm_build_filter(files, 1, &fopt);
    if (!hf) { std::fprintf(stderr, "filter build failed\n"); std::remove(tmpl); return 1; }

    hifiasm_sketch_ctx_t* ctx = hifiasm_sketch_ctx_init();

    struct Run { const char* label; bool use_hf; int sd; };
    Run runs[3] = {
        { "unfiltered           ", false, 0 },
        { "filtered (hf only)   ", true,  0 },
        { "filtered (hf+sample) ", true,  sample_dist },
    };

    for (const Run& run : runs) {
        long long total = 0;
        auto t0 = std::chrono::steady_clock::now();
        for (const std::string& r : reads) {
            const hifiasm_minimizer_t* mz = nullptr; int n = 0;
            if (run.use_hf) {
                hifiasm_sketch_minimizers_ctx_filtered(
                    ctx, r.c_str(), (int)r.size(), w, k, 0, hf, run.sd, &mz, &n);
            } else {
                hifiasm_sketch_minimizers_ctx(
                    ctx, r.c_str(), (int)r.size(), w, k, 0, &mz, &n);
            }
            total += n;
        }
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
        std::printf("%s markers=%lld  time=%.2f ms  (%.3f us/read)\n",
                    run.label, total, ms, 1000.0 * ms / (double)reads.size());
    }

    hifiasm_sketch_ctx_destroy(ctx);
    hifiasm_filter_destroy(hf);
    std::remove(tmpl);
    return 0;
}
