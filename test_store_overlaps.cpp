/*
 * Parity test for the shared-read-store bridge path.
 *
 * The store path (hifiasm_reads_store_load + *_from_store) must produce the
 * SAME overlaps and read-name table as the file-based in-memory path
 * (hifiasm_detect_overlaps_mem), because it feeds hifiasm's pipeline from an
 * R_INF that is built to be bit-identical to what the file loader produces.
 *
 * The test:
 *   1. Generates a small synthetic read set with real overlaps (a reference is
 *      cut into reads with staggered, overlapping windows), plus a couple of
 *      reads containing ambiguous (N) bases to exercise the N_site path.
 *   2. Writes them to a temp FASTA and runs hifiasm_detect_overlaps_mem() over
 *      the file -> the oracle overlap set + name table.
 *   3. Loads the identical reads from memory (hifiasm_reads_store_load), builds
 *      a filter from the store, and runs hifiasm_detect_overlaps_from_store().
 *   4. Asserts the name tables match and the overlap sets are equal as sets.
 *
 * Build (from the submodule root, after `make lib`):
 *   g++ -std=c++11 -O2 -I. test_store_overlaps.cpp \
 *       libhifiasm_overlaps.a -lz -lpthread -lm -o test_store_overlaps
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <unistd.h>

#include "hifiasm_overlaps.h"

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

struct Read {
    std::string name;
    std::string seq;
};

/* Build overlapping reads by sliding a window across a long reference. Adjacent
 * reads share `ref_len_per - step` bases, guaranteeing detectable overlaps. */
static std::vector<Read> make_reads()
{
    std::vector<Read> reads;
    const std::string ref = random_seq(30000, 12345);
    const size_t win = 6000;
    const size_t step = 3000;
    int idx = 0;
    for (size_t s = 0; s + win <= ref.size(); s += step) {
        Read r;
        r.name = "read" + std::to_string(idx++);
        r.seq = ref.substr(s, win);
        reads.push_back(r);
    }
    /* Inject a handful of ambiguous bases into one read to exercise N_site. */
    if (!reads.empty()) {
        std::string &q = reads[1].seq;
        for (size_t i = 100; i < 110 && i < q.size(); ++i) q[i] = 'N';
    }
    return reads;
}

static bool write_fasta(const std::vector<Read> &reads, const char *path)
{
    FILE *fp = std::fopen(path, "w");
    if (!fp) return false;
    for (const auto &r : reads)
        std::fprintf(fp, ">%s\n%s\n", r.name.c_str(), r.seq.c_str());
    std::fclose(fp);
    return true;
}

/* A comparable key for an overlap, resolving q_id/t_id to names so the two
 * paths can be compared regardless of internal read ordering. */
struct OvKey {
    std::string q, t;
    uint32_t qs, qe, ts, te, nm, bl;
    uint8_t ss;
    bool operator<(const OvKey &o) const {
        if (q != o.q) return q < o.q;
        if (t != o.t) return t < o.t;
        if (qs != o.qs) return qs < o.qs;
        if (qe != o.qe) return qe < o.qe;
        if (ts != o.ts) return ts < o.ts;
        if (te != o.te) return te < o.te;
        if (nm != o.nm) return nm < o.nm;
        if (bl != o.bl) return bl < o.bl;
        return ss < o.ss;
    }
    bool operator==(const OvKey &o) const {
        return q == o.q && t == o.t && qs == o.qs && qe == o.qe &&
               ts == o.ts && te == o.te && nm == o.nm && bl == o.bl &&
               ss == o.ss;
    }
};

static std::string name_of(const char *names, const uint64_t *off, uint64_t i)
{
    return std::string(names + off[i], names + off[i + 1]);
}

static std::set<OvKey> to_keyset(const hifiasm_overlap_t *ov, uint64_t n,
                                 const char *names, const uint64_t *off)
{
    std::set<OvKey> s;
    for (uint64_t i = 0; i < n; ++i) {
        OvKey k;
        k.q = name_of(names, off, ov[i].q_id);
        k.t = name_of(names, off, ov[i].t_id);
        k.qs = ov[i].q_start; k.qe = ov[i].q_end;
        k.ts = ov[i].t_start; k.te = ov[i].t_end;
        k.nm = ov[i].n_match; k.bl = ov[i].block_len;
        k.ss = ov[i].is_same_strand;
        s.insert(k);
    }
    return s;
}

int main(void)
{
    std::vector<Read> reads = make_reads();
    CHECK(reads.size() >= 4, "need several reads for overlaps");

    char tmpl[] = "/tmp/hifiasm_store_test_XXXXXX.fa";
    /* mkstemps keeps the .fa suffix (kseq sniffs by content, but be tidy). */
    int fd = mkstemps(tmpl, 3);
    CHECK(fd >= 0, "mkstemps");
    if (fd >= 0) close(fd);
    CHECK(write_fasta(reads, tmpl), "write fasta");

    hifiasm_ovlp_opt_t opt;
    std::memset(&opt, 0, sizeof(opt));
    opt.threads = 1;
    /* Raw candidate set: deterministic and independent of the alignment
     * filter, so the two paths are directly comparable. */
    opt.raw_candidates = 1;

    /* ---- oracle: file-based in-memory path ---- */
    const char *files[1] = { tmpl };
    hifiasm_overlap_t *fov = NULL; uint64_t fn = 0;
    char *fnames = NULL; uint64_t *foff = NULL; uint64_t fnr = 0;
    int rc = hifiasm_detect_overlaps_mem(files, 1, &opt,
                                         &fov, &fn, &fnames, &foff, &fnr,
                                         NULL, NULL);
    CHECK(rc == 0, "file mem path rc==0");
    CHECK(fnr == reads.size(), "file path read count matches input");

    /* ---- store path ---- */
    std::vector<hifiasm_read_t> hr(reads.size());
    for (size_t i = 0; i < reads.size(); ++i) {
        hr[i].seq = reads[i].seq.data();
        hr[i].seq_len = reads[i].seq.size();
        hr[i].name = reads[i].name.data();
        hr[i].name_len = (uint32_t)reads[i].name.size();
    }
    rc = hifiasm_reads_store_load(hr.data(), hr.size());
    CHECK(rc == 0, "store load rc==0");

    /* Build a filter from the store (HiFi overlap params: HPC, k=w=51). This
     * exercises hifiasm_build_filter_from_store; the overlap detection below
     * builds its own filter internally, so this is an independent smoke test. */
    hifiasm_filter_opt_t fopt;
    std::memset(&fopt, 0, sizeof(fopt));
    fopt.threads = 1; fopt.k_mer_length = 51; fopt.mz_win = 51; fopt.is_hpc = 1;
    fopt.min_read_len = -1;
    hifiasm_filter_t *hf = hifiasm_build_filter_from_store(&fopt);
    CHECK(hf != NULL, "build_filter_from_store");
    if (hf) hifiasm_filter_destroy(hf);

    hifiasm_overlap_t *sov = NULL; uint64_t sn = 0;
    char *snames = NULL; uint64_t *soff = NULL; uint64_t snr = 0;
    rc = hifiasm_detect_overlaps_from_store(&opt,
                                            &sov, &sn, &snames, &soff, &snr,
                                            NULL, NULL, NULL, NULL);
    CHECK(rc == 0, "store overlap path rc==0");
    CHECK(snr == fnr, "store path read count matches file path");

    /* Names must match position-for-position: the store preserves input order,
     * and the file loader loads reads in file order, which is the same. */
    bool names_ok = (snr == fnr);
    for (uint64_t i = 0; names_ok && i < snr; ++i)
        names_ok = (name_of(fnames, foff, i) == name_of(snames, soff, i));
    CHECK(names_ok, "name tables identical between paths");

    /* Overlap sets must be identical (compared by name-resolved key). */
    std::set<OvKey> fs = to_keyset(fov, fn, fnames, foff);
    std::set<OvKey> ss = to_keyset(sov, sn, snames, soff);
    CHECK(fs == ss, "overlap sets identical between file and store paths");
    std::fprintf(stderr,
                 "[info] file overlaps=%llu store overlaps=%llu (unique f=%zu s=%zu)\n",
                 (unsigned long long)fn, (unsigned long long)sn,
                 fs.size(), ss.size());

    /* Tear down the store and its buffers before running any further file-mem
     * pipeline. The file-mem path drives the same global R_INF and frees it
     * internally (candidates.cpp: `if (!from_store) destory_All_reads`), so a
     * loaded store must be released first or the trailing store_release would
     * double-free R_INF. */
    hifiasm_overlaps_mem_free(sov, snames, soff, NULL);
    sov = NULL; snames = NULL; soff = NULL;
    hifiasm_reads_store_release();

    /* ---- CIGAR passthrough: alignment (non-raw) path ---- */
    /* The raw candidate set carries no base-level CIGAR, so re-run the file
     * path WITHOUT raw_candidates to exercise gen_hc_r_alin_ea and assert that
     * every returned overlap carries a CIGAR whose op spans match its
     * coordinates, in the SAM/PAF convention (op2='I' consumes query, op3='D'
     * consumes target). This is what dinara consumes with no op remapping. */
    {
        hifiasm_ovlp_opt_t aopt;
        std::memset(&aopt, 0, sizeof(aopt));
        aopt.threads = 1;
        aopt.raw_candidates = 0;   /* run the base-level alignment/filter */

        hifiasm_overlap_t *aov = NULL; uint64_t an = 0;
        char *anames = NULL; uint64_t *aoff = NULL; uint64_t anr = 0;
        uint16_t *acig = NULL; uint64_t acig_n = 0;
        int arc = hifiasm_detect_overlaps_mem(files, 1, &aopt,
                                              &aov, &an, &anames, &aoff, &anr,
                                              &acig, &acig_n);
        CHECK(arc == 0, "alignment path rc==0");

        uint64_t withCigar = 0, spanOk = 0, spanBad = 0;
        uint64_t maxTokEnd = 0;
        for (uint64_t i = 0; i < an; ++i) {
            if (aov[i].cigar_len == 0) continue;
            ++withCigar;
            uint64_t end = aov[i].cigar_offset + aov[i].cigar_len;
            if (end > maxTokEnd) maxTokEnd = end;
            CHECK(end <= acig_n, "cigar slice within arena");

            /* Decode tokens: op = tok>>14, len = tok & 0x3FFF. Sum the query
             * and target spans and compare to the overlap coordinates. */
            uint64_t qSpan = 0, tSpan = 0;
            for (uint32_t j = 0; j < aov[i].cigar_len; ++j) {
                uint16_t tok = acig[aov[i].cigar_offset + j];
                uint8_t  op  = (uint8_t)(tok >> 14);
                uint32_t len = (uint32_t)(tok & 0x3FFF);
                switch (op) {
                    case 0: case 1: qSpan += len; tSpan += len; break; /* =,X */
                    case 2: qSpan += len; break; /* I: query only */
                    case 3: tSpan += len; break; /* D: target only */
                }
            }
            uint64_t qCoord = aov[i].q_end - aov[i].q_start;
            uint64_t tCoord = aov[i].t_end - aov[i].t_start;
            if (qSpan == qCoord && tSpan == tCoord) ++spanOk; else ++spanBad;
        }
        CHECK(withCigar > 0, "alignment path returned at least one CIGAR");
        CHECK(spanBad == 0, "every CIGAR span matches its overlap coordinates");
        CHECK(acig == NULL || maxTokEnd <= acig_n, "arena length consistent");
        std::fprintf(stderr,
                     "[info] alignment overlaps=%llu withCigar=%llu spanOk=%llu "
                     "spanBad=%llu tokens=%llu\n",
                     (unsigned long long)an, (unsigned long long)withCigar,
                     (unsigned long long)spanOk, (unsigned long long)spanBad,
                     (unsigned long long)acig_n);

        hifiasm_overlaps_mem_free(aov, anames, aoff, acig);
    }

    /* Release everything. The store was already released above (before the
     * CIGAR/file-mem block); calling release() again must be a safe no-op. */
    hifiasm_overlaps_mem_free(fov, fnames, foff, NULL);
    hifiasm_reads_store_release();

    std::remove(tmpl);

    if (g_failures) {
        std::fprintf(stderr, "%d CHECK(s) FAILED\n", g_failures);
        return 1;
    }
    std::fprintf(stderr, "all store-parity checks passed\n");
    return 0;
}
