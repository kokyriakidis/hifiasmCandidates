/*
 * hifiasm_overlaps.cpp
 *
 * File-based bridge implementation. Builds an argv equivalent to the CLI and
 * reuses hifiasm's own option processing (so all presets, coverage logic and
 * validation behave exactly like the command-line tool), then runs the
 * overlap-detection pipeline.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hifiasm_overlaps.h"
#include "CommandLines.h"
#include "Process_Read.h"
#include "htab.h"

/* Defined in candidates.cpp (C++ linkage): runs index build + candidate
 * detection and, unless asm_opt.dbg_ovec_cal is set, the alignment/filter
 * step. Writes the PAF using asm_opt.output_file_name as the prefix. */
int ha_detect_candidates(void);

/* Reset the process-global read store between calls so the bridge can be
 * invoked more than once in a single process. init_opt() memsets asm_opt, so
 * options are fully reset by CommandLine_process below; R_INF is reset here. */
static void reset_read_store(void)
{
    memset(&R_INF, 0, sizeof(R_INF));
}

/* -----------------------------------------------------------------------
 * Frequency-filter handle
 * -----------------------------------------------------------------------
 * ha_ft_gen() returns the raw yak_ft_t as a void*; ha_ft_cnt()/ha_ft_destroy()
 * take that same void*. We wrap it in an opaque struct so the public API type
 * is distinct and can grow later (e.g. carry the k/w it was built with for
 * validation) without changing the ABI. */
struct hifiasm_filter_s {
    void *raw; /* yak_ft_t* from ha_ft_gen(); owned by this handle */
    int   k;   /* k the table was built with (for reference/debugging) */
    int   w;   /* w the table was built with */
    int   is_hpc;
};

/* Unwrap for passing to mz1_ha_sketch (which takes const void*). NULL-safe. */
static const void *hifiasm_filter_raw(const hifiasm_filter_t *hf)
{
    return hf ? hf->raw : NULL;
}

int hifiasm_detect_overlaps(const char *const *read_files,
                            int n_read_files,
                            const char *output_prefix,
                            const hifiasm_ovlp_opt_t *opt,
                            char **out_paf_path)
{
    if (out_paf_path) *out_paf_path = NULL;
    if (read_files == NULL || n_read_files < 1 || output_prefix == NULL) {
        fprintf(stderr, "[hifiasm_detect_overlaps] invalid arguments\n");
        return 1;
    }

    /* -------- build an argv equivalent to the CLI invocation -------- */
    /* Fixed slots: prog, -o, <prefix>, -t, <threads> ; optional: --ont,
     * -k <k>, -w <w>, -D <maxdiff>, --dbg-ovec ; then the read files. */
    int    argc = 0;
    int    max_argc = 16 + n_read_files;
    char **argv = (char**)calloc(max_argc, sizeof(char*));

    char threads_buf[16], k_buf[16], w_buf[32], d_buf[32];
    int  threads = (opt && opt->threads > 0) ? opt->threads : 1;
    snprintf(threads_buf, sizeof(threads_buf), "%d", threads);

    argv[argc++] = (char*)"hifiasm";
    argv[argc++] = (char*)"-o";
    argv[argc++] = (char*)output_prefix;
    argv[argc++] = (char*)"-t";
    argv[argc++] = threads_buf;

    if (opt && opt->is_ont)          argv[argc++] = (char*)"--ont";
    if (opt && opt->k_mer_length > 0) {
        snprintf(k_buf, sizeof(k_buf), "%d", opt->k_mer_length);
        argv[argc++] = (char*)"-k"; argv[argc++] = k_buf;
    }
    if (opt && opt->mz_win > 0) {
        snprintf(w_buf, sizeof(w_buf), "%d", opt->mz_win);
        argv[argc++] = (char*)"-w"; argv[argc++] = w_buf;
    }
    if (opt && opt->max_ov_diff > 0.0) {
        snprintf(d_buf, sizeof(d_buf), "%g", opt->max_ov_diff);
        argv[argc++] = (char*)"--max-od-ec"; argv[argc++] = d_buf;
    }
    /* raw_candidates re-uses the --dbg-ovec switch, which the driver maps to
     * "emit pre-alignment candidates" (.candidates.paf). */
    if (opt && opt->raw_candidates) argv[argc++] = (char*)"--dbg-ovec";

    for (int i = 0; i < n_read_files; ++i)
        argv[argc++] = (char*)read_files[i];

    /* -------- run the pipeline -------- */
    reset_read_store();
    yak_reset_realtime();
    init_opt(&asm_opt);
    int ok = CommandLine_process(argc, argv, &asm_opt);
    if (!ok) {
        free(argv);
        fprintf(stderr, "[hifiasm_detect_overlaps] option processing failed\n");
        return 2;
    }

    int ret = ha_detect_candidates();

    int raw = asm_opt.dbg_ovec_cal;
    destory_opt(&asm_opt);
    free(argv);

    if (ret != 0) {
        fprintf(stderr, "[hifiasm_detect_overlaps] pipeline failed (%d)\n", ret);
        return ret;
    }

    /* -------- report the output path -------- */
    if (out_paf_path) {
        const char *suffix = raw ? ".candidates.paf" : ".ovlp.paf";
        size_t n = strlen(output_prefix) + strlen(suffix) + 1;
        char *p = (char*)malloc(n);
        snprintf(p, n, "%s%s", output_prefix, suffix);
        *out_paf_path = p;
    }
    return 0;
}

/*
 * Build a high-occurrence k-mer filter table over the given read files.
 *
 * Reuses the same argv/CommandLine_process shell as hifiasm_detect_overlaps so
 * all preset defaults, k/w handling and validation match the CLI, then forces
 * the k, w, HPC mode and read-length cutoff the caller asked for and calls
 * ha_ft_gen(read_from_store=0), which streams the read files, records read
 * lengths/names in R_INF, and builds the filter table.
 *
 * Only the filter table is kept; R_INF (lengths/names) is dropped before
 * returning, so nothing genome-scale outlives the call except the table
 * itself.
 */
hifiasm_filter_t *hifiasm_build_filter(const char *const *read_files,
                                       int n_read_files,
                                       const hifiasm_filter_opt_t *opt)
{
    if (read_files == NULL || n_read_files < 1) {
        fprintf(stderr, "[hifiasm_build_filter] invalid arguments\n");
        return NULL;
    }

    int k       = (opt && opt->k_mer_length > 0) ? opt->k_mer_length : 50;
    int w       = (opt && opt->mz_win       > 0) ? opt->mz_win       : 50;
    int is_hpc  = (opt && opt->is_hpc) ? 1 : 0;
    int threads = (opt && opt->threads > 0) ? opt->threads : 1;
    /* Default: keep every read (rl_cut < 0) so the counted set matches every
     * read the caller will later sketch. A caller that deliberately mirrors
     * hifiasm's own length cutoff can request it via min_read_len >= 0. */
    int64_t rl_cut = (opt && opt->min_read_len >= 0) ? opt->min_read_len : -1;

    /* -------- build an argv equivalent to a CLI invocation -------- */
    int    argc = 0;
    int    max_argc = 12 + n_read_files;
    char **argv = (char**)calloc(max_argc, sizeof(char*));
    if (argv == NULL) {
        fprintf(stderr, "[hifiasm_build_filter] out of memory\n");
        return NULL;
    }

    char threads_buf[16], k_buf[16], w_buf[16];
    snprintf(threads_buf, sizeof(threads_buf), "%d", threads);
    snprintf(k_buf, sizeof(k_buf), "%d", k);
    snprintf(w_buf, sizeof(w_buf), "%d", w);

    argv[argc++] = (char*)"hifiasm";
    argv[argc++] = (char*)"-o";
    argv[argc++] = (char*)"hifiasm_filter"; /* prefix; no files are written */
    argv[argc++] = (char*)"-t";
    argv[argc++] = threads_buf;
    argv[argc++] = (char*)"-k"; argv[argc++] = k_buf;
    argv[argc++] = (char*)"-w"; argv[argc++] = w_buf;
    for (int i = 0; i < n_read_files; ++i)
        argv[argc++] = (char*)read_files[i];

    /* -------- run option processing + table build -------- */
    reset_read_store();
    yak_reset_realtime();
    init_opt(&asm_opt);
    int ok = CommandLine_process(argc, argv, &asm_opt);
    if (!ok) {
        free(argv);
        fprintf(stderr, "[hifiasm_build_filter] option processing failed\n");
        return NULL;
    }

    /* Force the exact sketch parameters. CommandLine_process may have applied
     * preset k/w; override so the filter's k/w/HPC match what we will sketch
     * with. HA_F_NO_HPC controls HPC mode inside ha_ft_gen via ha_count. */
    asm_opt.k_mer_length = k;
    asm_opt.mz_win       = w;
    if (is_hpc) asm_opt.flag &= ~HA_F_NO_HPC;
    else        asm_opt.flag |=  HA_F_NO_HPC;
    asm_opt.rl_cut = rl_cut;

    int hom_cov = -1;
    /* read_from_store=0: stream the files, writing read lengths (not
     * sequences) into R_INF and building the high-occurrence table. */
    void *raw = ha_ft_gen(&asm_opt, &R_INF, &hom_cov, /*is_hp_mode*/ 0,
                          /*read_from_store*/ 0);

    /* Drop the transient store. ha_ft_gen with read_from_store=0 uses
     * HAF_RS_WRITE_LEN, so init_All_reads ran (allocating read_length and
     * name_index) but malloc_All_reads did NOT (N_site, read_sperate, paf,
     * name, ... are still NULL). destory_All_reads would dereference those
     * NULL per-read arrays while total_reads>0, so free only what was actually
     * allocated here, then zero the store. */
    free(R_INF.read_length);
    free(R_INF.name_index);
    reset_read_store();
    destory_opt(&asm_opt);
    free(argv);

    if (raw == NULL) {
        fprintf(stderr, "[hifiasm_build_filter] filter table build failed\n");
        return NULL;
    }

    hifiasm_filter_t *hf = (hifiasm_filter_t*)calloc(1, sizeof(*hf));
    if (hf == NULL) {
        ha_ft_destroy(raw);
        fprintf(stderr, "[hifiasm_build_filter] out of memory\n");
        return NULL;
    }
    hf->raw = raw; hf->k = k; hf->w = w; hf->is_hpc = is_hpc;
    return hf;
}

void hifiasm_filter_destroy(hifiasm_filter_t *hf)
{
    if (hf == NULL) return;
    ha_ft_destroy(hf->raw);
    free(hf);
}

int32_t hifiasm_filter_count(const hifiasm_filter_t *hf, uint64_t hash)
{
    if (hf == NULL || hf->raw == NULL) return 0;
    return ha_ft_cnt(hf->raw, hash);
}

const void *hifiasm_filter_raw_for_test(const hifiasm_filter_t *hf)
{
    return hifiasm_filter_raw(hf);
}

/* Order minimizers by ascending START position, breaking ties by hash so the
 * output is deterministic regardless of the sketcher's emission order. */
static int cmp_minimizer_by_pos(const void *pa, const void *pb)
{
    const hifiasm_minimizer_t *a = (const hifiasm_minimizer_t*)pa;
    const hifiasm_minimizer_t *b = (const hifiasm_minimizer_t*)pb;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    if (a->hash != b->hash) return a->hash < b->hash ? -1 : 1;
    return 0;
}

/*
 * Reusable sketch context. Keeps the sketcher scratch (mz, mt) and the
 * converted-output buffer (out) alive across calls so steady-state sketching
 * does no per-read allocation once the buffers have grown to fit.
 */
struct hifiasm_sketch_ctx_s {
    ha_mz1_v            mz;    /* raw sketch output (END positions, hashed key) */
    st_mt_t             mt;    /* mz1_ha_sketch scratch (dereferenced always)   */
    hifiasm_minimizer_t *out;  /* converted, START-anchored, sorted+deduped     */
    size_t              out_m; /* capacity of `out` in elements                 */
};

hifiasm_sketch_ctx_t *hifiasm_sketch_ctx_init(void)
{
    hifiasm_sketch_ctx_t *ctx =
        (hifiasm_sketch_ctx_t*)calloc(1, sizeof(hifiasm_sketch_ctx_t));
    return ctx; /* NULL on failure is fine; callers check. */
}

void hifiasm_sketch_ctx_destroy(hifiasm_sketch_ctx_t *ctx)
{
    if (ctx == NULL) return;
    /* mz.a / mt.a were grown by mz1_ha_sketch with km=NULL -> stdlib malloc. */
    free(ctx->mz.a);
    free(ctx->mt.a);
    free(ctx->out);
    free(ctx);
}

/*
 * Core sketch into ctx-owned storage. On success sets *out_mz (pointer into
 * ctx->out, owned by the ctx) and *out_n. Shared by the context and one-shot
 * public entry points.
 */
/* hifiasm's overlap sketch passes asm_opt.mz_rewin as the select_mz_h window
 * (see anchor.cpp; sketch.cpp threads it in as the sketch `ws` argument). Its
 * default is 1000. We reproduce that here so the distance-subsampling stage
 * matches the overlap path exactly. Only consulted when subsampling is active
 * (sample_dist > w); ignored otherwise. */
#define HIFIASM_MZ_REWIN_DEFAULT 1000

static int sketch_into_ctx(hifiasm_sketch_ctx_t *ctx,
                           const char *seq, int len, int w, int k, int is_hpc,
                           const void *hf, int sample_dist, int ws,
                           const hifiasm_minimizer_t **out_mz, int *out_n)
{
    *out_mz = NULL;
    *out_n  = 0;
    if (seq == NULL) {
        fprintf(stderr, "[hifiasm_sketch_minimizers] invalid arguments\n");
        return 1;
    }
    /* Mirror mz1_ha_sketch's own asserts as graceful errors so callers get a
     * return code instead of an abort on out-of-range parameters. */
    if (len <= 0 || w <= 0 || w >= 256 || k <= 0 || k > 63) {
        fprintf(stderr, "[hifiasm_sketch_minimizers] parameter out of range "
                        "(len=%d, w=%d, k=%d)\n", len, w, k);
        return 1;
    }

    /* Reuse buffers: reset logical length, keep capacity/pointer. mz1_ha_sketch
     * appends via kv_push_km and grows only when n exceeds m. */
    ctx->mz.n = 0;
    ctx->mt.n = 0;

    /* Filter parameters (see hifiasm_sketch_minimizers_ctx_filtered):
     *   hf          -> frequency filter; NULL means no filtering (cnt always 0).
     *   sample_dist -> distance subsampling via select_mz_h, active only when
     *                  sample_dist > w. <=0 disables it.
     *   ws          -> select_mz_h window; hifiasm's overlap path uses
     *                  asm_opt.mz_rewin (default 1000). Only consulted when
     *                  subsampling is active.
     * pt=NULL, dp_min_len=-1 -> refine_sketch never runs (matches the overlap
     * sketch, which also passes pt=NULL). is_unique=0, dp_e=-1.
     * With hf=NULL and sample_dist<=0 the result depends only on
     * (seq, len, w, k, is_hpc), i.e. the raw unfiltered sketch. */
    mz1_ha_sketch(seq, len, w, k, /*rid*/ 0, is_hpc ? 1 : 0, &ctx->mz,
                  hf, sample_dist, /*k_flag*/ NULL,
                  /*dbg_ct*/ NULL, /*pt*/ NULL, /*min_freq*/ -1,
                  /*dp_min_len*/ -1, /*dp_e*/ -1.0f, &ctx->mt, ws,
                  /*is_unique*/ 0, /*km*/ NULL);

    const uint32_t n = ctx->mz.n;
    if (n == 0) return 0;

    /* Grow the output buffer if needed (keeps capacity across calls). */
    if ((size_t)n > ctx->out_m) {
        size_t m = ctx->out_m ? ctx->out_m : 64;
        while (m < (size_t)n) m <<= 1;
        hifiasm_minimizer_t *p =
            (hifiasm_minimizer_t*)realloc(ctx->out, m * sizeof(*p));
        if (p == NULL) {
            fprintf(stderr, "[hifiasm_sketch_minimizers] out of memory\n");
            return 1;
        }
        ctx->out = p;
        ctx->out_m = m;
    }

    /* ha_mz1_t stores the END position; convert to START = end + 1 - span.
     * With is_hpc=0, span == k. mz1_ha_sketch never emits a k-mer whose span
     * runs past the sequence, so start >= 0 and start + span <= len hold. */
    for (uint32_t i = 0; i < n; ++i) {
        const ha_mz1_t *m = &ctx->mz.a[i];
        const uint32_t span = (uint32_t)m->span;
        const uint32_t end  = (uint32_t)m->pos;
        ctx->out[i].pos  = end + 1u - span;
        ctx->out[i].span = span;
        ctx->out[i].rev  = (uint32_t)m->rev;
        ctx->out[i].hash = m->x;
    }

    /* Emitted order follows window scans and is largely position-ordered but
     * not strictly (identical-minimizer flushes and end handling can interleave
     * it). Sort by START for a clean monotonic list. */
    qsort(ctx->out, n, sizeof(hifiasm_minimizer_t), cmp_minimizer_by_pos);

    /* Deduplicate by START position. Exactly one k-mer starts at any position,
     * so equal positions are fully identical entries; collapsing them lets the
     * caller skip its own dedup pass. */
    uint32_t u = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (u == 0 || ctx->out[i].pos != ctx->out[u - 1].pos) {
            ctx->out[u++] = ctx->out[i];
        }
    }

    *out_mz = ctx->out;
    *out_n  = (int)u;
    return 0;
}

int hifiasm_sketch_minimizers_ctx(hifiasm_sketch_ctx_t *ctx,
                                  const char *seq,
                                  int len,
                                  int w,
                                  int k,
                                  int is_hpc,
                                  const hifiasm_minimizer_t **out_mz,
                                  int *out_n)
{
    if (out_mz) *out_mz = NULL;
    if (out_n)  *out_n  = 0;
    if (ctx == NULL || out_mz == NULL || out_n == NULL) {
        fprintf(stderr, "[hifiasm_sketch_minimizers_ctx] invalid arguments\n");
        return 1;
    }
    /* Unfiltered: hf=NULL, sample_dist=0 (subsampling off), ws unused. */
    return sketch_into_ctx(ctx, seq, len, w, k, is_hpc,
                           /*hf*/ NULL, /*sample_dist*/ 0, /*ws*/ w,
                           out_mz, out_n);
}

int hifiasm_sketch_minimizers_ctx_filtered(hifiasm_sketch_ctx_t *ctx,
                                           const char *seq,
                                           int len,
                                           int w,
                                           int k,
                                           int is_hpc,
                                           const hifiasm_filter_t *hf,
                                           int sample_dist,
                                           const hifiasm_minimizer_t **out_mz,
                                           int *out_n)
{
    if (out_mz) *out_mz = NULL;
    if (out_n)  *out_n  = 0;
    if (ctx == NULL || out_mz == NULL || out_n == NULL) {
        fprintf(stderr,
                "[hifiasm_sketch_minimizers_ctx_filtered] invalid arguments\n");
        return 1;
    }
    /* hifiasm_filter_t wraps the raw yak_ft_t that mz1_ha_sketch expects as a
     * const void*; unwrap it here. NULL hf is allowed and disables frequency
     * filtering (subsampling may still apply). */
    return sketch_into_ctx(ctx, seq, len, w, k, is_hpc,
                           hifiasm_filter_raw(hf), sample_dist,
                           /*ws*/ HIFIASM_MZ_REWIN_DEFAULT, out_mz, out_n);
}

int hifiasm_sketch_minimizers(const char *seq,
                              int len,
                              int w,
                              int k,
                              int is_hpc,
                              hifiasm_minimizer_t **out_mz,
                              int *out_n)
{
    if (out_mz) *out_mz = NULL;
    if (out_n)  *out_n  = 0;
    if (out_mz == NULL || out_n == NULL) {
        fprintf(stderr, "[hifiasm_sketch_minimizers] invalid arguments\n");
        return 1;
    }

    /* One-shot: sketch into a temporary context, then hand the caller an
     * owned copy sized exactly to the result so they can free() it. */
    hifiasm_sketch_ctx_t ctx = {{0, 0, NULL}, {0, 0, NULL}, NULL, 0};
    const hifiasm_minimizer_t *mz = NULL;
    int n = 0;
    int rc = sketch_into_ctx(&ctx, seq, len, w, k, is_hpc,
                             /*hf*/ NULL, /*sample_dist*/ 0, /*ws*/ w, &mz, &n);
    if (rc != 0 || n == 0) {
        free(ctx.mz.a); free(ctx.mt.a); free(ctx.out);
        return rc;
    }

    hifiasm_minimizer_t *result =
        (hifiasm_minimizer_t*)malloc((size_t)n * sizeof(hifiasm_minimizer_t));
    if (result == NULL) {
        free(ctx.mz.a); free(ctx.mt.a); free(ctx.out);
        fprintf(stderr, "[hifiasm_sketch_minimizers] out of memory\n");
        return 1;
    }
    memcpy(result, mz, (size_t)n * sizeof(hifiasm_minimizer_t));
    free(ctx.mz.a); free(ctx.mt.a); free(ctx.out);

    *out_mz = result;
    *out_n  = n;
    return 0;
}
