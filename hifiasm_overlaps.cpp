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
static int sketch_into_ctx(hifiasm_sketch_ctx_t *ctx,
                           const char *seq, int len, int w, int k, int is_hpc,
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

    /* hf=NULL           -> no frequency filter (cnt always 0, nothing dropped)
     * sample_dist=0     -> skips select_mz_h subsampling (needs sample_dist>w)
     * pt=NULL, dp_min_len=-1 -> skips refine_sketch
     * is_unique=0, dp_e=-1, ws ignored without subsampling.
     * So the result depends only on (seq, len, w, k, is_hpc). */
    mz1_ha_sketch(seq, len, w, k, /*rid*/ 0, is_hpc ? 1 : 0, &ctx->mz,
                  /*hf*/ NULL, /*sample_dist*/ 0, /*k_flag*/ NULL,
                  /*dbg_ct*/ NULL, /*pt*/ NULL, /*min_freq*/ -1,
                  /*dp_min_len*/ -1, /*dp_e*/ -1.0f, &ctx->mt, /*ws*/ w,
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
    return sketch_into_ctx(ctx, seq, len, w, k, is_hpc, out_mz, out_n);
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
    int rc = sketch_into_ctx(&ctx, seq, len, w, k, is_hpc, &mz, &n);
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
