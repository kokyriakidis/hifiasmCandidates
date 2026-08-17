// candidates.cpp
//
// Driver for INITIAL CANDIDATE OVERLAP DETECTION (pre error-correction).
//
// Two modes:
//
//  * default: build the k-mer filter table + position index, run candidate
//    detection (minimizer sketching + anchoring + chaining, h_ec_lchain) and
//    then the post-chaining base-level alignment/filter (gen_hc_r_alin_ea)
//    that the original hifiasm applies. Candidates that fail to align are
//    dropped; surviving overlaps are written per aligned window with CIGAR to
//    <prefix>.ovlp.paf. This matches the original hifiasm overlap output
//    exactly (byte-identical after sorting).
//
//  * --dbg-ovec: emit the raw pre-alignment candidate set instead, written to
//    <prefix>.candidates.paf. No base-level alignment is performed, so
//    coordinates are approximate and there is no CIGAR. This is the superset
//    of overlaps before the alignment filter removes non-aligning candidates.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "CommandLines.h"
#include "Process_Read.h"
#include "Hash_Table.h"
#include "Overlaps.h"
#include "htab.h"
#include "kthread.h"
#include "hetmer.h"
#include "ksort.h"
#include "hifiasm_overlaps_internal.h"

// ---------------------------------------------------------------------------
// Globals that the candidate path expects to be defined in the driver TU.
// (In the full assembler these live in Assembly.cpp, which is not part of the
// stripped-down build.)
// ---------------------------------------------------------------------------
All_reads R_INF;

// ---------------------------------------------------------------------------
// Small helpers pulled out of the (removed) Overlaps.cpp. They are the only
// symbols the candidate-detection objects reference from that translation
// unit. Kept here so the 39k-line Overlaps.cpp can be dropped entirely.
// ---------------------------------------------------------------------------
#define generic_key(x) (x)
#define ma_hit_key_tn(a) ((a).tn)
KRADIX_SORT_INIT(ovc_hit_tn, ma_hit_t, ma_hit_key_tn, member_size(ma_hit_t, tn))
KRADIX_SORT_INIT(ovc_arch64, uint64_t, generic_key, 8)

void ma_hit_sort_tn(ma_hit_t *a, long long n)
{
    radix_sort_ovc_hit_tn(a, a + n);
}

void sort_kvec_t_u64_warp(kvec_t_u64_warp* u_vecs, uint32_t is_descend)
{
    radix_sort_ovc_arch64(u_vecs->a.a, u_vecs->a.a + u_vecs->a.n);
    if (is_descend) {
        uint64_t i, t;
        for (i = 0; i < (u_vecs->a.n >> 1); ++i) {
            t = u_vecs->a.a[i];
            u_vecs->a.a[i] = u_vecs->a.a[u_vecs->a.n - i - 1];
            u_vecs->a.a[u_vecs->a.n - i - 1] = t;
        }
    }
}

void init_ma_hit_t_alloc(ma_hit_t_alloc* x)
{
    x->size = 0;
    x->buffer = NULL;
    x->length = 0;
    x->is_fully_corrected = 0;
    x->is_abnormal = 0;
}

// write_dbug/test_dbug are only reached from the unused UL-index persistence
// path in htab.cpp (guarded by ug != NULL, always NULL here). Provide inert
// stubs so the object links.
void write_dbug(ma_ug_t* /*ug*/, FILE* /*fp*/) {}
uint32_t test_dbug(ma_ug_t* /*ug*/, FILE* /*fp*/) { return 0; }

// ---------------------------------------------------------------------------
// Candidate-detection primitive (defined in anchor.cpp): minimizer sketch +
// anchoring against the position index + chaining. Fills overlap_list with
// candidate overlaps. No alignment / correction is performed.
// ---------------------------------------------------------------------------
#define HA_KMER_GOOD_RATIO 0.333
#define COV_W 3072

void h_ec_lchain(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k,
                 All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres,
                 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp,
                 uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t mcopy_num,
                 double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w);

// ---------------------------------------------------------------------------
// Per-thread buffers.
// ---------------------------------------------------------------------------
typedef struct {
    UC_Read self_read;
    Candidates_list clist;
    overlap_region_alloc olist;
    ha_abuf_t *ab;
    st_mt_t sp;
} cand_buf_t;

typedef struct {
    cand_buf_t *a;
    uint32_t n;
} cand_bufs_t;

static cand_bufs_t *cand_bufs_init(uint32_t n)
{
    cand_bufs_t *p = (cand_bufs_t*)calloc(1, sizeof(cand_bufs_t));
    p->n = n;
    p->a = (cand_buf_t*)calloc(n, sizeof(cand_buf_t));
    for (uint32_t k = 0; k < n; ++k) {
        cand_buf_t *z = &p->a[k];
        init_UC_Read(&z->self_read);
        init_Candidates_list(&z->clist);
        init_overlap_region_alloc(&z->olist);
        kv_init(z->sp);
        z->ab = ha_abuf_init();
    }
    return p;
}

static void cand_bufs_destroy(cand_bufs_t *p)
{
    for (uint32_t k = 0; k < p->n; ++k) {
        cand_buf_t *z = &p->a[k];
        destory_UC_Read(&z->self_read);
        destory_Candidates_list(&z->clist);
        destory_overlap_region_alloc(&z->olist);
        kv_destroy(z->sp);
        ha_abuf_destroy(z->ab);
    }
    free(p->a);
    free(p);
}

// ---------------------------------------------------------------------------
// Output: emit candidate overlaps as PAF. Chaining gives approximate query /
// target spans; no CIGAR (no base-level alignment at this stage).
// ---------------------------------------------------------------------------
static pthread_mutex_t g_out_mtx = PTHREAD_MUTEX_INITIALIZER;
static FILE *g_out_fp = NULL;

static void emit_candidates(overlap_region_alloc *ol)
{
    if (ol->length == 0) return;

    // In-memory path: push each overlap into the sink instead of writing PAF.
    // The sink serializes internally, so no g_out_mtx is needed here. Fields
    // match the PAF columns below exactly (same +1 half-open end convention).
    if (hifiasm_ovlp_sink_active()) {
        for (uint64_t k = 0; k < ol->length; ++k) {
            overlap_region *o = &ol->list[k];
            uint32_t span = (o->x_pos_e >= o->x_pos_s) ? (o->x_pos_e - o->x_pos_s + 1) : 0;
            // The raw candidate set carries no base-level CIGAR (that is only
            // produced by the later gen_hc_r_alin_ea alignment pass), so push a
            // zero-length CIGAR: (NULL, 0) tokens, anchor 0.
            hifiasm_ovlp_sink_push(
                o->x_id, o->y_id,
                o->x_pos_s, o->x_pos_e + 1,
                o->y_pos_s, o->y_pos_e + 1,
                o->shared_seed > 0 ? (uint32_t)o->shared_seed : 0,
                span,
                (uint8_t)(o->y_pos_strand == 0),
                NULL, 0, 0);
        }
        return;
    }

    pthread_mutex_lock(&g_out_mtx);
    for (uint64_t k = 0; k < ol->length; ++k) {
        overlap_region *o = &ol->list[k];
        int ql = (int)Get_READ_LENGTH(R_INF, o->x_id);
        int tl = (int)Get_READ_LENGTH(R_INF, o->y_id);
        // PAF: qname qlen qs qe strand tname tlen ts te matches alnlen mapq
        uint32_t span = (o->x_pos_e >= o->x_pos_s) ? (o->x_pos_e - o->x_pos_s + 1) : 0;
        fprintf(g_out_fp,
                "%.*s\t%d\t%u\t%u\t%c\t%.*s\t%d\t%u\t%u\t%u\t%u\t255\ttp:A:c\n",
                (int)Get_NAME_LENGTH(R_INF, o->x_id), Get_NAME(R_INF, o->x_id), ql,
                o->x_pos_s, o->x_pos_e + 1,
                "+-"[o->y_pos_strand],
                (int)Get_NAME_LENGTH(R_INF, o->y_id), Get_NAME(R_INF, o->y_id), tl,
                o->y_pos_s, o->y_pos_e + 1,
                o->shared_seed > 0 ? (uint32_t)o->shared_seed : 0,
                span);
    }
    pthread_mutex_unlock(&g_out_mtx);
}

// ---------------------------------------------------------------------------
// Per-read worker: candidate overlap detection only.
// ---------------------------------------------------------------------------
static void worker_candidates(void *data, long i, int tid)
{
    cand_buf_t *b = &(((cand_bufs_t*)data)->a[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ  = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;

    recover_UC_Read(&b->self_read, &R_INF, i);

    clear_Candidates_list(&b->clist);
    clear_overlap_region_alloc(&b->olist);

    h_ec_lchain(b->ab, i, b->self_read.seq, b->self_read.length, asm_opt.mz_win, asm_opt.k_mer_length,
                &R_INF, &b->olist, &b->clist, ((asm_opt.is_ont) ? (0.05) : (0.02)), asm_opt.max_n_chain,
                1, NULL, NULL, &b->sp, &high_occ, &low_occ, 1, 1, 3, 0.7, 2, 32, COV_W);

    emit_candidates(&b->olist);
}

// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
extern void *ha_flt_tab;
extern ha_pt_t *ha_idx;

// Post-chaining alignment + filter path (defined in ecovlp.cpp). Runs the same
// candidate detector as above and then base-level alignment; candidates that
// fail to align are dropped. Writes <prefix>.ovlp.paf. Only the alignment
// subset of ecovlp.cpp/Correct.cpp is retained in the binary (see Makefile
// --gc-sections); the error-correction/consensus/phasing code is discarded.
void cal_ec_r_dbg(uint64_t n_thre, uint64_t n_a);

// Raw candidate detection (no alignment): writes <prefix>.candidates.paf.
static int detect_candidates_raw(void)
{
    // In-memory path: no PAF file; emit_candidates() pushes into the sink.
    const int toMem = hifiasm_ovlp_sink_active();
    char *paf = NULL;
    if (!toMem) {
        // output file: <output_file_name>.candidates.paf
        paf = (char*)malloc(strlen(asm_opt.output_file_name) + 32);
        sprintf(paf, "%s.candidates.paf", asm_opt.output_file_name);
        g_out_fp = fopen(paf, "w");
        if (!g_out_fp) {
            fprintf(stderr, "[M::%s] ERROR: cannot open %s for writing\n", __func__, paf);
            free(paf);
            return 1;
        }
    }

    cand_bufs_t *bufs = cand_bufs_init(asm_opt.thread_num);
    kt_for(asm_opt.thread_num, worker_candidates, bufs, R_INF.total_reads);
    cand_bufs_destroy(bufs);

    if (!toMem) {
        fclose(g_out_fp);
        g_out_fp = NULL;
        fprintf(stderr, "[M::%s] candidate overlaps written to %s\n", __func__, paf);
        free(paf);
    } else {
        fprintf(stderr, "[M::%s] candidate overlaps collected in memory\n", __func__);
    }
    return 0;
}

// Shared body of the candidate pipeline. `from_store` selects how R_INF is
// populated:
//   from_store == 0 : the read files are streamed (ha_ft_gen writes read
//                     *lengths* into R_INF, then ha_pt_gen reads the sequences);
//                     R_INF is torn down at the end. This is the original CLI /
//                     file-based path.
//   from_store != 0 : R_INF was ALREADY loaded from memory by the bridge
//                     (hifiasm_reads_store_load), so ha_ft_gen/ha_pt_gen must
//                     read sequences from the store (read_from_store=1) rather
//                     than from files, and R_INF is left intact for the caller
//                     to release with hifiasm_reads_store_release().
static int ha_detect_candidates_impl(int from_store)
{
    int hom_cov = -1, het_cov = -1;
    int ret = 0;

    ha_flt_tab = NULL;
    ha_idx = NULL;

    // 1) high-occurrence k-mer filter table.
    //   file path   : ha_ft_gen(read_from_store=0) scans the files and records
    //                 read lengths into R_INF (HAF_RS_WRITE_LEN).
    //   store path  : ha_ft_gen(read_from_store=1) reads sequences back from the
    //                 pre-loaded R_INF (HAF_RS_READ); it writes nothing into it.
    if (!(asm_opt.flag & HA_F_NO_KMER_FLT)) {
        ha_flt_tab = ha_ft_gen(&asm_opt, &R_INF, &hom_cov, 0, from_store);
        ha_opt_update_cov(&asm_opt, hom_cov);
    }

    // 2) position index over the raw reads.
    //   file path   : read_from_store=0. When the filter step ran it wrote only
    //                 read *lengths*, so ha_pt_gen loads the sequences; when the
    //                 filter was skipped R_INF.total_reads==0 and ha_pt_gen does
    //                 the length+sequence load itself.
    //   store path  : read_from_store=1. Sequences are already in R_INF, so
    //                 ha_pt_gen reads them from the store in both of its passes.
    ha_idx = ha_pt_gen(&asm_opt, ha_flt_tab, from_store, 0, &R_INF,
                       &hom_cov, &het_cov);
    asm_opt.hom_cov = hom_cov;
    asm_opt.het_cov = het_cov;
    if (ha_flt_tab == NULL) ha_opt_update_cov(&asm_opt, hom_cov);

    fprintf(stderr, "[M::%s] indexed %lu reads; hom_cov=%d het_cov=%d\n",
            __func__, (unsigned long)R_INF.total_reads, hom_cov, het_cov);
    if (g_hetmer_set != NULL)
        fprintf(stderr, "[M::%s] het-mer seeds force-emitted during indexing: %llu\n",
                __func__, (unsigned long long)g_hetmer_emitted);

    // 3) detect overlaps.
    // Default: candidate detection + base-level alignment/filter, matching the
    // original hifiasm overlap output exactly (writes <prefix>.ovlp.paf).
    // With --dbg-ovec: emit the raw pre-alignment candidate set instead
    // (writes <prefix>.candidates.paf), i.e. before the alignment filter drops
    // non-aligning candidates.
    if (asm_opt.dbg_ovec_cal) {
        ret = detect_candidates_raw();
    } else {
        cal_ec_r_dbg(asm_opt.thread_num, R_INF.total_reads);
        fprintf(stderr, "[M::%s] overlaps written to %s.ovlp.paf\n",
                __func__, asm_opt.output_file_name);
    }

    // Snapshot the read-name table for the in-memory path before R_INF is torn
    // down, so the caller can resolve overlap read indices to its own ids by
    // name. No-op when the sink is inactive (file path).
    hifiasm_ovlp_sink_capture_names(R_INF.name, R_INF.name_index,
                                    R_INF.total_reads, R_INF.total_name_length);

    // 4) cleanup. Always drop the indices/filter. R_INF is torn down only on
    // the file path; on the store path the caller owns R_INF and releases it
    // with hifiasm_reads_store_release() so it can be reused across calls (the
    // same loaded store serves both the filter build and overlap detection).
    ha_pt_destroy(ha_idx); ha_idx = NULL;
    ha_ft_destroy(ha_flt_tab); ha_flt_tab = NULL;
    if (!from_store) destory_All_reads(&R_INF);

    return ret;
}

int ha_detect_candidates(void)
{
    return ha_detect_candidates_impl(/*from_store*/ 0);
}

// Store-fed variant: R_INF must already be loaded (hifiasm_reads_store_load).
// Reads nothing from disk and leaves R_INF intact for the caller to release.
int ha_detect_candidates_from_store(void)
{
    return ha_detect_candidates_impl(/*from_store*/ 1);
}
