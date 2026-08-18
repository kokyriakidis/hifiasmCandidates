#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Correct.h"
#include "Process_Read.h"
#include "ecovlp.h"
#include "kthread.h"
#include "htab.h"
#include "hifiasm_overlaps_internal.h"
#define HA_KMER_GOOD_RATIO 0.333
#define E_KHIT 31
#define CNS_DEL_E (0x7fffffffu)
#define del_cns_arc(z, arc_i) ((z).arc.a[(arc_i)].v == CNS_DEL_E)
#define CNS_DEL_V (0x1fffffffu)
#define del_cns_nn(z, nn_i) ((z).a[(nn_i)].sc == CNS_DEL_V)
#define REFRESH_N 128 
#define COV_W 3072
#define RES_K 19
#define RES_W 19

KDQ_INIT(uint32_t)

typedef struct {
	uint32_t v:31, f:1;
	uint32_t sc;
} cns_arc;
typedef struct {size_t n, m, nou; cns_arc *a; } cns_arc_v;

typedef struct {
	// uint16_t c:2, t:2, f:1, sc:3;
	uint32_t c:2, f:1, sc:29;
	cns_arc_v arc;
}cns_t;

typedef struct {
	size_t n, m; 
	cns_t *a;
	uint32_t si, ei, off, bn, bb0, bb1, cns_g_wl;
	kdq_t(uint32_t) *q;
}cns_gfa;

typedef struct {
	// chaining and overlapping related buffers
	UC_Read self_read, ovlp_read;
	Candidates_list clist;
	overlap_region_alloc olist;
	ha_abuf_t *ab;
	// int64_t num_read_base, num_correct_base, num_recorrect_base;
	uint64_t cnt[6], rr;
	haplotype_evdience_alloc hap;
	bit_extz_t exz;
	kv_ul_ov_t pidx;
	asg64_v v64;
	asg32_v v32;
	asg16_v v16;
	asg8_v v8q, v8t;

    kvec_t_u8_warp k_flag;
	st_mt_t sp;
	cns_gfa cns;	
} ec_ovec_buf_t0;

typedef struct {
	ec_ovec_buf_t0 *a;
	uint32_t n, rev;
    uint8_t *cr;
} ec_ovec_buf_t;

typedef struct {
	ec_ovec_buf_t *p;
    asg64_v idx;
    ma_ug_t *ug;
} ec_polish_buf_t;

typedef struct {
	uint32_t n_thread, n_a, chunk_size, cn;
    FILE *fp;
} cal_ec_r_dbg_t;

// One emitted overlap (whole overlap_region, not per-window) for the in-memory
// sink: full box in the query-forward / target-alignment frame plus a slice into
// this thread's dense-chain arena (och). Collected in the worker where the
// overlap_region (and its z->chain) is in scope, then pushed once per overlap in
// the dump step. The per-window ma_hit_t path below is retained only for the
// file-PAF output.
typedef struct {
    uint32_t q_id, t_id;
    uint32_t q_start, q_end;   // query forward coords (x_pos_s .. x_pos_e+1)
    uint32_t t_start, t_end;   // target forward coords (already rev-adjusted)
    uint32_t n_match, block_len;
    uint8_t  rev;              // 1 => reverse (is_same_strand = !rev)
    uint32_t cig_t_start;      // target anchor in the alignment frame
    uint64_t chain_off;        // offset into och of this overlap's chain
    uint32_t chain_len;        // number of anchors
} r_dbg_ovlp_rec_t;

typedef struct { size_t n, m; r_dbg_ovlp_rec_t *a; } r_dbg_ovlp_vec_t;

typedef struct {
	ma_hit_t *a;
    size_t n, m;
    asg16_v ec;
    // Per-overlap records + their dense-chain anchor arena (only filled when the
    // in-memory sink is active; the file-PAF path leaves these empty).
    r_dbg_ovlp_vec_t ov;       // .a/.n/.m records
    asg64_v och;               // packed (q_start<<32)|t_start anchors
} r_dbg_step_res_t;

typedef struct { // data structure for each step in kt_pipeline()
    ec_ovec_buf_t *buf;
    r_dbg_step_res_t *res;
    uint32_t si, ei;
} cal_ec_r_dbg_step_t;

ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n);
void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p);


#define generic_key(x) (x)
KRADIX_SORT_INIT(ec16, uint16_t, generic_key, 2)
KRADIX_SORT_INIT(ec32, uint32_t, generic_key, 4)
KRADIX_SORT_INIT(ec64, uint64_t, generic_key, 8)

#define kdq_clear(q) ((q)->count = (q)->front = 0)

typedef struct {size_t n, m; asg16_v *a; uint8_t *f; } cc_v;
cc_v scc = {0, 0, NULL, NULL};
cc_v scb = {0, 0, NULL, NULL};
cc_v sca = {0, 0, NULL, NULL};

typedef struct {size_t n, m; char *a; UC_Read z; asg8_v q;} sl_v;


void h_ec_lchain(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t mcopy_num, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w);
void h_ec_lchain_amz(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int max_n_chain, int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t chain_cutoff, uint32_t mcopy_khit_cut, uint64_t ocv_w);
void h_ec_lchain_re_gen(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, ha_pt_t *ha_idx, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								 int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, UC_Read *tu, asg64_v *oidx, asg16_v *scc);
void h_ec_lchain_re_gen3(ha_abuf_t *ab, uint32_t rid, char* rs, uint64_t rl, uint64_t mz_w, uint64_t mz_k, ha_pt_t *ha_idx, All_reads *rref, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								 int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, UC_Read *tu, asg64_v *oidx, asg16_v *scc);
uint64_t get_mz1(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, ha_abuf_t *ab, const void *hf, ha_pt_t *ha_idx, int sample_dist, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, ha_pt_t *pt, int min_freq, int32_t dp_min_len, float dp_e, st_mt_t *mt, int32_t ws, int32_t is_unique, void *km, uint64_t beg_i);
void get_pi_ec_chain(ha_abuf_t *ab, uint64_t rid, uint64_t rl, uint32_t tid, char* ts, uint64_t tl, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, 
								int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, /**uint32_t is_accurate,**/ uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, 
								int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip);
void set_lchain_dp_op(uint32_t is_accurate, uint32_t mz_k, int64_t *max_skip, int64_t *max_iter, int64_t *max_dis, double *chn_pen_gap, double *chn_pen_skip, int64_t *quick_check);
void h_ec_lchain_re_gen_srt(ha_abuf_t *ab, ha_pt_t *ha_idx, overlap_region_alloc *olst, Candidates_list *cl);
uint64_t h_ec_lchain_re_gen_qry(ha_abuf_t *ab, uint64_t *k, uint64_t *l, uint64_t *i, uint64_t *idx_a, uint64_t idx_n, uint64_t *tid, uint64_t *trev);
uint64_t h_ec_lchain_re_chn(ha_abuf_t *ab, uint64_t si, uint64_t ei, uint32_t rid, char* rs, uint64_t rl, uint64_t tid, char* ts, uint64_t tl, uint64_t trev, uint64_t mz_w, uint64_t mz_k, overlap_region_alloc *olst, Candidates_list *cl, double bw_thres, 
								 int apend_be, uint64_t max_cnt, uint64_t min_cnt, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, int64_t max_skip, int64_t max_iter, int64_t max_dis, int64_t quick_check, double chn_pen_gap, double chn_pen_skip, tiny_queue_t *tq, asg16_v *scc, int64_t *n, int64_t *zn);
overlap_region* h_ec_lchain_fast(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, uint64_t mz_w, uint64_t mz_k, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, double bw_thres, 
								 int apend_be, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* dbg_ct, st_mt_t *sp, uint32_t *high_occ, uint32_t *low_occ, uint32_t is_accurate, uint32_t gen_off, int64_t enable_mcopy, double mcopy_rate, uint32_t mcopy_khit_cut, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, double sh);
void h_ec_lchain_fast_new(ha_abuf_t *ab, uint32_t rid, UC_Read *qu, UC_Read *tu, All_reads *rref, overlap_region_alloc *ol, Candidates_list *cl, bit_extz_t *exz, asg16_v *buf, asg64_v *srt_i, ma_hit_t_alloc *in0, ma_hit_t_alloc *in1, double sh);

ec_ovec_buf_t* gen_ec_ovec_buf_t(uint32_t n)
{
    uint32_t k; ec_ovec_buf_t0 *z = NULL;
    ec_ovec_buf_t *p = NULL; CALLOC(p, 1);
    p->n = n; CALLOC(p->a, p->n);
    for (k = 0; k < p->n; k++) {
        z = &(p->a[k]);
        init_UC_Read(&z->self_read);
	    init_UC_Read(&z->ovlp_read);
	    init_Candidates_list(&z->clist);
	    init_overlap_region_alloc(&z->olist);

        // init_fake_cigar(&(z->tmp.f_cigar));
        // memset(&(z->tmp.w_list), 0, sizeof(z->tmp.w_list));
        // CALLOC(z->tmp.w_list.a, 1); z->tmp.w_list.n = z->tmp.w_list.m = 1;

        // kv_init(z->b_buf.a);
        // kv_init(z->r_buf.a);
        kv_init(z->k_flag.a);
        kv_init(z->sp);
        kv_init(z->pidx);
	    kv_init(z->v64);
        kv_init(z->v32);
        kv_init(z->v16);
        kv_init(z->v8q);
        kv_init(z->v8t);
        init_bit_extz_t(&(z->exz), 31);

        z->ab = ha_abuf_init();
        
        InitHaplotypeEvdience(&z->hap);
        z->cns.q = kdq_init(uint32_t);
    }
    
    return p;
}

void destroy_cns_gfa(cns_gfa *p)
{
    size_t k;
    for (k = 0; k < p->m; k++) {
        kv_destroy(p->a[k].arc);
    }
    free(p->a); kdq_destroy(uint32_t, p->q);
}

void destroy_ec_ovec_buf_t(ec_ovec_buf_t *p)
{
    uint32_t k; ec_ovec_buf_t0 *z = NULL;
    for (k = 0; k < p->n; k++) {
        z = &(p->a[k]); z->rr = 0;
        destory_UC_Read(&z->self_read);
        destory_UC_Read(&z->ovlp_read);
        destory_Candidates_list(&z->clist);
	    destory_overlap_region_alloc(&z->olist);

        // destory_fake_cigar(&(z->tmp.f_cigar));
        // free(z->tmp.w_list.a); free(z->tmp.w_list.c.a);

        // kv_destroy(z->r_buf.a);
        kv_destroy(z->k_flag.a);
        kv_destroy(z->sp);
        kv_destroy(z->pidx);
	    kv_destroy(z->v64);
        kv_destroy(z->v32);
        kv_destroy(z->v16);
        kv_destroy(z->v8q);
        kv_destroy(z->v8t);
        destroy_bit_extz_t(&(z->exz));

        ha_abuf_destroy(z->ab);
        
        destoryHaplotypeEvdience(&z->hap);
        destroy_cns_gfa(&(z->cns));

    }
    free(p->a); free(p->cr); free(p);

    // fprintf(stderr, "[M::%s-chains] #->%lld\n", __func__, asm_opt.num_bases);
    // fprintf(stderr, "[M::%s-passed-chains-0] #->%lld\n", __func__, asm_opt.num_corrected_bases);
    // fprintf(stderr, "[M::%s-cis-chains-1] #->%lld\n", __func__, asm_opt.num_recorrected_bases);
}

inline void refresh_ec_ovec_buf_t0(ec_ovec_buf_t0 *z, uint64_t n)
{
    z->rr++;
    if((z->rr%n) == 0) {
        free(z->self_read.seq); memset(&(z->self_read), 0, sizeof(z->self_read));
        free(z->ovlp_read.seq); memset(&(z->ovlp_read), 0, sizeof(z->ovlp_read));

        destory_Candidates_list(&z->clist); memset(&(z->clist), 0, sizeof(z->clist));
        destory_overlap_region_alloc(&z->olist); memset(&(z->olist), 0, sizeof(z->olist)); init_overlap_region_alloc(&z->olist);

        kv_destroy(z->k_flag.a); kv_init(z->k_flag.a);
        kv_destroy(z->sp); kv_init(z->sp);
        kv_destroy(z->pidx); kv_init(z->pidx);
        kv_destroy(z->v64); kv_init(z->v64);
        kv_destroy(z->v32); kv_init(z->v32);
        kv_destroy(z->v16); kv_init(z->v16);
        kv_destroy(z->v8q); kv_init(z->v8q);
        kv_destroy(z->v8t); kv_init(z->v8t);

        destroy_bit_extz_t(&(z->exz)); init_bit_extz_t(&(z->exz), 31);

        ha_abuf_destroy(z->ab); z->ab = ha_abuf_init();

        destoryHaplotypeEvdience(&z->hap); memset(&(z->hap), 0, sizeof(z->hap)); InitHaplotypeEvdience(&z->hap);

        destroy_cns_gfa(&(z->cns)); memset(&(z->cns), 0, sizeof(z->cns)); z->cns.q = kdq_init(uint32_t);

        // z->rr = 1;
    }
}



overlap_region *fetch_aux_ovlp(overlap_region_alloc* ol) /// exactly same to gen_aux_ovlp
{
	if (ol->length + 1 >= ol->size) {
		uint64_t sl = ol->size;
        ol->size = ol->length + 1;
        kroundup64(ol->size);
        REALLOC(ol->list, ol->size);
        /// need to set new space to be 0
        memset(ol->list + sl, 0, sizeof(overlap_region)*(ol->size - sl));
	}
    ///debug for memory
    // if(ol->length + 1 >= ol->size) {
    //     fprintf(stderr, "[M::%s] length::%lu, size::%lu\n", __func__, ol->length, ol->size);
    // }
	return &(ol->list[ol->length+1]);
}

typedef struct {
	ul_ov_t *c_idx;
    asg64_v *idx;
    int64_t i, i0, srt_n, rr, ru;
    uint64_t mms, mme;
} cc_idx_t;


///[s, e)

#define simp_vote_len 6

///[s, e)

typedef struct {
	All_reads *rref;
    UC_Read *tu;    
    uint64_t s, e, n0, n1, id, rev;
} rr_seq_t;

inline void insert_cns_arc(cns_gfa *cns, uint32_t src, uint32_t des, uint32_t is_ou, uint32_t plus0, uint32_t rid)
{
    if(src >= cns->n) {
        fprintf(stderr, "[M::%s] rid::%u, src::%u, des::%u, (*cns).n::%u\n", __func__, rid, src, des, (uint32_t)(*cns).n);
        exit(1);
    }
    cns_arc *p, t; kv_pushp(cns_arc, (*cns).a[src].arc, &p);
    p->f = 0; p->sc = plus0; p->v = des;
    if(is_ou) {
        (*cns).a[src].arc.nou++;
        if((*cns).a[src].arc.nou < (*cns).a[src].arc.n) {
            t = (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1]; 
            (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1] = *p;
            *p = t;
        }
    } 
}

inline uint32_t insert_cns_node(cns_gfa *cns)
{
    cns_t *p; uint32_t m0; 
    if (((*cns)).n == ((*cns)).m) { 
        m0 = ((*cns)).m;
        ((*cns)).m = ((*cns)).m? ((*cns)).m<<1 : 2; 
        ((*cns)).a = (cns_t*)realloc(((*cns)).a, sizeof(cns_t) * ((*cns)).m); 
        if(((*cns)).m > m0) {
            memset(((*cns)).a + m0, 0, sizeof(cns_t)*(((*cns)).m-m0));
        }
    } 
    *(&p) = &((*cns)).a[((*cns)).n++];
    p->arc.n = p->arc.nou = 0;
    p->c = p->f = p->sc = 0;
    return ((*cns)).n - 1;
}

inline uint32_t add_cns_arc(cns_gfa *cns, uint32_t src, uint32_t des, uint32_t is_ou, uint32_t plus)
{
    uint32_t k, s, e;
    if(is_ou) {
        s = 0; e = (*cns).a[src].arc.nou;
    } else {
        s = (*cns).a[src].arc.nou; e = (*cns).a[src].arc.n;
    }

    for (k = s; k < e; k++) {
        if((*cns).a[src].arc.a[k].v == des) {
            (*cns).a[src].arc.a[k].sc += plus;
            break;
        }
    }
    
    return ((k < e)?(1):(0));
}

inline void prt_cns_arc(cns_gfa *cns, uint32_t src, const char* cmd)
{
    uint32_t k;
    fprintf(stderr, "\n%s\t[M::%s] src::%u, sc::%u, c::%u\n", cmd, __func__, src, (*cns).a[src].sc, (*cns).a[src].c);
    for (k = 0; k < (*cns).a[src].arc.n; k++) {
        fprintf(stderr, "%s\t[M::%s] des::%u, sc::%u, is_ou::%u\n", cmd, __func__, (*cns).a[src].arc.a[k].v, (*cns).a[src].arc.a[k].sc, k<(*cns).a[src].arc.nou?1:0);
    }
}

inline uint32_t get_cns_arc_bp(cns_gfa *cns, uint32_t src, uint32_t bp, uint32_t is_ou, uint32_t av_bp)
{
    uint32_t k, s, e;
    if(is_ou) {
        s = 0; e = (*cns).a[src].arc.nou;
    } else {
        s = (*cns).a[src].arc.nou; e = (*cns).a[src].arc.n;
    }

    for (k = s; k < e; k++) {
        if((*cns).a[src].arc.a[k].v == 0 || (*cns).a[src].arc.a[k].v == 1) continue;
        if(av_bp && (*cns).a[src].arc.a[k].v >= (*cns).bb0 && (*cns).a[src].arc.a[k].v < (*cns).bb1) continue;///no backbone
        if((*cns).a[(*cns).a[src].arc.a[k].v].c == bp) {
            return k;
        }
    }
    
    return ((uint32_t)-1);
}

inline uint32_t add_cns_arc_bp(cns_gfa *cns, uint32_t src, uint32_t bp, uint32_t plus0, uint32_t rid, uint32_t av_bp)
{
    uint32_t rr, des;
    rr = get_cns_arc_bp(cns, src, bp, 1, av_bp);
    if(rr != ((uint32_t)-1)) {///find an existing node
        des = (*cns).a[src].arc.a[rr].v;
        (*cns).a[des].sc++;
        (*cns).a[src].arc.a[rr].sc += plus0;

        rr = add_cns_arc(cns, des, src, 0, plus0); 
        // if(rr == 0) {
        //     fprintf(stderr, "[M::%s] src::%u -> des::%u\n", __func__, src, des);
        //     prt_cns_arc(cns, src);
        //     prt_cns_arc(cns, des);
        // }
        assert(rr);

        return des;
    } else {///create a new node
        des = insert_cns_node(cns);
        (*cns).a[des].sc++; (*cns).a[des].c = bp;
        insert_cns_arc(cns, src, des, 1, plus0, rid); 
        insert_cns_arc(cns, des, src, 0, plus0, rid);
    }

    return des;
}


///[s, e)


///[s, e)




///[s, e)







inline void gen_mm_cns_arc(cns_gfa *cns, uint32_t src, uint32_t des, uint32_t sc, uint32_t f)
{
    cns_t *av = &((*cns).a[src]), *aw;
    uint32_t vk, wk;
    for (vk = 0; vk < av->arc.nou; vk++) {///out-edge of src
        if((av->arc.a[vk].v != des) || (del_cns_arc((*av), vk))) continue;
        // av->arc.a[vk].f = 1; ///not sure if we should set these edges as visited
        av->arc.a[vk].f = f;
        av->arc.a[vk].sc += sc;

        aw = &((*cns).a[des]);
        for (wk = aw->arc.nou; wk < aw->arc.n; wk++) {///in-edge of des
            if((aw->arc.a[wk].v != src) || (del_cns_arc((*aw), wk))) continue;
            // aw->arc.a[wk].f = 1; ///not sure if we should set these edges as visited
            aw->arc.a[wk].f = f;
            aw->arc.a[wk].sc += sc; 
            break;
        }

        assert(wk < aw->arc.n);
        return;
    }

    cns_arc *p, t;
    ///src -> des
    kv_pushp(cns_arc, (*cns).a[src].arc, &p);
    p->sc = sc; p->v = des;
    // p->f = 1; ///not sure if we should set these edges as visited
    p->f = f; 
    ///ou-edge
    (*cns).a[src].arc.nou++;
    if((*cns).a[src].arc.nou < (*cns).a[src].arc.n) {
        t = (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1]; 
        (*cns).a[src].arc.a[(*cns).a[src].arc.nou-1] = *p;
        *p = t;
    }

    ///src <- des; in-edge
    kv_pushp(cns_arc, (*cns).a[des].arc, &p);
    p->sc = sc; p->v = src;
    // p->f = 1; ///not sure if we should set these edges as visited
    p->f = f; 
}









///no e_end since e always covers end; s may not have end






















inline uint64_t exact_ec_check(char *qstr, uint64_t ql, char *tstr, uint64_t tl, int64_t qs, int64_t qe, int64_t ts, int64_t te)
{
    if(qe - qs != te - ts) return 0;
    if(memcmp(qstr + qs, tstr + ts, qe - qs) == 0) return 1;
    return 0;
}

void gen_hc_r_alin_ea(overlap_region_alloc* ol, Candidates_list *cl, All_reads *rref, UC_Read* qu, UC_Read* tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t wl, int64_t rid, int64_t khit, int64_t move_gap, asg16_v *buf, asg64_v *srt, ma_hit_t_alloc *in, uint8_t chem_drop, double align_gap_rate, int64_t align_gap_max)
{
    if(ol->length <= 0) return;


    // uint64_t k, l, i, s, m, mm_k, *ei, en, *oi, on, tid, trev, nec; int64_t sc, mm_sc, plus, minus; overlap_region *z, t; ma_hit_t *p;
    uint64_t k, i, m, *ei, en, *oi, on, tid, trev, nec; overlap_region *z; ma_hit_t *p;
    srt->n = 0;
    for (k = 0; k < in->length; k++) {
        if(in->buffer[k].el) {
            m = in->buffer[k].tn; m <<= 1; m |= in->buffer[k].rev; 
            m <<= 32; m |= k; kv_push(uint64_t, (*srt), m);
        }
    }

    if(!(srt->n)) {
        gen_hc_r_alin(ol, cl, rref, qu, tu, exz, aux_o, e_rate, wl, rid, khit, move_gap, buf, chem_drop, align_gap_rate, align_gap_max);
    } else {
        ///debug for memory
        // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

        kv_resize(uint64_t, *srt, (srt->n + ol->length));
        ei = srt->a; en = srt->n; oi = srt->a + srt->n; on = ol->length;
        for (k = 0; k < on; k++) {
            z = &(ol->list[k]); z->is_match = z->strong = z->without_large_indel = 0;
            oi[k] = z->y_id; oi[k] <<= 1; oi[k] |= z->y_pos_strand;
            oi[k] <<= 32; oi[k] |= k;
        }

        radix_sort_ec64(ei, ei + en); radix_sort_ec64(oi, oi + on);
        for (k = i = nec = 0; k < on; k++) {
            z = &(ol->list[(uint32_t)oi[k]]); tid = z->y_id; trev = z->y_pos_strand;
            for (; (i < en) && ((ei[i]>>32) < ((tid<<1)|trev)); i++);
            if((i < en) && ((ei[i]>>32) == ((tid<<1)|trev))) {
                p = &(in->buffer[(uint32_t)ei[i]]);
                if((z->x_pos_s == ((uint32_t)p->qns)) && (z->x_pos_e + 1 == p->qe) && 
                            (z->y_pos_s == p->ts) && (z->y_pos_e + 1 == p->te)) {
                    resize_UC_Read(tu, p->te - p->ts); recover_UC_Read_sub_region(tu->seq, p->ts, p->te - p->ts, trev, rref, tid);
                    if(exact_ec_check(qu->seq, qu->length, tu->seq, p->te - p->ts, ((uint32_t)p->qns), p->qe, 0, p->te - p->ts)) {
                        z->is_match = 1; z->shared_seed = z->non_homopolymer_errors;///for index
                        z->non_homopolymer_errors = 0; z->strong = z->without_large_indel = 0;
                        set_exact_exz(exz, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1); push_alnw(z, exz);
                        nec++;
                    }
                }
            }
        }
        ///debug for memory
        // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

        if(on > nec) {
            gen_hc_r_alin_nec(ol, cl, rref, qu, tu, exz, aux_o, e_rate, wl, rid, khit, move_gap, buf, chem_drop, align_gap_rate, align_gap_max);
        }

        // fprintf(stderr, "[M::%s] srt->n::%u, nec::%lu, on::%lu\n", __func__, (uint32_t)srt->n, nec, on);
        ///debug for memory
        // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);
    }

    /**
    if(ol->length > 1) {///for duplicated chains
        overlap_region_sort_y_id(ol->list, ol->length);
        for (k = 1, l = m = 0; k <= ol->length; k++) {
            if(k == ol->length || ol->list[k].y_id != ol->list[l].y_id) {
                // fprintf(stderr, "\n[M::%s::tid->%u] n->%lu\n", __func__, ol->list[l].y_id, k - l);
                mm_k = l;
                if(k - l > 1) {
                    for (s = l, mm_sc = INT32_MIN, mm_k = ((uint64_t)-1); s < k; s++) {
                        z = &(ol->list[s]);
                        plus = z->x_pos_e + 1 - z->x_pos_s; minus = (z->non_homopolymer_errors) * 12;
                        sc = plus - minus;
                        if((sc > mm_sc) || ((sc == mm_sc) && ((ol->list[mm_k].x_pos_e+1-ol->list[mm_k].x_pos_s) < (z->x_pos_e+1-z->x_pos_s)))) {
                            mm_sc = sc; mm_k = s;
                        }
                        // fprintf(stderr, "[M::%s::%c] q::[%u, %u), t::[%u, %u), sc::%ld, err::%u, s::%lu\n", __func__, "+-"[z->y_pos_strand], z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1, sc, z->non_homopolymer_errors, s);
                    }
                }
                // fprintf(stderr, "[M::%s::tid->%u] mm_k::%lu\n", __func__, ol->list[l].y_id, mm_k);
                if(mm_k != ((uint64_t)-1)) {
                    if(mm_k != m) {
                        t = ol->list[mm_k];
                        ol->list[mm_k] = ol->list[m];
                        ol->list[m] = t;
                    }
                    m++;
                }
                l = k;
            }
        }
        ol->length = m;
    }
    **/
}











static void worker_hap_ec_dbg_paf(void *data, long i, int tid)
{
	ec_ovec_buf_t0 *b = &(((cal_ec_r_dbg_step_t*)data)->buf->a[tid]);
    r_dbg_step_res_t *rr = &(((cal_ec_r_dbg_step_t*)data)->res[tid]);
    uint32_t high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
    uint32_t low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
    overlap_region *aux_o = NULL; i += ((cal_ec_r_dbg_step_t*)data)->si;

    // debug_retrive_bqual(D, &b->v8t, i, 256); return;

    recover_UC_Read(&b->self_read, &R_INF, i); 

    h_ec_lchain(b->ab, i, b->self_read.seq, b->self_read.length, asm_opt.mz_win, asm_opt.k_mer_length, &R_INF, &b->olist, &b->clist, ((asm_opt.is_ont)?(0.05):(0.02)), asm_opt.max_n_chain, 1, NULL, NULL, &(b->sp), &high_occ, &low_occ, 1, 1, 3, 0.7, 2, 32, COV_W);///ONT high error

    aux_o = fetch_aux_ovlp(&b->olist);///must be here

    // stderr_phase_ovlp(&b->olist);
    gen_hc_r_alin_ea(&b->olist, &b->clist, &R_INF, &b->self_read, &b->ovlp_read, &b->exz, aux_o, asm_opt.max_ov_diff_ec, (asm_opt.is_ont)?(WINDOW_OHC):(WINDOW_HC), i, E_KHIT, 1, &b->v16, &b->v64, &(R_INF.paf[i]), 0, -1, -1);

    // When the in-memory sink is active, dinara consumes one record per overlap
    // (full box + native dense chain), not per window. Collect those here where
    // z->chain is in scope; the per-window ma_hit_t loop below still runs for the
    // file-PAF path.
    const int toMem = hifiasm_ovlp_sink_active();

    uint32_t k, m, tl; overlap_region *z; bit_extz_t ez; ma_hit_t *t; 
    for (k = 0; k < b->olist.length; k++) {
        z = &(b->olist.list[k]);
        if(!(z->w_list.n)) continue;

        tl = Get_READ_LENGTH((R_INF), z->y_id);

        if (toMem) {
            // One per-overlap record for the sink. Box: forward query
            // (x_pos_s..x_pos_e+1) vs forward-strand target (rev-adjusted like
            // the window path). The native chain (z->chain) is in the alignment
            // frame (forward query, alignment-orientation target); dinara
            // reframes it exactly like the CIGAR using cig_t_start = y_pos_s
            // (the overlap's alignment-orientation target start).
            r_dbg_ovlp_rec_t *ov;
            kv_pushp(r_dbg_ovlp_rec_t, rr->ov, &ov);
            ov->q_id = z->x_id; ov->t_id = z->y_id;
            ov->q_start = z->x_pos_s; ov->q_end = z->x_pos_e + 1;
            ov->rev = (uint8_t)z->y_pos_strand;
            if (!ov->rev) {
                ov->t_start = z->y_pos_s;
                ov->t_end   = z->y_pos_e + 1;
            } else {
                ov->t_start = tl - (z->y_pos_e + 1);
                ov->t_end   = tl - z->y_pos_s;
            }
            ov->cig_t_start = z->y_pos_s; // alignment-orientation target start
            ov->n_match = z->shared_seed > 0 ? (uint32_t)z->shared_seed : 0;
            ov->block_len = ov->q_end - ov->q_start;
            // Copy this overlap's dense chain into the thread-local arena.
            ov->chain_off = rr->och.n;
            ov->chain_len = z->chain.length;
            if (z->chain.length) {
                kv_resize(uint64_t, rr->och, rr->och.n + z->chain.length);
                memcpy(rr->och.a + rr->och.n, z->chain.buffer,
                       z->chain.length * sizeof(uint64_t));
                rr->och.n += z->chain.length;
            }
        }

        for (m = 0; m < z->w_list.n; m++) {
            if(is_ualn_win(z->w_list.a[m])) continue;
            set_bit_extz_t(ez, (*z), m);
            kv_pushp(ma_hit_t, *rr, &t);

            t->qns = z->x_id; t->qns = t->qns << 32;
            t->tn = z->y_id;

            t->qns = t->qns | (uint64_t)(z->w_list.a[m].x_start);
            t->qe = z->w_list.a[m].x_end + 1;

            t->ts = z->w_list.a[m].y_start;
            t->te = z->w_list.a[m].y_end + 1;

            t->rev = z->y_pos_strand;
            
            t->bl = rr->ec.n;
            kv_resize(uint16_t, rr->ec, rr->ec.n + ez.cigar.n);
            memcpy(rr->ec.a + rr->ec.n, ez.cigar.a, ez.cigar.n * sizeof((*(rr->ec.a))));
            rr->ec.n += ez.cigar.n;
            t->cc = rr->ec.n - t->bl;

            if(t->rev) {
                t->ts = tl - z->w_list.a[m].y_end - 1;
                t->te = tl - z->w_list.a[m].y_start;
            }
        }
    }
}
























































void print_ov_dbg_paf(FILE *fp, char *ref_str, char *ref_id, int32_t ref_id_n, char *qry_str, char *qry_id, int32_t qry_id_n, uint64_t rs, uint64_t re, uint64_t rl, uint64_t qs, uint64_t qe, uint64_t ql, uint64_t rev, bit_extz_t *ez, char *ezh)
{
    // Emit a standard-compliant PAF record. ezh maps cigar ops to characters
    // as {'=' match, 'X' mismatch, 'I' insertion, 'D' deletion}. The number of
    // residue matches (col 10) is the total '=' length; the alignment block
    // length (col 11) is the total over all ops. The extended CIGAR is carried
    // in the optional cg:Z: tag (col 13), so downstream tools that read the
    // standard 12 columns work unchanged.
    uint64_t ci = 0; uint16_t c; uint32_t cl;
    uint64_t n_match = 0, blk_len = 0;
    for (ci = 0; ci < ez->cigar.n; ) {
        ci = pop_trace(&(ez->cigar), ci, &c, &cl);
        if (ezh[c] == '=') n_match += cl;
        blk_len += cl;
    }

    fprintf(fp, "%.*s\t%lu\t%lu\t%lu\t", qry_id_n, qry_id, ql, qs, qe);
    fprintf(fp, "%c\t", "+-"[rev]);
    fprintf(fp, "%.*s\t%lu\t%lu\t%lu\t", ref_id_n, ref_id, rl, rs, re);
    // col 10 matches, col 11 alignment block length, col 12 mapping quality.
    fprintf(fp, "%lu\t%lu\t255\tcg:Z:", n_match, blk_len);
    for (ci = 0; ci < ez->cigar.n; ) {
        ci = pop_trace(&(ez->cigar), ci, &c, &cl);
        fprintf(fp, "%u%c", cl, ezh[c]);
    }
    fprintf(fp, "\n");
}

static void *worker_ov_dbg_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    cal_ec_r_dbg_t *p = (cal_ec_r_dbg_t*)data; char cm[4]; cm[0] = '='; cm[1] = 'X'; cm[2] = 'I'; cm[3] = 'D'; ///op0=match('='), op1=mismatch('X'), op2=ins('I'), op3=del('D')
    // cal_ec_r_dbg_step_t
    if (step == 0) { // step 1: read a block of sequences
        cal_ec_r_dbg_step_t *s; CALLOC(s, 1);
        s->si = p->cn; p->cn += p->chunk_size; 
        if(p->cn > p->n_a) p->cn = p->n_a; s->ei = p->cn;
        if(s->si >= s->ei) free(s);
        else return s;
    } else if (step == 1) { // step 2: alignment
        cal_ec_r_dbg_step_t *s = (cal_ec_r_dbg_step_t*)in;
        s->buf = gen_ec_ovec_buf_t(p->n_thread); CALLOC(s->res, p->n_thread);  
        kt_for(p->n_thread, worker_hap_ec_dbg_paf, s, (s->ei - s->si));
        destroy_ec_ovec_buf_t(s->buf);
        return s;
    } else if (step == 2) { // step 3: dump
        cal_ec_r_dbg_step_t *s = (cal_ec_r_dbg_step_t*)in; uint64_t k, z; bit_extz_t ez; memset(&ez, 0, sizeof(ez));
        // UC_Read qu; UC_Read tu; init_UC_Read(&qu); init_UC_Read(&tu); 
        const int toMem = hifiasm_ovlp_sink_active();
        if (toMem) {
            // Per-overlap emission: one sink record per overlap_region, carrying
            // the full box and hifiasm's native dense chain. The per-window
            // ma_hit_t records (s->res[k].a) are not used on this path; they are
            // built only for the file-PAF branch below. No base CIGAR is pushed
            // (dinara re-derives alignment from the chain), so cigar=(NULL,0).
            for (k = 0; k < p->n_thread; k++) {
                r_dbg_step_res_t *r = &s->res[k];
                for (z = 0; z < r->ov.n; z++) {
                    r_dbg_ovlp_rec_t *o = &r->ov.a[z];
                    const uint64_t *chain =
                        o->chain_len ? (r->och.a + o->chain_off) : NULL;
                    hifiasm_ovlp_sink_push(
                        o->q_id, o->t_id,
                        o->q_start, o->q_end,
                        o->t_start, o->t_end,
                        o->n_match, o->block_len,
                        (uint8_t)(o->rev == 0),
                        /*cigar*/ NULL, /*cigar_len*/ 0, o->cig_t_start,
                        chain, o->chain_len);
                }
                free(r->a); free(r->ec.a);
                free(r->ov.a); free(r->och.a);
            }
            free(s->res); free(s);
            return 0;
        }
        for (k = 0; k < p->n_thread; k++) {
            for (z = 0; z < s->res[k].n; z++) {
                ez.cigar.a = s->res[k].ec.a + s->res[k].a[z].bl; ez.cigar.n = ez.cigar.m = s->res[k].a[z].cc;

                // UC_Read_resize(qu, (s->res[k].a[z].qe - ((uint32_t)s->res[k].a[z].qns))); 
                // recover_UC_Read_sub_region(qu.seq, ((uint32_t)s->res[k].a[z].qns), (s->res[k].a[z].qe - ((uint32_t)s->res[k].a[z].qns)), 0, &R_INF, (s->res[k].a[z].qns>>32));

                // UC_Read_resize(tu, (s->res[k].a[z].te - s->res[k].a[z].ts));
                // recover_UC_Read_sub_region(tu.seq, s->res[k].a[z].ts, s->res[k].a[z].te - s->res[k].a[z].ts, 0, &R_INF, s->res[k].a[z].tn); 

                // File-PAF path only (the in-memory sink path returned early
                // above and emits per-overlap, not per-window).
                print_ov_dbg_paf(p->fp, NULL/**qu.seq**/, Get_NAME(R_INF, (s->res[k].a[z].qns>>32)), Get_NAME_LENGTH(R_INF, (s->res[k].a[z].qns>>32)), 
                            NULL/**tu.seq**/, Get_NAME(R_INF, (s->res[k].a[z].tn)), Get_NAME_LENGTH(R_INF, (s->res[k].a[z].tn)), (uint32_t)s->res[k].a[z].qns, s->res[k].a[z].qe, Get_READ_LENGTH(R_INF, (s->res[k].a[z].qns>>32)), s->res[k].a[z].ts, s->res[k].a[z].te, Get_READ_LENGTH(R_INF, (s->res[k].a[z].tn)), s->res[k].a[z].rev, &ez, cm);
            }
            free(s->res[k].a); free(s->res[k].ec.a);
        }
        free(s->res); free(s); ///destory_UC_Read(&qu); destory_UC_Read(&tu);
    }
    return 0;
}

void cal_ec_r_dbg(uint64_t n_thre, uint64_t n_a)
{
    // In-memory path: the dump step pushes into the sink, so no PAF file is
    // opened and sl.fp stays NULL.
    const int toMem = hifiasm_ovlp_sink_active();
    char *paf = NULL;

    cal_ec_r_dbg_t sl; memset(&sl, 0, sizeof(sl));
    sl.n_thread = n_thre; sl.n_a = n_a; sl.chunk_size = 2000; sl.cn = 0;
    if (!toMem) {
        MALLOC(paf, (strlen(asm_opt.output_file_name)+64));
        sprintf(paf, "%s.ovlp.paf", asm_opt.output_file_name);
        sl.fp = fopen(paf, "w");
    }

    kt_pipeline(3, worker_ov_dbg_pipeline, &sl, 3);

    if (!toMem) { fclose(sl.fp); free(paf); }
}








