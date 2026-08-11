#define __STDC_LIMIT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include "Correct.h"
#include "Levenshtein_distance.h"
#include "Assembly.h"
#include "CommandLines.h"
// #include "ksw2.h"
#include "ksort.h"
#include "kalloc.h"
#include "htab.h"
#include "Overlaps.h"
#include "inter.h"
#include "Process_Read.h"
#define A_L 16
#define ext_w 6
#define r_simi_w 0.05
#define rphase_thres 4

#define generic_key(x) (x)
KRADIX_SORT_INIT(b32, uint32_t, generic_key, 4)
KRADIX_SORT_INIT(bc64, uint64_t, generic_key, 8)

#define haplotype_evdience_key(x) ((x).site)
KRADIX_SORT_INIT(haplotype_evdience_srt, haplotype_evdience, haplotype_evdience_key, member_size(haplotype_evdience, site))

#define haplotype_evdience_id_key(x) ((x).overlapID)
KRADIX_SORT_INIT(haplotype_evdience_id_srt, haplotype_evdience, haplotype_evdience_id_key, member_size(haplotype_evdience, overlapID))

#define haplotype_evdience_os_key(x) ((x).overlapSite)
KRADIX_SORT_INIT(haplotype_evdience_os_srt, haplotype_evdience, haplotype_evdience_os_key, member_size(haplotype_evdience, overlapSite))

#define overlap_region_dp_key(x) ((x).x_pos_e)
KRADIX_SORT_INIT(overlap_region_dp_srt, overlap_region, overlap_region_dp_key, member_size(overlap_region, x_pos_e))

#define window_list_xs_key(x) ((x).x_start)
KRADIX_SORT_INIT(window_list_xs_srt, window_list, window_list_xs_key, member_size(window_list, x_start))

#define uov_qs_key(p) ((p).qs)
KRADIX_SORT_INIT(uov_srt_qs, ul_ov_t, uov_qs_key, member_size(ul_ov_t, qs))

#define k_mer_hit_self_key(p) ((p).self_offset)
KRADIX_SORT_INIT(k_mer_hit_self, k_mer_hit, k_mer_hit_self_key, member_size(k_mer_hit, self_offset))

#define k_mer_hit_off_key(p) ((p).offset)
KRADIX_SORT_INIT(k_mer_hit_off, k_mer_hit, k_mer_hit_off_key, member_size(k_mer_hit, offset))

#define ul_ov_srt_qs1_key(p) ((p).qs)
KRADIX_SORT_INIT(ul_ov_srt_qs1, ul_ov_t, ul_ov_srt_qs1_key, member_size(ul_ov_t, qs))

#define ul_ov_srt_tn1_key(p) ((p).tn)
KRADIX_SORT_INIT(ul_ov_srt_tn1, ul_ov_t, ul_ov_srt_tn1_key, member_size(ul_ov_t, tn))

#define ul_ov_srt_qn1_key(p) ((p).qn)
KRADIX_SORT_INIT(ul_ov_srt_qn1, ul_ov_t, ul_ov_srt_qn1_key, member_size(ul_ov_t, qn))

#define ul_ov_srt_qe1_key(p) ((p).qe)
KRADIX_SORT_INIT(ul_ov_srt_qe1, ul_ov_t, ul_ov_srt_qe1_key, member_size(ul_ov_t, qe))

#define ul_ov_srt_ts1_key(p) ((p).ts)
KRADIX_SORT_INIT(ul_ov_srt_ts1, ul_ov_t, ul_ov_srt_ts1_key, member_size(ul_ov_t, ts))

#define MAX_SEC_ERR (0x3fffffffU)


int ha_ov_type(const overlap_region *r, uint32_t len);
void set_lchain_dp_op(uint32_t is_accurate, uint32_t mz_k, int64_t *max_skip, int64_t *max_iter, int64_t *max_dis, double *chn_pen_gap, double *chn_pen_skip, int64_t *quick_check);







inline int get_interval(long long window_start, long long window_end, overlap_region_alloc* overlap_list, Correct_dumy* dumy, long long blockLen)
{
    uint64_t i, fud = 0;
    long long Len;
    if(window_start == 0) dumy->start_i = 0;
    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///this interval is smaller than all overlaps
        ///in this case, the next interval should start from 0
        if (window_end < (long long)overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            dumy->lengthNT = 0;
            return 0;
        }
        else ///if window_end >= overlap_list->list[i].x_pos_s，this overlap might be overlapped with current interval
        {
            dumy->start_i = i;
            break;
        }
    }


    ///this interval is larger than all overlaps, so we don't need to scan next overlap
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        dumy->lengthNT = 0;
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;
    fud = 0;

    for (; i < overlap_list->length; i++)
    {
        if((Len = OVERLAP(window_start, window_end, (long long)overlap_list->list[i].x_pos_s, (long long)overlap_list->list[i].x_pos_e)) > 0)
        {
            ///sometimes the length of window > WINDOW, but overlap length == WINDOW
            // if (Len == WINDOW && window_end - window_start + 1 == WINDOW)
            if (Len == blockLen && window_end - window_start + 1 == blockLen)
            {
                dumy->overlapID[dumy->length] = i;
                dumy->length++; 
            }
            else
            {
                dumy->lengthNT++;
                dumy->overlapID[dumy->size - dumy->lengthNT] = i;
            }
            if(fud == 0) fud = 1, dumy->start_i = i;
        }
        
        if((long long)overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    if ( dumy->length + dumy->lengthNT == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


inline int get_available_interval(long long window_start, long long window_end, overlap_region_alloc* overlap_list, Correct_dumy* dumy)
{
    uint64_t i, fud = 0;
    long long Len;
    if(window_start == 0) dumy->start_i = 0;
    for (i = dumy->start_i; i < overlap_list->length; i++)
    {
        ///this interval is smaller than all overlaps
        ///in this case, the next interval should start from 0
        if (window_end < (long long)overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            dumy->length = 0;
            dumy->lengthNT = 0;
            return 0;
        }
        else ///if window_end >= overlap_list->list[i].x_pos_s，this overlap might be overlapped with current interval
        {
            dumy->start_i = i;
            break;
        }
    }

    ///this interval is larger than all overlaps, so we don't need to scan next overlap
    if (i >= overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        dumy->length = 0;
        dumy->lengthNT = 0;
        return -2;
    }
    
    dumy->length = 0;
    dumy->lengthNT = 0;
    fud = 0;

    
    long long fake_length = 0;

    for (; i < overlap_list->length; i++)
    {
        ///check if the interval is overlapped with current overlap   
        if((Len = OVERLAP(window_start, window_end, (long long)overlap_list->list[i].x_pos_s, (long long)overlap_list->list[i].x_pos_e)) > 0)
        {
            ///number of overlaps
            fake_length++;

            ///check if this overlap is available
            if (overlap_list->list[i].is_match == 1)
            {
                dumy->overlapID[dumy->length] = i;
                dumy->length++; 
            }
            if(fud == 0) fud = 1, dumy->start_i = i;
        }
        
        if((long long)overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    ///fake_length is the number of overlaps, instead of the number of available overlaps
    if (fake_length == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

///Len = OVERLAP(window_start, window_end, overlap_list->list[i].x_pos_s, overlap_list->list[i].x_pos_e))





void fill_subregion(char* r, long long start_pos, long long length, uint8_t strand, All_reads* R_INF, long long ID, 
int extra_begin, int extra_end)
{
    
    recover_UC_Read_sub_region(r+extra_begin, start_pos, length, strand, R_INF, ID);
    memset(r, 'N', extra_begin);
    memset(r+extra_begin+length, 'N', extra_end);
}

void fill_subregion_ul(char* r, long long start_pos, long long length, uint8_t strand, const ul_idx_t *uref, long long ID, 
int extra_begin, int extra_end)
{
    retrieve_u_seq(NULL, r+extra_begin, &(uref->ug->u.a[ID]), strand, start_pos, length, NULL);
    memset(r, 'N', extra_begin);
    memset(r+extra_begin+length, 'N', extra_end);
}

int determine_overlap_region(int threshold, long long y_start, long long y_ID, long long Window_Len, /**All_reads* R_INF**/long long y_len,
int* r_extra_begin, int* r_extra_end, long long* r_y_start, long long* r_y_length)
{
    int extra_begin;
    int extra_end;
    long long currentIDLen;
    long long o_len;

    ///the length of y
    // currentIDLen = Get_READ_LENGTH((*R_INF), y_ID);
    currentIDLen = y_len;
    
    ///since Window_Len == x_len + (threshold << 1)
    if(y_start < 0 || currentIDLen <= y_start || 
    currentIDLen - y_start + 2 * threshold + THRESHOLD_MAX_SIZE < Window_Len)
    {
        return 0;
    }
    
    extra_begin = extra_end = 0;
    ///y maybe less than 0
    y_start = y_start - threshold;
    o_len = MIN(Window_Len, currentIDLen - y_start);
    extra_end = Window_Len - o_len;

    if (y_start < 0)
    {
        extra_begin = -y_start;
        y_start = 0;
        o_len = o_len - extra_begin;
    }

    (*r_extra_begin) = extra_begin;
    (*r_extra_end) = extra_end;
    (*r_y_start) = y_start;
    (*r_y_length) = o_len;

    return 1;
}






int32_t init_waln(int64_t err, int64_t s, int64_t l, int64_t w_l, int64_t* aux_beg, int64_t* aux_end, int64_t* r_s, int64_t* r_l)
{
    (*aux_beg) = (*aux_end) = (*r_s) = (*r_l) = -1;
    ///since w_l == x_len + (err << 1)
    if((s < 0) || (s >= l) || ((l-s+(2*err)+THRESHOLD_MAX_SIZE) < w_l)) return 0;
    (*aux_beg) = (*aux_end) = 0;
    ///s might be less than 0
    (*r_s) = s - err; 
    (*r_l) = l-(*r_s); if((*r_l) > w_l) (*r_l) = w_l;
    (*aux_end) = w_l - (*r_l);

    if ((*r_s) < 0) {
        (*aux_beg) = -(*r_s); (*r_s) = 0; (*r_l) -= (*aux_beg);
    }
    return 1;
}

///[s, e)
int64_t get_num_wins(int64_t s, int64_t e, int64_t block_s)
{
    int64_t nl = e - ((s/block_s)*block_s), nw;
    nw = (nl/block_s); if((nl%block_s)>0) nw++;
    return nw;
}



///error_rate should be 30%

inline int double_error_threshold(int pre_threshold, int x_len)
{

    pre_threshold = Adjust_Threshold(pre_threshold, x_len);
    int threshold = pre_threshold * 2;
    ///may have some bugs
    if(x_len >= 300 && threshold < THRESHOLD_MAX_SIZE)
    {
        threshold = THRESHOLD_MAX_SIZE;
    }
    
    if(threshold > THRESHOLD_MAX_SIZE)
    {
        threshold = THRESHOLD_MAX_SIZE;
    }

    return threshold;
}

inline int double_ul_error_threshold(int pre_threshold, int x_len)
{

    pre_threshold = Adjust_Threshold(pre_threshold, x_len);
    int threshold = THRESHOLD_UL_MAX * x_len;
    if(threshold < pre_threshold) threshold = pre_threshold;
    if(threshold > THRESHOLD_MAX_SIZE) threshold = THRESHOLD_MAX_SIZE;
    return threshold;
}


inline int verify_sub_window(All_reads* R_INF, Correct_dumy* dumy, UC_Read* g_read, 
long long x_beg, long long xLen, long long y_beg, long long yLen, uint64_t y_id, 
uint64_t y_pos_strand, int threshold, int alignment_strand,
unsigned int* get_error, int* get_y_end, int* get_x_end, int* get_aligned_xLen)
{
    (*get_aligned_xLen) = 0;
    (*get_y_end) = -1;
    (*get_x_end) = -1;
    (*get_error) = (unsigned int)-1;

    int extra_begin, extra_end, r_x_end, r_y_end, aligned_xLen;
    long long o_len;
    unsigned int r_error;
    if(!determine_overlap_region(threshold, y_beg, y_id, yLen, 
    Get_READ_LENGTH((*R_INF), y_id), &extra_begin, &extra_end, &y_beg, &o_len))
    {
        return 0;
    }

    fill_subregion(dumy->overlap_region, y_beg, o_len, y_pos_strand, R_INF, 
    y_id, extra_begin, extra_end);

    char* x_string = g_read->seq + x_beg;
    char* y_string = dumy->overlap_region;

    aligned_xLen = 0;

    alignment_extension(y_string, yLen, x_string, xLen, threshold, 
                    alignment_strand, &r_error, &r_y_end, &r_x_end, &aligned_xLen);
    
    (*get_error) = r_error;
    (*get_y_end) = r_y_end;
    (*get_x_end) = r_x_end;
    (*get_aligned_xLen) = aligned_xLen;

    if(aligned_xLen == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

inline int verify_ul_sub_window(const ul_idx_t *uref, Correct_dumy* dumy, UC_Read* g_read, 
long long x_beg, long long xLen, long long y_beg, long long yLen, uint64_t y_id, 
uint64_t y_pos_strand, int threshold, int alignment_strand,
unsigned int* get_error, int* get_y_end, int* get_x_end, int* get_aligned_xLen)
{
    (*get_aligned_xLen) = 0;
    (*get_y_end) = -1;
    (*get_x_end) = -1;
    (*get_error) = (unsigned int)-1;

    int extra_begin, extra_end, r_x_end, r_y_end, aligned_xLen;
    long long o_len;
    unsigned int r_error;
    if(!determine_overlap_region(threshold, y_beg, y_id, yLen, 
    uref->ug->u.a[y_id].len, &extra_begin, &extra_end, &y_beg, &o_len))
    {
        return 0;
    }

    fill_subregion_ul(dumy->overlap_region, y_beg, o_len, y_pos_strand, uref, 
    y_id, extra_begin, extra_end);

    // if(y_id == 6) {
    //     fprintf(stderr, "-[M::%s::aln_dir->%d] qs->%lld, ts->%lld, thres->%d, aux_beg->%d, aux_end->%d, t_pri_l->%lld\n", 
    //     __func__, alignment_strand, x_beg, y_beg, threshold, extra_begin, extra_end, o_len);
    // }

    char* x_string = g_read->seq + x_beg;
    char* y_string = dumy->overlap_region;

    aligned_xLen = 0;

    alignment_extension(y_string, yLen, x_string, xLen, threshold, 
                    alignment_strand, &r_error, &r_y_end, &r_x_end, &aligned_xLen);
    
    (*get_error) = r_error;
    (*get_y_end) = r_y_end;
    (*get_x_end) = r_x_end;
    (*get_aligned_xLen) = aligned_xLen;

    if(aligned_xLen == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

inline int64_t get_init_err_thres(int64_t len, double e_rate, int64_t block_s, int64_t block_err)
{
    if(len >= block_s) return block_err;
    int64_t thres = len * e_rate;
    thres = Adjust_Threshold(thres, len);
    if(thres > THRESHOLD_MAX_SIZE) thres = THRESHOLD_MAX_SIZE;
    return thres;
}


uint32_t get_init_paras(All_reads* rref, const ul_idx_t *uref, overlap_region *z, int64_t x_s, int64_t x_e, double e_rate, int64_t block_s, 
                        int64_t *r_ys, int64_t *r_ex_beg, int64_t *r_ex_end, int64_t *r_err_thre)
{
    int e, ex_beg, ex_end; long long y_s, o_len, Window_Len;
    e = get_init_err_thres(x_e+1-x_s, e_rate, block_s, rref?THRESHOLD:THRESHOLD_MAX_SIZE);
    y_s = (x_s-z->x_pos_s) + z->y_pos_s; y_s += y_start_offset(x_s, &(z->f_cigar));
    Window_Len = (x_e+1-x_s) + (e<<1);

    if(!determine_overlap_region(e, y_s, z->y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), z->y_id)):(uref->ug->u.a[z->y_id].len)),
    &ex_beg, &ex_end, &y_s, &o_len)) {
        return 0;
    }
    (*r_ys) = y_s; (*r_ex_beg) = ex_beg; (*r_ex_end) = ex_end; (*r_err_thre) = e;
    return 1;
}

int64_t check_coverage_gap(uint64_t *v_idx, uint64_t w_s, uint64_t w_e, int64_t block_s)
{
    int64_t wid = w_s/block_s, a_n = (uint32_t)(v_idx[wid]), k;
    uint64_t *a = v_idx + (v_idx[wid]>>32);
    for (k = 0; k < a_n; k++) {
        if(((a[k]>>32) == w_s) && (((uint32_t)(a[k])) == w_e)) return 1;
    }
    return 0;
}

inline double non_trim_error_rate(overlap_region *z, All_reads* rref, const ul_idx_t *uref, const kvec_t_u64_warp* v_idx, Correct_dumy* dumy, UC_Read* g_read, double e_rate, int64_t block_s)
{
    int64_t nw, aw = z->w_list.n, k, m, w_id, wn_id, w_s, w_e, idx_e, tErr = 0, tLen = 0, y_s, ex_beg, ex_end, err_thre, p_err_thre;
    int64_t x_len, Window_Len, y_beg_left, y_beg_right;
    unsigned int r_error_left, r_error_right; int32_t r_x_end_left, r_y_end_left, aligned_xLen_left, r_x_end_right, r_y_end_right, aligned_xLen_right;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s);
    assert(nw >= aw && aw > 0);

    for (k = aw-1, idx_e = nw; k >= 0; k--) {
        w_id = get_win_id_by_e(z, z->w_list.a[k].x_end, block_s, &w_s);
        assert(w_s == z->w_list.a[k].x_start && w_id < idx_e && k <= w_id);
        tLen += z->w_list.a[k].x_end + 1 - z->w_list.a[k].x_start;
        tErr += z->w_list.a[k].error;///matched window
        // if(z->y_id == 1) {
        //     fprintf(stderr, "+[M::%s] ws->%d, we->%d, tot_l->%ld, tot_e->%ld\n", 
        //                         __func__, z->w_list.a[k].x_start, z->w_list.a[k].x_end, tLen, tErr);
        // }
        // if(k != w_id) z->w_list.a[w_id] = z->w_list.a[k];
        ///from mapped window w_list.a[k] to the following unmapped windows
        for (m = w_id+1, w_e = z->w_list.a[k].x_end; m < idx_e; m++) {
            w_s = w_e + 1;
            wn_id = get_win_id_by_s(z, w_s, block_s, &w_e);
            assert(wn_id == m); x_len = w_e + 1 - w_s; tLen += x_len;
            ///check if there are some windows that cannot be algined by any overlaps/unitigs
            ///if no, it is likely that the UL read itself has issues
            if(uref && v_idx && z->is_match == 4) {
                if(check_coverage_gap(v_idx->a.a, w_s, w_e, block_s)) {
                    tErr += THRESHOLD_MAX_SIZE;
                    // if(z->y_id == 1) {
                    //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                    // }
                    continue;
                } 
            }
            if(!get_init_paras(rref, uref, z, w_s, w_e, e_rate, block_s, &y_s, &ex_beg, &ex_end, &err_thre)) {
                tErr += x_len;
                // if(z->y_id == 1) {
                //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                // }
                continue;
            }
            p_err_thre = err_thre;

            if(rref) {
                err_thre = double_error_threshold(err_thre, x_len);
            } else {
                err_thre = double_ul_error_threshold(err_thre, x_len);
            }
            Window_Len = x_len + (err_thre << 1); 
            r_error_left = r_error_right = 0;
            aligned_xLen_left = aligned_xLen_right = 0;
            y_beg_left = y_beg_right = -1;

            if(m == w_id+1) { ///if the previous window is mapped
                y_beg_left = z->w_list.a[k].y_end + 1;///incorrect
            }

            if(m+1 == idx_e && k+1 < aw) { ///if the next window is mapped
                y_beg_right = z->w_list.a[k+1].y_start-x_len;///incorrect
            }

            if(y_beg_left == -1 && y_beg_right == -1) {
                y_beg_left = y_s;
                if(ex_beg >= 0) y_beg_left = y_beg_left + p_err_thre - ex_beg;
                y_beg_right = y_beg_left;
            }

            if(y_beg_left == -1 && y_beg_right != -1) y_beg_left = y_beg_right;
            if(y_beg_right == -1 && y_beg_left != -1) y_beg_right = y_beg_left;

            if(y_beg_left != -1) {///note: this function will change tstr/qstr
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                }
            }

            if(y_beg_right != -1) {
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                }
            }

            ///aligned in both directions
            if(aligned_xLen_left != 0 && aligned_xLen_right != 0) {
                if(aligned_xLen_left + aligned_xLen_right <= x_len) {
                    tErr += r_error_left + r_error_right + (x_len - aligned_xLen_left - aligned_xLen_right);
                } else {
                    float E_rate = (float)(x_len)/(float)(aligned_xLen_left + aligned_xLen_right);
                    tErr += (r_error_left + r_error_right)*E_rate;
                }
            }///not aligned in both directions
            else if(aligned_xLen_left == 0 && aligned_xLen_right == 0) {
                tErr += x_len;
            }///only aligned in left
            else if(aligned_xLen_left != 0) {
                tErr += r_error_left + (x_len - aligned_xLen_left);
            }///only aligned in right
            else if(aligned_xLen_right != 0) {
                tErr += r_error_right + (x_len - aligned_xLen_right);
            }
            // if(z->y_id == 1) {
            //     fprintf(stderr, "*[M::%s] qs->%ld, ts->%ld, tb[0]->%ld, tb[1]->%ld, di[0]->%u, di[1]->%u, al[0]->%d, al[1]->%d, err_thre->%ld\n", __func__, 
            //     w_s, y_s, y_beg_left, y_beg_right, r_error_left, r_error_right, aligned_xLen_left, aligned_xLen_right, err_thre);
            // }
            // if(z->y_id == 1) {
            //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
            // }
        }
        idx_e = w_id;
    }

    if(idx_e > 0) {
        for (m = 0, w_e = (int64_t)z->x_pos_s-1; m < idx_e; m++) {
            w_s = w_e + 1;
            wn_id = get_win_id_by_s(z, w_s, block_s, &w_e);
            assert(wn_id == m); x_len = w_e + 1 - w_s; tLen += x_len;
            ///check if there are some windows that cannot be algined by any overlaps/unitigs
            ///if no, it is likely that the UL read itself has issues
            if(uref && v_idx && z->is_match == 4) {
                if(check_coverage_gap(v_idx->a.a, w_s, w_e, block_s)) {
                    tErr += THRESHOLD_MAX_SIZE;
                    // if(z->y_id == 1) {
                    //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                    // }
                    continue;
                } 
                // else {
                //     if(z->y_id == 575) {
                //         fprintf(stderr, "---[M::%s::] z::y_id->%u, w_s->%ld, w_e->%ld\n", __func__, z->y_id, w_s, w_e);
                //     } 
                // }
            }
            if(!get_init_paras(rref, uref, z, w_s, w_e, e_rate, block_s, &y_s, &ex_beg, &ex_end, &err_thre)) {
                tErr += x_len;
                // if(z->y_id == 1) {
                //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
                // }
                continue;
            }
            p_err_thre = err_thre;
            if(rref) {
                err_thre = double_error_threshold(err_thre, x_len);
            } else {
                err_thre = double_ul_error_threshold(err_thre, x_len);
            }
            Window_Len = x_len + (err_thre << 1); 
            r_error_left = r_error_right = 0;
            aligned_xLen_left = aligned_xLen_right = 0;
            y_beg_left = y_beg_right = -1;

            ///impossible that the previous window is mapped
            // if(m == w_id+1) { ///if the previous window is mapped
            //     y_beg_left = z->w_list.a[k].y_end + 1;
            // }

            if(m+1 == idx_e && k+1 < aw) { ///if the next window is mapped
                y_beg_right = z->w_list.a[k+1].y_start-x_len;
            }

            if(y_beg_left == -1 && y_beg_right == -1) {
                y_beg_left = y_s;
                if(ex_beg >= 0) y_beg_left = y_beg_left + p_err_thre - ex_beg;
                y_beg_right = y_beg_left;
            }

            if(y_beg_left == -1 && y_beg_right != -1) y_beg_left = y_beg_right;
            if(y_beg_right == -1 && y_beg_left != -1) y_beg_right = y_beg_left;

            if(y_beg_left != -1) {
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_left, Window_Len, 
                    z->y_id, z->y_pos_strand, err_thre, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
                }
            }

            if(y_beg_right != -1) {
                if(rref) {
                    verify_sub_window(rref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                } else {
                    verify_ul_sub_window(uref, dumy, g_read, w_s, x_len, y_beg_right, Window_Len, z->y_id, z->y_pos_strand,
                    err_thre, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
                }
            }

            ///aligned in both directions
            if(aligned_xLen_left != 0 && aligned_xLen_right != 0) {
                if(aligned_xLen_left + aligned_xLen_right <= x_len) {
                    tErr += r_error_left + r_error_right + (x_len - aligned_xLen_left - aligned_xLen_right);
                } else {
                    float E_rate = (float)(x_len)/(float)(aligned_xLen_left + aligned_xLen_right);
                    tErr += (r_error_left + r_error_right)*E_rate;
                }
            }///not aligned in both directions
            else if(aligned_xLen_left == 0 && aligned_xLen_right == 0) {
                tErr += x_len;
            }///only aligned in left
            else if(aligned_xLen_left != 0) {
                tErr += r_error_left + (x_len - aligned_xLen_left);
            }///only aligned in right
            else if(aligned_xLen_right != 0) {
                tErr += r_error_right + (x_len - aligned_xLen_right);
            }
            // if(z->y_id == 1) {
            //     fprintf(stderr, "-[M::%s] ws->%ld, we->%ld, tot_l->%ld, tot_e->%ld\n", __func__, w_s, w_e, tLen, tErr);
            // }
        }
    }
    
    assert(tLen == z->x_pos_e + 1 - z->x_pos_s);
    return (double)(tErr)/(double)(tLen);
}


/**
inline double non_trim_ul_error_rate(overlap_region_alloc* overlap_list, long long ID,
const ul_idx_t *uref, Correct_dumy* dumy, UC_Read* g_read)
{
    long long tLen, tError,i, subWinLen, subWinNum;
    
    tLen = 0;
    tError = 0;

    subWinNum = overlap_list->list[ID].w_list_length;

    
    for (i = 0; i < subWinNum; i++)
    {
        subWinLen = overlap_list->list[ID].w_list[i].x_end - overlap_list->list[ID].w_list[i].x_start + 1;
        tLen += subWinLen;

        if(overlap_list->list[ID].w_list[i].y_end != -1)
        {
            tError += overlap_list->list[ID].w_list[i].error;
        }
        else
        {
            int x_len = subWinLen;
            int threshold = double_ul_error_threshold(overlap_list->list[ID].w_list[i].error_threshold, x_len);
            int Window_Len = x_len + (threshold << 1);
            unsigned int r_error_left = 0;
            int r_x_end_left, r_y_end_left, aligned_xLen_left;
            unsigned int r_error_right = 0;
            int r_x_end_right, r_y_end_right, aligned_xLen_right;
            long long y_beg_left, y_beg_right;
            
            aligned_xLen_left = aligned_xLen_right = 0;
            y_beg_left = y_beg_right = -1;

            if(overlap_list->list[ID].w_list[i].y_start == -1)
            {
                tError += x_len;
                continue;
            }

            ///if the previous window is mapped
            if(i > 0 && overlap_list->list[ID].w_list[i - 1].y_end != -1)
            {
                y_beg_left = overlap_list->list[ID].w_list[i - 1].y_end + 1;
            }

            ///if the next window is mapped
            if(i < (long long)(overlap_list->list[ID].w_list_length - 1) && overlap_list->list[ID].w_list[i + 1].y_end != -1)
            {
                y_beg_right = 1 + overlap_list->list[ID].w_list[i + 1].y_start - 1 - x_len;
            }


            if(y_beg_left == -1 && y_beg_right == -1)
            {
                y_beg_left = overlap_list->list[ID].w_list[i].y_start;
                if(overlap_list->list[ID].w_list[i].extra_begin >= 0)
                {
                    y_beg_left = y_beg_left + overlap_list->list[ID].w_list[i].error_threshold - 
                    overlap_list->list[ID].w_list[i].extra_begin;
                }
                y_beg_right = y_beg_left;
            }

            if(y_beg_left == -1 && y_beg_right != -1)
            {
                y_beg_left = y_beg_right;
            }

            if(y_beg_right == -1 && y_beg_left != -1)
            {
                y_beg_right = y_beg_left;
            }

            
            if(y_beg_left != -1)
            {
                verify_ul_sub_window(uref, dumy, g_read, overlap_list->list[ID].w_list[i].x_start, 
                x_len, y_beg_left, Window_Len, overlap_list->list[ID].y_id, overlap_list->list[ID].y_pos_strand, 
                threshold, 0, &r_error_left, &r_y_end_left, &r_x_end_left, &aligned_xLen_left);
            }

            if(y_beg_right != -1)
            {
                verify_ul_sub_window(uref, dumy, g_read, overlap_list->list[ID].w_list[i].x_start, 
                x_len, y_beg_right, Window_Len, overlap_list->list[ID].y_id, overlap_list->list[ID].y_pos_strand,
                threshold, 1, &r_error_right, &r_y_end_right, &r_x_end_right, &aligned_xLen_right);
            }

            ///aligned in both direction
            if(aligned_xLen_left != 0 && aligned_xLen_right != 0)
            {
                if(aligned_xLen_left + aligned_xLen_right <= x_len)
                {
                    tError = tError + r_error_left + r_error_right + 
                    (x_len - aligned_xLen_left - aligned_xLen_right);
                }
                else
                {
                    float E_rate = (float)(x_len)/(float)(aligned_xLen_left + aligned_xLen_right);
                    tError = tError + (r_error_left + r_error_right)*E_rate;
                }
            }///not aligned in both direction
            else if(aligned_xLen_left == 0 && aligned_xLen_right == 0)
            {
                tError += x_len;
            }///only aligned in left
            else if(aligned_xLen_left != 0)
            {
                tError = tError + r_error_left + (x_len - aligned_xLen_left);
            }///only aligned in right
            else if(aligned_xLen_right != 0)
            {
                tError = tError + r_error_right + (x_len - aligned_xLen_right);
            }
        }
    }

    double error_rate = (double)(tError)/(double)(tLen);

    return error_rate;
}

int calculate_hpm_errors(char* x, int x_len, char* y, int y_len, CIGAR* cigar, int error)
{
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;
    int hpm_error = 0;


    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        if (operation == 0)
        {
            x_i = x_i + operationLen;
            y_i = y_i + operationLen;
        }
        else if (operation == 1)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {
                if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
                {
                    hpm_error++;
                }

                x_i++;
                y_i++;
            }
        }
        else if (operation == 2)
        {

            if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
            {
                hpm_error++;
            }
            cigar_error += operationLen;
            y_i += operationLen;
        }
        else if (operation == 3)
        {

            if(if_is_homopolymer_repeat(x_i, x, x_len) || if_is_homopolymer_repeat(y_i, y, y_len))
            {
                hpm_error++;
            }

            cigar_error += operationLen;
            x_i += operationLen;
        }
        
        cigar_i++;
    }
    return hpm_error;
}

int verify_cigar(char* x, int x_len, char* y, int y_len, CIGAR* cigar, int error)
{
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    int cigar_error = 0;
    int flag_error = 0;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 means there are more y, 3 means there are more x
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]!=y[y_i])
                {
                    ///fprintf(stderr, "error match\n");
                    flag_error = 1;
                }
                x_i++;
                y_i++;
            }
        }
        else if (operation == 1)
        {
            cigar_error += operationLen;
            for (i = 0; i < operationLen; i++)
            {

                if (x[x_i]==y[y_i])
                {
                    ///fprintf(stderr, "error mismatch, cigar_i: %d, x_i: %d, y_i: %d\n",cigar_i, x_i, y_i);
                    flag_error = 1;
                }
                x_i++;
                y_i++;
            }
        }
        else if (operation == 2)
        {
            cigar_error += operationLen;
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            cigar_error += operationLen;
            x_i += operationLen;
        }
        
        cigar_i++;
    }

    
    if (cigar_error != error)
    {
        
        // fprintf(stderr, "error cigar_error: cigar_error: %d, error: %d\n", cigar_error, error);
        // for (i = 0; i < cigar->length; i++)
        // {
        //     fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        // }
        
        
       flag_error = 1;
        
    }
    
   
    if (flag_error == 1)
    {
        
        // print_string(x, x_len);
        // print_string(y, y_len);
        // fprintf(stderr, "x_len: %d, y_len: %d, cigar_len: %d, error: %d\n", x_len, y_len, cigar->length, error);
        // for (i = 0; i < cigar->length; i++)
        // {
        //     fprintf(stderr, "%u: %u\n", cigar->C_L[i], cigar->C_C[i]);
        // }
        
        
    }
    
    
    return flag_error;
    
}
**/

int32_t scan_cigar(window_list *idx, window_list_alloc *cc, int64_t* get_error, int64_t scanXLen, int64_t direction)
{
    uint8_t c = (uint8_t)-1; uint32_t cl = (uint32_t)-1;
    (*get_error) = -1;
    if(idx->clen == 1) {
        get_cigar_cell(idx, cc, 0, &c, &cl);
        if(c == 0) {
            (*get_error) = 0;
            return 1;
        }
    }
    int32_t x_i = 0, y_i = 0, c_i, c_n = idx->clen, c_err = 0;
    uint32_t i;

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2: there are more bases at y, 3: there are more bases at x
    if(direction == 0) {
        for (c_i = 0; c_i < c_n; c_i++) {
            get_cigar_cell(idx, cc, c_i, &c, &cl);
            if (c == 0) { //match
                x_i += cl; y_i += cl;
                if(x_i >= scanXLen) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
            else if (c == 1) {
                for (i = 0; i < cl; i++) {
                    x_i++; y_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            }
            else if (c == 2) {///y has more bases than x
                c_err += cl; y_i += cl;
            }
            else if (c == 3) {///x has more bases than y
                for (i = 0; i < cl; i++) {
                    x_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            }
        }
    } else {
        for (c_i = c_n-1; c_i >= 0; c_i--) {
            get_cigar_cell(idx, cc, c_i, &c, &cl);
            if (c == 0) { //match
                x_i += cl; y_i += cl;
                if(x_i >= scanXLen) {
                    (*get_error) = c_err;
                    return 1;
                }
            } else if (c == 1) { //mismatch
                for (i = 0; i < cl; i++) {
                    x_i++; y_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            } else if (c == 2) {///y has more bases than x
                c_err += cl; y_i += cl;
            } else if (c == 3) {///x has more bases than y
                for (i = 0; i < cl; i++) {
                    x_i++; c_err++;
                    if(x_i >= scanXLen) {
                        (*get_error) = c_err;
                        return 1;
                    }
                }
            }
        }
    }

    (*get_error) = c_err;
    return 0;
}

///[scanXbeg, scanXend]
int scan_cigar_interval(window_list *idx, window_list_alloc *cc, int64_t* get_error, int64_t scanXbeg, int64_t scanXend)
{
    uint8_t c; uint32_t cl;
    (*get_error) = -1;
     if(idx->clen == 1) {
        get_cigar_cell(idx, cc, 0, &c, &cl);
        if(c == 0) {
            (*get_error) = 0;
            return 1;
        }
    }



    int32_t x_i = 0, y_i = 0, c_i, c_n = idx->clen, c_err = 0;
    uint32_t i;
    

    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2: there are more bases at y, 3: there are more bases at x
    for (c_i = 0; c_i < c_n; c_i++) {
        get_cigar_cell(idx, cc, c_i, &c, &cl);
        if (c == 0) {//match
            for (i = 0; i < cl; i++) {
                if(x_i == scanXbeg) c_err = 0;
                x_i++; y_i++;

                if(x_i == scanXend + 1) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
        } else if (c == 1) {//mismatch
            for (i = 0; i < cl; i++) {
                if(x_i == scanXbeg) c_err = 0;
                x_i++; y_i++; c_err++;

                if(x_i == scanXend + 1) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
        } else if (c == 2) {///y has more bases than x
            c_err += cl; y_i += cl;
        }
        else if (c == 3) {
            for (i = 0; i < cl; i++) {
                if(x_i == scanXbeg) c_err = 0;
                x_i++; c_err++;

                if(x_i == scanXend + 1) {
                    (*get_error) = c_err;
                    return 1;
                }
            }
        }
    }

    (*get_error) = c_err;
    return 0;
}

inline int move_gap_greedy(char* path, int path_i, int path_length, char* x, int x_i, char* y, int y_i, unsigned int* new_error)
{
    if(path[path_i] < 2)
    {
        return 0;
    }

    /**
         * 
    GGCG-TGTGCCTGT
        *
    GGCAATGTGCCTGT
        *
    00013000000000
    **/

    int flag = 0;

    char oper = path[path_i];
    

    if(oper == 3)///there are more x
    {
        path_i++;
        y_i--;
        for (; path_i < path_length && x_i >= 0 && y_i >= 0; path_i++, x_i--, y_i--)
        {
            if(path[path_i] == 2 || path[path_i] == 3 || (path[path_i] == 0 && x[x_i] != y[y_i]))
            {
                break;
            }
            else ///path[path_i] = 1 || path[path_i] = 0, exchange path[path_i] with path[path_i-1]
            {
                if(path[path_i] == 1 && x[x_i] == y[y_i])
                {
                    path[path_i - 1] = 0;
                    (*new_error)--;
                }
                else
                {
                    path[path_i - 1] = path[path_i];
                }

                path[path_i] = oper;
                

                flag = 1;
            }
            
        }
    }
    else if(oper == 2)///there are more y
    {
        path_i++;
        x_i--;
        for (; path_i < path_length && x_i >= 0 && y_i >= 0; path_i++, x_i--, y_i--)
        {
            if(path[path_i] == 2 || path[path_i] == 3 || (path[path_i] == 0 && x[x_i] != y[y_i]))
            {
                break;
            }
            else
            {

                if(path[path_i] == 1 && x[x_i] == y[y_i])
                {
                    path[path_i - 1] = 0;
                    (*new_error)--;
                }
                else
                {
                    path[path_i - 1] = path[path_i];
                }


                path[path_i] = oper;
                flag = 1;
            }
            
        }
    }

    return flag;
}

inline void generate_cigar(char* path, int path_length, window_list *idx, window_list_alloc *res, int* start, int* end, unsigned int* old_error,
    char* x, int x_len, char* y)
{
    // uint8_t debug_c; uint32_t debug_c_len;
    idx->cidx = res->c.n;
    if ((*old_error) == 0) {
        push_cigar_cell(res, 0, idx->x_end + 1 - idx->x_start);
        idx->clen = res->c.n - idx->cidx;
        // get_cigar_cell(idx, res, idx->clen-1, &debug_c, &debug_c_len);
        // assert(debug_c==0 && debug_c_len==(idx->x_end + 1 - idx->x_start));
        return;
    }

    ///0 is match, 1 is mismatch, 2 is up (more y), 3 is left (more x)
	int32_t i = 0, pre_cl = 0, trem_p = -1; char pre_c = 5;    
    for (i = 0; i < path_length; i++) {
        if(path[i] == 1) {
            path[i] = 3;(*end)--; trem_p = i;
        } else {
            break;
        }
    }

    for (i = path_length - 1; i >= 0; i--) {
        if(path[i] == 1) {
            path[i] = 3; (*start)++;
        }
        else {
            break;
        }
    }


    // for (i = path_length - 1; i >= 0; i--)
    // {

    //     if (pre_ciga != path[i])
    //     {
    //         if (pre_ciga_length != 0)
    //         {
    //             result->cigar.C_L[result->cigar.length] = pre_ciga_length;
    //             result->cigar.C_C[result->cigar.length] = pre_ciga;
    //             result->cigar.length++;
    //         }

    //         pre_ciga = path[i];
    //         pre_ciga_length = 1;
    //     }
    //     else
    //     {
    //         pre_ciga_length++;
    //     }
    // }

    // if (pre_ciga_length != 0)
    // {
    //     result->cigar.C_L[result->cigar.length] = pre_ciga_length;
    //     result->cigar.C_C[result->cigar.length] = pre_ciga;
    //     result->cigar.length++;
    // }

    ///verify_cigar(x, x_len, y + (*start), (*end) - (*start) + 1, &(result->cigar), error);
    
    
    y = y + (*start);
    int32_t x_i = 0, y_i = 0;
    ///terminate_site = -1 in default
    for (i = path_length - 1; i > trem_p; i--) {
        if(path[i] == 0) {
            x_i++; y_i++;
        }
        else if(path[i] == 1) {
            x_i++; y_i++;
        }
        else if(path[i] == 2) {///there are more y
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, old_error);
            y_i++;
        }
        else if(path[i] == 3) {///there are more x
            move_gap_greedy(path, i, path_length, x, x_i, y, y_i, old_error);
            x_i++;
        }
    }



    pre_c = 5; pre_cl = 0;
    for (i = path_length - 1; i >= 0; i--) {
        if (pre_c != path[i]) {
            if (pre_cl != 0) {
                push_cigar_cell(res, pre_c, pre_cl);
                // get_cigar_cell(idx, res, res->c.n - idx->cidx - 1, &debug_c, &debug_c_len);
                // assert(debug_c==pre_c && debug_c_len==pre_cl);
            }
            pre_c = path[i]; pre_cl = 1;
        }
        else {
            pre_cl++;
        }
    }

    if (pre_cl != 0) {
        push_cigar_cell(res, pre_c, pre_cl);   
        // get_cigar_cell(idx, res, res->c.n - idx->cidx -1, &debug_c, &debug_c_len);
        // assert(debug_c==pre_c && debug_c_len==pre_cl);
    }

    idx->clen = res->c.n - idx->cidx;
    // if(verify_cigar(x, x_len, y, (*end) - (*start) + 1, &(result->cigar), *old_error))
    // {
    //     fprintf(stderr, "error\n");
    // }
}


inline int fix_ul_boundary(char* x_string, long long x_len, int threshold,
long long total_y_start, long long local_y_start, long long local_y_end,
long long old_extra_begin, long long old_extra_end,
long long y_ID, long long Window_Len, const ul_idx_t *uref,
Correct_dumy* dumy, int y_strand, unsigned int old_error,
long long* r_total_y_start, int* r_start_site, int* r_end_site,
int* r_extra_begin, int* r_extra_end, unsigned int* r_error)
{

    
    int new_extra_begin, new_extra_end;
    long long new_y_start, new_y_length;
    int new_end_site, new_start_site;
    unsigned int new_error;
    char* y_string;


    int path_length;

    ///if the start pos at the left boundary 
    if(local_y_start == 0)
    {
        total_y_start = total_y_start + local_y_start;
        ///if local_y_start == 0 and old_extra_begin != 0
        ///this means total_y_start == 0, so shift to the left cannot get a new start pos
        if(old_extra_begin != 0)
        {
            return 0;
        }

        ///if the begining of alignment is 0, we should try to shift the window to find a better result
        ///shift to the left by threshold-1 bases
        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, uref->ug->u.a[y_ID].len,
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        ///if new_y_start is equal to total_y_start, recalculate makes no sense
        if(new_y_start == total_y_start)
        {
            return 0;
        }

        fill_subregion_ul(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, uref, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);
        
        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }
    }
    else if(local_y_end == Window_Len - 1)
    {
        ///if local_y_end == Window_Len - 1 and old_extra_end > 0
        ///this means local_y_end is the end of the y
        ///so shit to the right makes no sense
        if(old_extra_end != 0)
        {
            return 0;
        }
        long long total_y_end = total_y_start + local_y_end;

        total_y_start = total_y_end - x_len + 1;

        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, uref->ug->u.a[y_ID].len,
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        if(new_y_start == total_y_end - local_y_end)
        {
            return 0;
        }

        fill_subregion_ul(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, uref, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);

        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }


    }
    return 0;
}


inline int fix_boundary(char* x_string, long long x_len, int threshold,
long long total_y_start, long long local_y_start, long long local_y_end,
long long old_extra_begin, long long old_extra_end,
long long y_ID, long long Window_Len, All_reads* R_INF,
Correct_dumy* dumy, int y_strand, unsigned int old_error,
long long* r_total_y_start, int* r_start_site, int* r_end_site,
int* r_extra_begin, int* r_extra_end, unsigned int* r_error)
{

    
    int new_extra_begin, new_extra_end;
    long long new_y_start, new_y_length;
    int new_end_site, new_start_site;
    unsigned int new_error;
    char* y_string;


    int path_length;

    ///if the start pos at the left boundary 
    if(local_y_start == 0)
    {
        total_y_start = total_y_start + local_y_start;
        ///if local_y_start == 0 and old_extra_begin != 0
        ///this means total_y_start == 0, so shift to the left cannot get a new start pos
        if(old_extra_begin != 0)
        {
            return 0;
        }

        ///if the begining of alignment is 0, we should try to shift the window to find a better result
        ///shift to the left by threshold-1 bases
        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, Get_READ_LENGTH((*R_INF), y_ID),
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        ///if new_y_start is equal to total_y_start, recalculate makes no sense
        if(new_y_start == total_y_start)
        {
            return 0;
        }

        fill_subregion(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, R_INF, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);
        
        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }
    }
    else if(local_y_end == Window_Len - 1)
    {
        ///if local_y_end == Window_Len - 1 and old_extra_end > 0
        ///this means local_y_end is the end of the y
        ///so shit to the right makes no sense
        if(old_extra_end != 0)
        {
            return 0;
        }
        long long total_y_end = total_y_start + local_y_end;

        total_y_start = total_y_end - x_len + 1;

        if(!determine_overlap_region(threshold, total_y_start, y_ID, Window_Len, Get_READ_LENGTH((*R_INF), y_ID),
        &new_extra_begin, &new_extra_end, &new_y_start, &new_y_length))
        {
            return 0;
        }

        if(new_y_start == total_y_end - local_y_end)
        {
            return 0;
        }

        fill_subregion(dumy->overlap_region_fix, new_y_start, new_y_length, y_strand, R_INF, y_ID, 
        new_extra_begin, new_extra_end);

        y_string = dumy->overlap_region_fix;

        new_end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &new_error, &new_start_site,
                    &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);

        if (new_error != (unsigned int)-1 && new_error < old_error)
        {
            (*r_total_y_start) = new_y_start;
            (*r_start_site) = new_start_site;
            (*r_end_site) = new_end_site;
            (*r_extra_begin) = new_extra_begin;
            (*r_extra_end) = new_extra_end;
            (*r_error) = new_error;

            dumy->path_length = path_length;
            memcpy(dumy->path, dumy->path_fix, path_length);
            memcpy(dumy->overlap_region, dumy->overlap_region_fix, Window_Len);
            return 1;
        }


    }
    return 0;
}

inline char *return_str_seq(char *buf, int64_t s, int64_t pri_l, uint8_t rev, hpc_t *hpc_g, const ul_idx_t *uref, int64_t id, int64_t aux_beg, int64_t aux_end)
{
    if(!hpc_g) {
        memset(buf, 'N', aux_beg);
        retrieve_u_seq(NULL, buf+aux_beg, &(uref->ug->u.a[id]), rev, s, pri_l, NULL);
        memset(buf+aux_beg+pri_l, 'N', aux_end);
        return buf;
    } else {
        char *z = hpc_str(*hpc_g, id, rev);
        if((aux_beg == 0) && (aux_end == 0)) {
            return z+s;
        } else {
            memset(buf, 'N', aux_beg);
            memcpy(buf+aux_beg, z+s, pri_l);
            memset(buf+aux_beg+pri_l, 'N', aux_end);
            return buf;
        }
    }
}

inline char *return_str_seq_exz(char *buf, int64_t s, int64_t pri_l, uint8_t rev, hpc_t *hpc_g, const ul_idx_t *uref, int64_t id)
{
    if(!hpc_g) {
        retrieve_u_seq(NULL, buf, &(uref->ug->u.a[id]), rev, s, pri_l, NULL);
        return buf;
    } else {
        return hpc_str(*hpc_g, id, rev) + s;
    }
}

///cannot use tstr in-place
inline int recal_boundary(char* qstr, char* tstr1, int64_t ql, int64_t thres,
int64_t global_ts0, int64_t local_ts0, int64_t local_te0,
int64_t aux_beg0, int64_t aux_end0, unsigned int err0, 
int64_t tid, int64_t aln_l, uint32_t rev, 
Correct_dumy* dumy, All_reads* rref, hpc_t *hpc_g, const ul_idx_t *uref,
int64_t* global_ts1, int* local_ts1, int* local_te1,
int64_t* aux_beg1, int64_t* aux_end1, unsigned int* err1)
{
    int64_t ts, t_tot_l, aux_beg, aux_end, t_pri_l, t_end; 
    char *q_string = qstr, *t_string; unsigned int error = (unsigned int)-1; 
    int r_ts = 0, path_length = 0;
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, tid);
    else if(uref) t_tot_l = uref->ug->u.a[tid].len;
    else t_tot_l = Get_READ_LENGTH((*rref), tid);

    if(local_ts0 == 0) {//left boundary
        if(aux_beg0 > 0) return 0;///shift to the left cannot get a new start pos
        ts = global_ts0;
    } else if((local_te0 + 1) == aln_l) {//right boundary
        if(aux_end0 > 0) return 0;///shift to the right cannot get a new start pos
        ts = global_ts0 + local_te0 - ql + 1;
    } else {
        return 0;
    }
    if(!init_waln(thres, ts, t_tot_l, aln_l, &aux_beg, &aux_end, &ts, &t_pri_l)) return 0;
    if(ts == global_ts0) return 0;//unchanged, make no sense

    if(rref) {
        fill_subregion(tstr1, ts, t_pri_l, rev, rref, tid, aux_beg, aux_end); t_string = tstr1;
    } else {
        t_string = return_str_seq(tstr1, ts, t_pri_l, rev, hpc_g, uref, tid, aux_beg, aux_end);
    }

    t_end = Reserve_Banded_BPM_PATH(t_string, aln_l, q_string, ql, thres, &error, &r_ts,
                &path_length, dumy->matrix_bit, dumy->path_fix, -1, -1);

    if (error != (unsigned int)-1 && error < err0) {
        (*global_ts1) = ts;
        (*local_ts1) = r_ts;
        (*local_te1) = t_end;
        (*aux_beg1) = aux_beg;
        (*aux_end1) = aux_end;
        (*err1) = error;

        dumy->path_length = path_length;
        memcpy(dumy->path, dumy->path_fix, path_length);
        // memcpy(tstr0, t_string, aln_l);
        return 1;
    }
    return 0;
}

///cannot use tstr in-place
inline int recal_boundary_exz(char* qstr, char* tstr, int64_t ql0, int64_t tl0, int64_t thres,
int64_t toff, int64_t ts0, int64_t te0, int64_t err0, 
int64_t tid, uint32_t rev, bit_extz_t *exz,
All_reads* rref, hpc_t *hpc_g, const ul_idx_t *uref, 
int64_t *ts_r, int64_t *aux_beg_r, int64_t *aux_end_r)
{
    int64_t ts, tl, t_tot_l, aux_beg, aux_end, t_pri_l, aln_l = ql0 + (thres << 1);
    char *q_string = qstr, *t_string; 
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, tid);
    else if(uref) t_tot_l = uref->ug->u.a[tid].len;
    else t_tot_l = Get_READ_LENGTH((*rref), tid);

    if(ts0 == 0) {//left boundary
        ts = toff;
    } else if((te0 + 1) == tl0) {//right boundary
        ts = toff + te0 - ql0 + 1;
    } else {
        return 0;
    }
    if(!init_waln(thres, ts, t_tot_l, aln_l, &aux_beg, &aux_end, &ts, &t_pri_l)) return 0;
    if(ts == toff && tl0 == t_pri_l) return 0;//unchanged, make no sense

    tl = t_pri_l;
    if(rref) {
        recover_UC_Read_sub_region(tstr, ts, tl, rev, rref, tid); t_string = tstr;
    } else {
        t_string = return_str_seq_exz(tstr, ts, tl, rev, hpc_g, uref, tid);
    }

    clear_align(*exz); 
    ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql0, thres, aux_beg, exz);

    if(is_align(*exz) && exz->err < err0) {
        (*aux_beg_r) = aux_beg;
        (*aux_end_r) = aux_end;
        (*ts_r) = ts;
        return 1;
    }
    return 0;
}

///cannot use tstr in-place
inline int recal_boundary_non_retrieve_exz(char* qstr, char* tstr, int64_t t_tot_l, 
int64_t ql0, int64_t tl0, int64_t thres,
int64_t toff, int64_t ts0, int64_t te0, int64_t err0, 
int64_t tid, uint32_t rev, bit_extz_t *exz,
int64_t *ts_r, int64_t *aux_beg_r, int64_t *aux_end_r)
{
    int64_t ts, tl, aux_beg, aux_end, t_pri_l, aln_l = ql0 + (thres << 1);
    char *q_string = qstr, *t_string; 

    if(ts0 == 0) {//left boundary
        ts = toff;
    } else if((te0 + 1) == tl0) {//right boundary
        ts = toff + te0 - ql0 + 1;
    } else {
        return 0;
    }
    if(!init_waln(thres, ts, t_tot_l, aln_l, &aux_beg, &aux_end, &ts, &t_pri_l)) return 0;
    if(ts == toff && tl0 == t_pri_l) return 0;//unchanged, make no sense

    tl = t_pri_l; t_string = tstr + ts;

    clear_align(*exz); 
    ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql0, thres, aux_beg, exz);

    if(is_align(*exz) && exz->err < err0) {
        (*aux_beg_r) = aux_beg;
        (*aux_end_r) = aux_end;
        (*ts_r) = ts;
        return 1;
    }
    return 0;
}

inline char *update_des_str(char *des, int64_t s, int64_t pri_l, uint8_t rev, All_reads *rref, hpc_t *hpc_g, 
    const ul_idx_t *uref, int64_t id, int64_t aux_beg, int64_t aux_end, char *src)
{
    if(src) {
        // memcpy(des, src, (pri_l+aux_beg+aux_end));
        // return des;
        return src;
    } else {
        if(rref) {
            fill_subregion(des, s, pri_l, rev, rref, id, aux_beg, aux_end); 
            return des;
        } else {
            return return_str_seq(des, s, pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
        }
    }
}

/**
void debug_scan_cigar(overlap_region* sub_list)
{
    long long i;
    int f_err, b_err, fLen, xLen;
    for (i = 0; i < (long long)sub_list->w_list_length; i++)
    {
        if(sub_list->w_list[i].y_end == -1 || sub_list->w_list[i].cigar.length == -1)
        {
            continue;
        }
        xLen = sub_list->w_list[i].x_end - sub_list->w_list[i].x_start + 1;

        scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        xLen, 1);
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        xLen, 0);

        if(b_err != sub_list->w_list[i].error || f_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "error\n");
        }

        scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        WINDOW, 1);
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        WINDOW, 0);

        if(b_err != sub_list->w_list[i].error || f_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "error\n");
        }

        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, 0, xLen-1);
        if(b_err != sub_list->w_list[i].error)
        {
            fprintf(stderr, "error\n");
        }

        fLen = xLen / 3;
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, fLen, 0);
        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, 0, fLen-1);
        if(f_err != b_err)
        {
            fprintf(stderr, "error\n");
        }

        fLen = xLen / 3;
        scan_cigar(&(sub_list->w_list[i].cigar), &f_err, fLen, 1);
        scan_cigar_interval(&(sub_list->w_list[i].cigar), &b_err, xLen-fLen, xLen-1);
        if(f_err != b_err)
        {
            fprintf(stderr, "\nerror\n");
            fprintf(stderr, "b_err: %d, f_err: %d\n",b_err, f_err);
            long long j;
            for (j = 0; j < sub_list->w_list[i].cigar.length; j++)
            {
                fprintf(stderr, "len: %d, opera: %d\n", 
                sub_list->w_list[i].cigar.C_L[j], sub_list->w_list[i].cigar.C_C[j]);
            }
        }

        // bLen = xLen / 3;
        // fLen = xLen - bLen;

        // scan_cigar(&(sub_list->w_list[i].cigar), &b_err, 
        // bLen, 1);
        // scan_cigar(&(sub_list->w_list[i].cigar), &f_err, 
        // fLen, 0);

        // if(b_err + f_err != sub_list->w_list[i].error)
        // {
        //     fprintf(stderr, "\nsub_list->w_list[i].error: %d, bLen: %d, b_err: %d, fLen: %d, f_err: %d\n",
        //     sub_list->w_list[i].error, bLen, b_err, fLen, f_err);
        //     long long j;
        //     for (j = 0; j < sub_list->w_list[i].cigar.length; j++)
        //     {
        //         fprintf(stderr, "len: %d, opera: %d\n", 
        //         sub_list->w_list[i].cigar.C_L[j], sub_list->w_list[i].cigar.C_C[j]);
        //     }
            
        // }
    }
}
**/

void calculate_boundary_cigars(overlap_region* z, All_reads* R_INF, Correct_dumy* dumy, UC_Read* g_read, double e_rate)
{
    assert(z->w_list.n > 0);
    int64_t nw = z->w_list.n;
    resize_window_list_alloc(&(z->boundary_cigars), nw - 1);
    int64_t y_id = z->y_id, y_strand = z->y_pos_strand;
    int64_t y_readLen = Get_READ_LENGTH((*R_INF), y_id);
    int64_t i, y_distance, f_err = -1, b_err = -1, m_error;
    int64_t scanLen = 10, boundaryLen = 200;
    int64_t single_sideLen = boundaryLen/2;
    int64_t force_useless_side = single_sideLen/2;
    int64_t L_useless_side, R_useless_side, alpha = 1;
    long long y_start, x_start, x_end, yLen, xLen, leftLen, rightLen, threshold, o_len;
    char *x_string = NULL, *y_string = NULL;
    int end_site, real_y_start, extra_begin, extra_end;
    unsigned int error;
    z->boundary_cigars.n = nw - 1;
    ///the (i)-th boundary between the (i)-th window and the (i+1)-th window
    ///that means it includes (the tail of (i)-th window) and (the header of (i+1)-th window)
    ///note the (i)-th boundary is calculated at the (i)-th window
    for (i = 0; i + 1 < nw; i++) {
        ///if both of the two windows are not aligned
        ///it is not necessary to calculate the boundary
        if(z->w_list.a[i].y_end == -1 || z->w_list.a[i+1].y_end == -1) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }
        ///y_distance can be less than 0, or larger than 0
        y_distance = (int64_t)z->w_list.a[i+1].y_start - (int64_t)z->w_list.a[i].y_end - 1;

        ///if two windows are aligned
        if(z->w_list.a[i].y_end != -1 && z->w_list.a[i+1].y_end != -1 && y_distance == 0) {
            ///scan backward
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, scanLen, 1);
            ///scan forward
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, scanLen, 0);
            if(b_err == 0 && f_err == 0) {
                z->boundary_cigars.a[i].error = -2; z->boundary_cigars.a[i].y_end = -1;
                continue;
            }
        }

        
        if(z->w_list.a[i].y_end != -1) {
            y_start = z->w_list.a[i].y_end; x_start = z->w_list.a[i].x_end;
        }///if the (i)-th window is not matched, have a look at the (i+1)-th window
        else if(z->w_list.a[i+1].y_end != -1) {
            y_start = z->w_list.a[i+1].y_start; x_start = z->w_list.a[i+1].x_start;
        }///if both of these two windows are not matched, directly skip
        else {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        ///it seems we don't need to record x_start and y_start
        z->boundary_cigars.a[i].extra_begin = x_start;
        z->boundary_cigars.a[i].extra_end = y_start;

        ///leftLen and rightLen are used for x        
        ///x should be at [sub_list->w_list[i].x_start, sub_list->w_list[i+1].x_end]
        ///y shouldn't have limitation
        ///note that the x_start and x_end should not be -1 in any case
        ///up to now, x_start and y_start are not -1
        ///leftLen does not include x_start itself, rightLen does
        ///gnerally speaking, rightLen should be always larger than leftLen
        leftLen = MIN(MIN((x_start - (int64_t)z->w_list.a[i].x_start), y_start), single_sideLen);
        rightLen = MIN(MIN(((int64_t)z->w_list.a[i+1].x_end + 1 - x_start), y_readLen - y_start), single_sideLen);

        ///xLen should be the sum length of two windows
        xLen = leftLen + rightLen;
        x_start = x_start - leftLen;
        x_end = x_start + xLen - 1;
        y_start = y_start - leftLen;

        ///if we don't have enough leftLen and rightLen
        // if(leftLen <= useless_side || rightLen <= useless_side)
        // {
        //     sub_list->boundary_cigars.buffer[i].error = -1;
        //     sub_list->boundary_cigars.buffer[i].y_end = -1;
        //     continue;
        // }


        threshold = xLen * e_rate/**asm_opt.max_ov_diff_ec**/;
        threshold = Adjust_Threshold(threshold, xLen);
        threshold = double_error_threshold(threshold, xLen);

        yLen = xLen + (threshold << 1);
        if(!determine_overlap_region(threshold, y_start, y_id, yLen, Get_READ_LENGTH((*R_INF), y_id), 
                    &extra_begin, &extra_end, &y_start, &o_len)) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        if(o_len < xLen) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, R_INF, y_id, extra_begin, extra_end);
        
        x_string = g_read->seq + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM_PATH(y_string, yLen, x_string, xLen, threshold, &error, 
        &real_y_start, &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

        ///means this window is matched
        if (error!=(unsigned int)-1) {
            z->boundary_cigars.a[i].x_start = x_start;
            z->boundary_cigars.a[i].x_end = x_end;

            generate_cigar(dumy->path, dumy->path_length, &(z->boundary_cigars.a[i]), &(z->boundary_cigars),
                        &real_y_start, &end_site, &error, x_string, xLen, y_string);
            ///should not adjust cigar here, adjust cigar may cause problem
            ///that is not what we want

            ///y_distance can be less than 0, or larger than 0
            ///please if one of the two windows is not matched, 
            ///y_distance may have potential problems 
            if(y_distance < 0) y_distance = y_distance * (-1);
            ///leftLen, rightLen
            // if(leftLen <= useless_side || rightLen <= useless_side)
            // {
            //     sub_list->boundary_cigars.buffer[i].error = -1;
            //     sub_list->boundary_cigars.buffer[i].y_end = -1;
            //     continue;
            // }
            L_useless_side = R_useless_side = force_useless_side;

            ///first window
            if((i == 0) && (x_start == (int64_t)z->w_list.a[0].x_start)) {
                L_useless_side = 0;
            }
            ///last window
            if((i == (int64_t)(z->w_list.n) - 2) && 
                    (x_end == (long long)(z->w_list.a[(int64_t)(z->w_list.n)-1].x_end))) {
                R_useless_side = 0;
            }

            if(leftLen <= L_useless_side || rightLen <= R_useless_side) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }



            ///up to now, if we require (i)-th window and (i+1)-th window are matched
            ///boundary_cigars.buffer[i].cigar, w_list[i].cigar and w_list[i+1].cigar are avaiable
            ///get the error excluding the first and the last useless_side bases
            scan_cigar_interval(&(z->boundary_cigars.a[i]), &(z->boundary_cigars), &m_error, L_useless_side, xLen-R_useless_side-1);
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, leftLen-L_useless_side, 1);
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, rightLen-R_useless_side, 0);

            if(f_err + b_err + y_distance + alpha < m_error) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }

            z->boundary_cigars.a[i].error = error;
            z->boundary_cigars.a[i].y_start = y_start + real_y_start - extra_begin;
            z->boundary_cigars.a[i].y_end = y_start + end_site - extra_begin;

            z->boundary_cigars.a[i].x_start = x_start;
            z->boundary_cigars.a[i].x_end = x_end;
            ///sub_list->boundary_cigars.buffer[i].error_threshold = useless_side;
            z->boundary_cigars.a[i].extra_begin = L_useless_side;
            z->boundary_cigars.a[i].extra_end = R_useless_side;
        }
        else {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }     
    }

}


void calculate_ul_boundary_cigars(overlap_region* z, const ul_idx_t *uref, Correct_dumy* dumy, 
UC_Read* g_read, double max_ov_diff_ec, long long blockLen)
{
    assert(z->w_list.n > 0);
    int64_t nw = z->w_list.n;
    resize_window_list_alloc(&(z->boundary_cigars), nw - 1);
    int64_t y_id = z->y_id;
    int64_t y_strand = z->y_pos_strand;
    int64_t y_readLen = uref->ug->u.a[y_id].len;
    int64_t i, y_distance;
    int64_t f_err, b_err, m_error, scanLen = 10;
    int64_t boundaryLen = WINDOW_UL_BOUND_RATE*blockLen;
    boundaryLen >>= 2; boundaryLen <<= 2;
    if(boundaryLen < WINDOW_UL_BOUND) boundaryLen = WINDOW_UL_BOUND;
    int64_t single_sideLen = boundaryLen/2;
    int64_t force_useless_side = single_sideLen/2;
    int64_t L_useless_side, R_useless_side;
    int64_t alpha = 1;
    long long y_start, x_start, x_end, yLen, xLen, leftLen, rightLen, threshold, o_len;
    int extra_begin, extra_end, end_site, real_y_start;
    char* x_string;
    char* y_string;
    unsigned int error;
    z->boundary_cigars.n = nw - 1;
    ///the (i)-th boundary between the (i)-th window and the (i+1)-th window
    ///that means it includes (the tail of (i)-th window) and (the header of (i+1)-th window)
    ///note the (i)-th boundary is calculated at the (i)-th window
    for (i = 0; i + 1 < nw; i++) {
        ///if both of the two windows are not aligned
        ///it is not necessary to calculate the boundary
        ///if(sub_list->w_list[i].y_end == -1 && sub_list->w_list[i+1].y_end == -1)
        if(z->w_list.a[i].y_end == -1 || z->w_list.a[i+1].y_end == -1) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        ///note if w_list[i+1].y_start or sub_list->w_list[i].y_end is -1
        ///y_distance might have some problems at the last of this function
        ///we need to deal with it carefully
        y_distance = (int64_t)z->w_list.a[i+1].y_start - (int64_t)z->w_list.a[i].y_end - 1;

        ///if two windows are aligned
        if(z->w_list.a[i].y_end != -1 && z->w_list.a[i+1].y_end != -1 && y_distance == 0) {
            ///scan backward
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, scanLen, 1);
            ///scan forward
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, scanLen, 0);
            if(b_err == 0 && f_err == 0) {
                z->boundary_cigars.a[i].error = -2; z->boundary_cigars.a[i].y_end = -1;
                continue;
            }
        }


        ///y_distance can be less than 0, or larger than 0
        if(z->w_list.a[i].y_end != -1) {
            y_start = z->w_list.a[i].y_end; x_start = z->w_list.a[i].x_end;
        }///if the (i)-th window is not matched, have a look at the (i+1)-th window
        else if(z->w_list.a[i+1].y_end != -1) {
            y_start = z->w_list.a[i+1].y_start; x_start = z->w_list.a[i+1].x_start;
        }///if both of these two windows are not matched, directly skip
        else {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        ///it seems we don't need to record x_start and y_start
        z->boundary_cigars.a[i].extra_begin = x_start;
        z->boundary_cigars.a[i].extra_end = y_start;

        ///leftLen and rightLen are used for x        
        ///x should be at [sub_list->w_list[i].x_start, sub_list->w_list[i+1].x_end]
        ///y shouldn't have limitation
        ///note that the x_start and x_end should not be -1 in any case
        ///up to now, x_start and y_start are not -1
        ///leftLen does not include x_start itself, rightLen does
        ///gnerally speaking, rightLen should be always larger than leftLen
        leftLen = MIN(MIN((x_start - (long long)z->w_list.a[i].x_start), y_start), single_sideLen);
        rightLen = MIN(MIN(((long long)z->w_list.a[i+1].x_end + 1 - x_start), y_readLen - y_start), single_sideLen);

        ///xLen should be the sum length of two windows
        xLen = leftLen + rightLen;
        x_start = x_start - leftLen;
        x_end = x_start + xLen - 1;
        y_start = y_start - leftLen;

        ///if we don't have enough leftLen and rightLen
        // if(leftLen <= useless_side || rightLen <= useless_side)
        // {
        //     sub_list->boundary_cigars.buffer[i].error = -1;
        //     sub_list->boundary_cigars.buffer[i].y_end = -1;
        //     continue;
        // }


        threshold = xLen * max_ov_diff_ec;
        threshold = Adjust_Threshold(threshold, xLen);
        threshold = double_ul_error_threshold(threshold, xLen);

        yLen = xLen + (threshold << 1);
        if(!determine_overlap_region(threshold, y_start, y_id, yLen, uref->ug->u.a[y_id].len, 
                    &extra_begin, &extra_end, &y_start, &o_len)) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        if(o_len < xLen) {
            z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
            continue;
        }

        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, 
                    uref, y_id, extra_begin, extra_end);
        
        x_string = g_read->seq + x_start;
        y_string = dumy->overlap_region;

        end_site = Reserve_Banded_BPM_PATH(y_string, yLen, x_string, xLen, threshold, &error, 
        &real_y_start, &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

        ///means this window is matched
        if (error!=(unsigned int)-1) {
            z->boundary_cigars.a[i].x_start = x_start; z->boundary_cigars.a[i].x_end = x_end;

            generate_cigar(dumy->path, dumy->path_length, &(z->boundary_cigars.a[i]), &(z->boundary_cigars),
                        &real_y_start, &end_site, &error, x_string, xLen, y_string);
            ///should not adjust cigar here, adjust cigar may cause problem
            ///that is not what we want

            ///y_distance can be less than 0, or larger than 0
            ///please if one of the two windows is not matched, 
            ///y_distance may have potential problems 
            if(y_distance < 0) y_distance = y_distance * (-1);
            ///leftLen, rightLen
            // if(leftLen <= useless_side || rightLen <= useless_side)
            // {
            //     sub_list->boundary_cigars.buffer[i].error = -1;
            //     sub_list->boundary_cigars.buffer[i].y_end = -1;
            //     continue;
            // }
            L_useless_side = R_useless_side = force_useless_side;

            ///first window
            if((i == 0) && (x_start == (long long)z->w_list.a[0].x_start)) {
                L_useless_side = 0;
            }
            ///last window
            if((i == (int64_t)(z->w_list.n) - 2) && (x_end == (z->w_list.a[(int64_t)z->w_list.n - 1].x_end))) {
                R_useless_side = 0;
            }

            if(leftLen <= L_useless_side || rightLen <= R_useless_side) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }



            ///up to now, if we require (i)-th window and (i+1)-th window are matched
            ///boundary_cigars.buffer[i].cigar, w_list[i].cigar and w_list[i+1].cigar are avaiable
            ///get the error excluding the first and the last useless_side bases
            scan_cigar_interval(&(z->boundary_cigars.a[i]), &(z->boundary_cigars), &m_error, 
            L_useless_side, xLen-R_useless_side-1);
            scan_cigar(&(z->w_list.a[i]), &(z->w_list), &b_err, leftLen-L_useless_side, 1);
            scan_cigar(&(z->w_list.a[i+1]), &(z->w_list), &f_err, rightLen-R_useless_side, 0);

            if(f_err + b_err + y_distance + alpha < m_error) {
                z->boundary_cigars.a[i].error = -1; z->boundary_cigars.a[i].y_end = -1;
                z->boundary_cigars.c.n = z->boundary_cigars.a[i].cidx;
                continue;
            }

            z->boundary_cigars.a[i].error = error;
            z->boundary_cigars.a[i].y_start = y_start + real_y_start - extra_begin;
            z->boundary_cigars.a[i].y_end = y_start + end_site - extra_begin;

            z->boundary_cigars.a[i].x_start = x_start;
            z->boundary_cigars.a[i].x_end = x_end;
            ///sub_list->boundary_cigars.buffer[i].error_threshold = useless_side;
            z->boundary_cigars.a[i].extra_begin = L_useless_side;
            z->boundary_cigars.a[i].extra_end = R_useless_side;
        }
        else
        {
            z->boundary_cigars.a[i].error = -1;
            z->boundary_cigars.a[i].y_end = -1;
            continue;
        }
         
    }

}


/**
void debug_window_cigar(overlap_region_alloc* overlap_list, UC_Read* g_read, Correct_dumy* dumy,
All_reads* R_INF, int test_window, int test_boundary)
{
    uint64_t i, j, y_id, y_strand;
    char* x_string;
    char* y_string;
    long long x_start;
    long long x_end;
    long long x_len;
    long long y_start;
    long long y_end;
    long long y_len;

    for (j = 0; j < overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        if(overlap_list->list[j].is_match == 1)
        {

            if(test_window == 1)
            {
                for (i = 0; i < overlap_list->list[j].w_list_length; i++)
                {
                    if(overlap_list->list[j].w_list[i].y_end != -1)
                    {
                        ///there is no problem for x
                        x_start = overlap_list->list[j].w_list[i].x_start;
                        x_end = overlap_list->list[j].w_list[i].x_end;
                        x_len = x_end - x_start + 1;

                        x_string = g_read->seq + x_start;

                        y_start = overlap_list->list[j].w_list[i].y_start;
                        y_end = overlap_list->list[j].w_list[i].y_end;
                        y_len = y_end - y_start + 1;

                        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, R_INF, y_id);
                        y_string = dumy->overlap_region;


                        if(verify_cigar(x_string, x_len, y_string, y_len, &overlap_list->list[j].w_list[i].cigar, 
                        overlap_list->list[j].w_list[i].error))
                        {
                            fprintf(stderr, "error\n");
                        }
                    }
                }
            }
            

            
            if(test_boundary == 1)
            {
                for (i = 0; i < (uint64_t)overlap_list->list[j].boundary_cigars.length; i++)
                {
                    if(overlap_list->list[j].boundary_cigars.buffer[i].y_end != -1)
                    {
                        x_start = overlap_list->list[j].boundary_cigars.buffer[i].x_start;
                        x_end = overlap_list->list[j].boundary_cigars.buffer[i].x_end;
                        x_len = x_end - x_start + 1;
                        x_string = g_read->seq + x_start;

                        y_start = overlap_list->list[j].boundary_cigars.buffer[i].y_start;
                        y_end = overlap_list->list[j].boundary_cigars.buffer[i].y_end;
                        y_len = y_end - y_start + 1;

                        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_len, y_strand, 
                        R_INF, y_id);
                        y_string = dumy->overlap_region;


                        if(verify_cigar(x_string, x_len, y_string, y_len, 
                        &overlap_list->list[j].boundary_cigars.buffer[i].cigar, 
                        overlap_list->list[j].boundary_cigars.buffer[i].error))
                        {
                            fprintf(stderr, "error\n");
                        }
                    }
                }
            }


            if(test_window == 1 && test_boundary == 1)
            {
                if(overlap_list->list[j].w_list_length != 
                        (uint64_t)(overlap_list->list[j].boundary_cigars.length + 1))
                {
                    fprintf(stderr, "error\n");
                }
            }
            
        }
    }
}
**/
int64_t get_adjust_winid(overlap_region *z, int64_t win_beg, int64_t win_len)
{   
    int64_t win_id, k;
    win_id = (win_beg-((z->x_pos_s/win_len)*win_len))/win_len;
    if((uint64_t)win_id < z->w_list.n && z->w_list.a[win_id].x_start == win_beg) return win_id;
    if(z->w_list.n == 0) return -1;
    // if(z->w_list.a[win_id].x_start <= win_beg) {
    //     fprintf(stderr, "z->w_list.n::%u, z->w_list.a[%ld].x_start::%d, win_beg::%ld\n", 
    //     (uint32_t)z->w_list.n, win_id, z->w_list.a[win_id].x_start, win_beg);
    // }
    if((uint64_t)win_id > z->w_list.n) win_id = z->w_list.n;
    // assert((z->w_list.a[win_id].x_start > win_beg);
    for (k = win_id - 1; k >= 0; k--) {
        // if(k < 0 || k >= (int64_t)z->w_list.n) fprintf(stderr, "win_id::%ld, k::%ld, z->w_list.n::%ld\n", win_id, k, (int64_t)z->w_list.n);
        if(z->w_list.a[k].x_start == win_beg) return k;
        if(z->w_list.a[k].x_start < win_beg) return -1;
    }
    return -1;
}

void set_herror_win(overlap_region_alloc* ovlp, Correct_dumy* du, kvec_t_u64_warp* v_idx, double max_ov_diff_ec, int64_t rLen, int64_t blockLen)
{
    Window_Pool w_inf; int32_t flag = 0; uint64_t cID, mm, fc, fw, idx_n, idx_i;
    init_Window_Pool(&w_inf, rLen, blockLen, (int)(1.0/max_ov_diff_ec));
    long long window_start, window_end; int64_t i, k, mLen, w_list_id, ws, we; 

    idx_n = get_num_wins(0, rLen, blockLen); idx_i = 0;
    kv_resize(uint64_t, v_idx->a, idx_n); v_idx->a.n = idx_n; 

    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2) {
        du->length = du->lengthNT = 0;
        flag = get_interval(window_start, window_end, ovlp, du, w_inf.window_length);
        switch (flag) {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        v_idx->a.a[idx_i++] = ((uint64_t)(v_idx->a.n))<<32;
        for (i = 0; i < (int64_t)du->length; i++) {
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            w_list_id = get_adjust_winid(&(ovlp->list[cID]), window_start, w_inf.window_length);
            if(w_list_id >= 0) break;///a matched window
        }
        if(i < (int64_t)du->length) continue;///if there is a matched window

        for (i = 0, mm = 0; i < (int64_t)du->length; i++) {///all windows are unmatched
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            ovlp->list[cID].is_match = 4; 
            ovlp->list[cID].align_length += window_end + 1 - window_start;
            mm++;
        }
        if(mm > 0) {
            kv_push(uint64_t, v_idx->a, (((uint64_t)window_start)<<32)|((uint64_t)window_end));
            v_idx->a.a[idx_i-1]++;
        }


        ///shorter than blockLen
        for (i = du->size-du->lengthNT, mLen = du->size-du->lengthNT, fc = 0; i < (int64_t)du->size; i++) {
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            get_win_se_by_normalize_xs(&(ovlp->list[cID]), window_start, blockLen, &ws, &we);
            w_list_id = get_adjust_winid(&(ovlp->list[cID]), ws, blockLen);
            if (w_list_id >= 0) {///matched
                cID = w_list_id; cID <<= 32; cID += (uint32_t)du->overlapID[i]; du->overlapID[i] = cID;
                if(mLen != i) {
                    mm = du->overlapID[i]; du->overlapID[i] = du->overlapID[mLen]; du->overlapID[mLen] = mm;
                }
                mLen++;
            } else {///unmatched
                cID = (uint32_t)-1; cID <<= 32; cID += (uint32_t)du->overlapID[i]; du->overlapID[i] = cID;
                fc++;
            }
        }
        // if(mLen == (int64_t)du->size) continue;///if all windows shorter than blockLen are matched
        if(fc == 0) continue;///no unmatched windows that are shorter than blockLen
        for (i = mLen; i < (int64_t)du->size; i++){///check the remaining unmatched windows that are shorter than blockLen
            cID = (uint32_t)du->overlapID[i];
            if(ovlp->list[cID].is_match!=3 && ovlp->list[cID].is_match!=4) continue;
            assert((du->overlapID[i]>>32)==(uint32_t)-1);
            get_win_se_by_normalize_xs(&(ovlp->list[cID]), window_start, blockLen, &ws, &we);
            for (k = du->size-du->lengthNT; k < mLen; k++) {///all matched windows
                fc = (uint32_t)du->overlapID[k]; fw = du->overlapID[k]>>32;
                assert(fw!=(uint32_t)-1); assert(ovlp->list[fc].is_match == 3 || ovlp->list[fc].is_match == 4);
                // if (ovlp->list[fc].w_list[fw].y_end == -1 || (ovlp->list[fc].is_match!=3 && ovlp->list[fc].is_match!=4)) fprintf(stderr, "ERROR\n");        
                ///if there is one matched window can cover the unmatched window
                if(ovlp->list[fc].w_list.a[fw].x_start<=ws && ovlp->list[fc].w_list.a[fw].x_end>=we) { 
                        break;
                }
            }

            if(k >= mLen) {///no matched window can cover the unmatched window
                ovlp->list[cID].is_match = 4; 
                ovlp->list[cID].align_length += we + 1 - ws;
                kv_push(uint64_t, v_idx->a, (((uint64_t)ws)<<32)|((uint64_t)we));
                v_idx->a.a[idx_i-1]++;
            }
        }
    }
}


inline void recalcate_window_advance(overlap_region_alloc* overlap_list, All_reads *rref, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, double e_rate_final)
{
    long long j, k, i;
    int threshold;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;
    int64_t nw, a_nw, w_id, w_s, w_e, is_srt;
    double error_rate;
    uint64_t *w_idx;
    overlap_region *z;
    window_list *p = NULL;


    overlap_list->mapped_overlaps_length = 0;
    for (j = 0; j < (long long)overlap_list->length; j++) {
        z = &(overlap_list->list[j]); z->is_match = 0; is_srt = 1;
        if(z->w_list.n == 0) continue;///no alignment
        nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s); a_nw = z->w_list.n;
        kv_resize(uint64_t, v_idx->a, (uint64_t)nw); memset(v_idx->a.a, -1, sizeof((*v_idx->a.a))*nw); w_idx = v_idx->a.a;
        for (i = 0; i < a_nw; i++) {
            assert(z->w_list.a[i].y_end != -1);
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, NULL);
            w_idx[w_id] = i;
        }
        // if(j == 248) {
        //     fprintf(stderr, "0-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu\n", __func__, 
        //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0]);
        // }

        y_id = z->y_id; y_strand = z->y_pos_strand; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        for (i = a_nw-1; i >= 0; i--) { //utilize the the end pos of pre-window in forward
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, &w_e);
            // if(z->w_list.a[i].x_end != w_e) {
            //     fprintf(stderr, "[M::%s] block_s->%ld, w_id->%ld, z::x_pos_s->%u, z::x_pos_e->%u, x_start->%d, x_end->%d, w_e->%ld\n", __func__, block_s, w_id, z->x_pos_s, z->x_pos_e,
            //     z->w_list.a[i].x_start, z->w_list.a[i].x_end, w_e);
            // }
            assert(z->w_list.a[i].x_end == w_e);
            total_y_start = z->w_list.a[i].y_end + 1 - z->w_list.a[i].extra_begin;
            for (k = w_id + 1; k < nw; k++) {
                if(w_idx[k] != (uint64_t)-1) break;
                w_s = w_e + 1;
                w_id = get_win_id_by_s(z, w_s, block_s, &w_e);
                assert(w_id == k);
                extra_begin = extra_end = 0;
                if (total_y_start >= y_readLen) break;
                ///there is no problem for x
                x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start; y_start = total_y_start; 
                ///there are two potiential reasons for unmatched window:
                ///1. this window has a large number of differences
                ///2. DP does not start from the right offset
                if(rref) {
                    threshold = double_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD), x_len);
                } else {
                    threshold = double_ul_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD_MAX_SIZE), x_len);
                }
                
                Window_Len = x_len + (threshold << 1);

                if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len)), 
                &extra_begin, &extra_end, &y_start, &o_len)) {
                    break;
                }
                if(o_len + threshold < x_len) break;
                
                if(rref) {
                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                } else {
                    fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                }

                x_string = g_read->seq + x_start; y_string = dumy->overlap_region;
                ///note!!! need notification
                end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);                    
                if (error!=(unsigned int)-1) {///unmatched
                    kv_pushp(window_list, z->w_list, &p);
                    p->x_start = x_start;
                    p->x_end = x_end;
                    p->y_start = y_start;
                    p->y_end = y_start + end_site;
                    p->error = error;
                    p->extra_begin = extra_begin;
                    p->extra_end = extra_end;
                    p->error_threshold = threshold;
                    p->cidx = p->clen = 0;

                    z->align_length += x_len; w_idx[k] = z->w_list.n - 1;

                    if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                }
                else {
                    break;
                }

                total_y_start = y_start + end_site + 1 - extra_begin;
            }
        }
        // if(j == 248) {
        //     fprintf(stderr, "1-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu\n", __func__, 
        //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0]);
        // }
        for (i = 0; i < nw; i++) { //utilize the the start pos of next window in backward
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(i > 0 && w_idx[i] != (uint64_t)-1 && w_idx[i-1] == (uint64_t)-1) {
                w_s = z->w_list.a[w_idx[i]].x_start;
                ///check if the start pos of this matched window has been calculated
                if(z->w_list.a[w_idx[i]].clen == 0) {
                    p = &(z->w_list.a[w_idx[i]]);
                    ///there is no problem for x
                    x_start = p->x_start; x_end = p->x_end; x_len = x_end + 1 - x_start; threshold = p->error_threshold;
                    /****************************may have bugs********************************/
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
                    Window_Len = x_len + (threshold << 1);
                    ///y_start is the real y_start
                    y_start = p->y_start; extra_begin = p->extra_begin; extra_end = p->extra_end;
                    o_len = Window_Len - extra_end - extra_begin;
                    if(rref) {
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                    } else {
                        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                    }
                    
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - y_start);
                    assert(error != (unsigned int)-1);

                    {
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0) {
                            if(rref) {
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, rref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            } else {
                                if(fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            }
                        }
                                                 
                        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);   

                        ///note!!! need notification
                        real_y_start = y_start + real_y_start - extra_begin;
                        p->y_start = real_y_start; 
                        ///I forget why don't reduce the extra_begin for y_end
                        ///it seems extra_begin will be reduced at the end of this function 
                        p->y_end = y_start + end_site;
                        p->error = error;                         
                    }
                } else {
                    real_y_start = p->y_start;
                }

                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && w_idx[k] == (uint64_t)-1; k--) {  
                    w_e = w_s - 1;
                    w_id = get_win_id_by_e(z, w_e, block_s, &w_s);
                    assert(w_id == k);
                    ///there is no problem in x
                    x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
                    ///there are two potiential reasons for unmatched window:
                    ///1. this window has a large number of differences
                    ///2. DP does not start from the right offset
                    if(rref) {
                        threshold = double_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD), x_len);
                    } else {
                        threshold = double_ul_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD_MAX_SIZE), x_len);
                    }
                    Window_Len = x_len + (threshold << 1);
                    if(total_y_end <= 0) break;
                    
                    ///y_start might be less than 0
                    y_start = total_y_end - x_len + 1;
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len)), 
                    &extra_begin, &extra_end, &y_start, &o_len)) {
                        break;
                    }

                    if(o_len + threshold < x_len) break;
                    
                    if(rref) {
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                    } else {
                        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                    }
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

                    if (error!=(unsigned int)-1) { 
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0) {
                            if(rref) {
                                fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, rref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error);
                            } else {
                                fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error);
                            }
                        }

                        kv_pushp(window_list, z->w_list, &p);
                        p->x_start = x_start; p->x_end = x_end;///must set x_start/x_end here
                        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);  
                        
                        ///y_start has no shift, but y_end has shift               
                        p->y_start = y_start + real_y_start - extra_begin;
                        p->y_end = y_start + end_site;
                        p->error = error;
                        p->extra_begin = extra_begin;
                        p->extra_end = extra_end;
                        p->error_threshold = threshold;
                        z->align_length += x_len; w_idx[k] = z->w_list.n - 1;

                        if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                    }
                    else {
                        break;
                    }

                    total_y_end = y_start + real_y_start - 1 - extra_begin;
                }
            }
        }

        // if(j == 248) {
        //     fprintf(stderr, "2-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu, w_idx[0]->cidx::%u, w_idx[0]->clen::%u, w_idx[0]->cigar[0]:%u\n", __func__, 
        //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0], z->w_list.a[w_idx[0]].cidx, z->w_list.a[w_idx[0]].clen, z->w_list.c.a[z->w_list.a[w_idx[0]].cidx]);
        // }

        if(uref) {
            z->is_match = 0;
            if((((z->x_pos_e + 1 - z->x_pos_s)*MIN_UL_ALIN_RATE) <= z->align_length) && (z->align_length >= MIN_UL_ALIN_LEN)){
                z->is_match = 3; overlap_list->mapped_overlaps_length += z->align_length;
                ///sort for set_herror_win
                if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            }
        }
    }

    // fprintf(stderr, "+++[M::%s::idx->%d::y_id->%u] z::align_length->%u, e_threshold->%f\n", 
    //         __func__, 27, overlap_list->list[27].y_id, overlap_list->list[27].align_length, e_rate);
    
    // fprintf(stderr, "+++[M::%s::idx->%d::y_id->%u] z::align_length->%u, e_threshold->%f\n", 
    //         __func__, 45, overlap_list->list[45].y_id, overlap_list->list[45].align_length, e_rate);

    // fprintf(stderr, "+++[M::%s::idx->%d::y_id->%u] z::align_length->%u, e_threshold->%f\n", 
    //         __func__, 277, overlap_list->list[277].y_id, overlap_list->list[277].align_length, e_rate);

    if(uref && overlap_list->mapped_overlaps_length > 0) {
        set_herror_win(overlap_list, dumy, v_idx, e_rate, g_read->length, block_s);
    }

    overlap_list->mapped_overlaps_length = 0;
    for (j = 0; j < (long long)overlap_list->length; j++) {
        z = &(overlap_list->list[j]);
        y_id = z->y_id; y_strand = z->y_pos_strand; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        overlap_length = z->x_pos_e + 1 - z->x_pos_s; //z->is_match = 0;
        // if(y_id == 0 || y_id == 1) {
        //     fprintf(stderr, "[M::%s::j->%lld] utg%.6dl(%c), align_length::%u, overlap_length::%lld\n", __func__, 
        //         j, (int32_t)z->y_id + 1, "+-"[z->y_pos_strand], z->align_length, overlap_length);
        // }
        // if(y_id == 24128) {
        //     fprintf(stderr, "[M::%s::idx->%lld::] x::[%u, %u), y::[%u, %u), ovl->%lld, aln->%u\n", 
        //     __func__, j, z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1, overlap_length, z->align_length);
        //     radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n); a_nw = z->w_list.n; 
        //     int64_t ss, ee;
        //     for (i = 0, ss = ee = -2; i < a_nw; i++) {
        //         p = &(z->w_list.a[i]);
        //         if(p->x_start == ee) {
        //             ee = p->x_end + 1;
        //         } else {
        //             if(ee > 0) {
        //                 fprintf(stderr, "[M::%s::] x::[%ld, %ld)\n", __func__, ss, ee);
        //             }
        //             ss = p->x_start; ee = p->x_end + 1;
        //         }
        //     }
        //     if(ee > 0) {
        //         fprintf(stderr, "[M::%s::] x::[%ld, %ld)\n", __func__, ss, ee);
        //     }

        //     for (i = 0; i < a_nw; i++) {
        //         p = &(z->w_list.a[i]); if(p->y_end == -1) continue;
        //         fprintf(stderr, "[M::%s::] x::[%d, %d), y::[%d, %d), error::%d\n", __func__, p->x_start, p->x_end + 1, p->y_start, p->y_end + 1, p->error);
        //     }
        // }

        ///debug_scan_cigar(&(overlap_list->list[j]));
        ///only calculate cigar for high quality overlaps
        if ((rref && (overlap_length*OVERLAP_THRESHOLD_HIFI_FILTER <= z->align_length)) || 
                                            (uref && (overlap_length*(1-e_rate) <= z->align_length))) {
            a_nw = z->w_list.n; 
            // int64_t tt = 0;
            for (i = 0, is_srt = 1; i < a_nw; i++) {
                p = &(z->w_list.a[i]);
                ///check if the cigar of this window has been got 
                if(p->clen == 0) {
                    ///there is no problem for x
                    x_start = p->x_start; x_end = p->x_end; x_len = x_end - x_start + 1;
                    /****************************may have bugs********************************/
                    ///threshold = x_len * asm_opt.max_ov_diff_ec;
                    threshold = p->error_threshold;
                    /****************************may have bugs********************************/
                    /****************************may have bugs********************************/
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    /****************************may have bugs********************************/
                    Window_Len = x_len + (threshold << 1);


                    ///y_start is the real y_start
                    ///for the window with cigar, y_start has already reduced extra_begin
                    y_start = p->y_start; extra_begin = p->extra_begin; extra_end = p->extra_end;
                    o_len = Window_Len - extra_end - extra_begin;
                    if(rref) {
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
                    } else {
                        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
                    }
                    x_string = g_read->seq + x_start; y_string = dumy->overlap_region;


                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - y_start);
                    // if(!(error != (unsigned int)-1)) {
                    //     fprintf(stderr, "[M::%s::]\tqid::%u\tqlen::%lu\tq::[%d,\t%d)\ttid::%u\ttlen::%lu\tt::[%d,\t%d)\te_beg::%d\te_end::%d\terr::%d\n", __func__, 
                    //     overlap_list->list[j].x_id, Get_READ_LENGTH((*rref), overlap_list->list[j].x_id), 
                    //     p->x_start, p->x_end+1, 
                    //     overlap_list->list[j].y_id, Get_READ_LENGTH((*rref), overlap_list->list[j].y_id), 
                    //     p->y_start, p->y_end+1, 
                    //     p->extra_begin, p->extra_end, p->error);
                    //     fprintf(stderr, "[M::%s::]\tqcal_len::%lld\ttcal_len::%lld\tthres::%d\n", __func__, 
                    //     x_len, Window_Len, threshold);
                    //     fprintf(stderr, "qid::%u\nqname::%.*s\n\t%.*s\n", overlap_list->list[j].x_id, 
                    //     (int32_t)Get_NAME_LENGTH((*rref), overlap_list->list[j].x_id), 
                    //     Get_NAME((*rref), overlap_list->list[j].x_id), (int32_t)x_len, x_string);

                    //     fprintf(stderr, "tid::%u\ntname::%.*s\n\t%.*s\n", overlap_list->list[j].y_id, 
                    //     (int32_t)Get_NAME_LENGTH((*rref), overlap_list->list[j].y_id), 
                    //     Get_NAME((*rref), overlap_list->list[j].y_id), (int32_t)Window_Len, y_string);

                    // }
                    assert(error != (unsigned int)-1);

                    {
                        if (end_site == Window_Len - 1 || real_y_start == 0) {
                            if(rref) {
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                                extra_begin, extra_end, y_id, Window_Len, rref, dumy, y_strand, error,
                                &y_start, &real_y_start, &end_site, &extra_begin, &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            } else {
                                if(fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                                &extra_end, &error)) {
                                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                                }
                            }
                        }

                        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);    

                        ///note!!! need notification
                        real_y_start = y_start + real_y_start - extra_begin;
                        p->y_start = real_y_start;  
                        p->y_end = y_start + end_site - extra_begin;
                        p->error = error;                              
                    }

                    // if(y_id == 4) {
                    //     fprintf(stderr, "+[M::idx->%lld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                else {
                    p->y_end -= p->extra_begin;
                    // if(y_id == 4) {
                    //     fprintf(stderr, "-[M::idx->%lld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                // tt += p->error;
                if(is_srt && i > 0 && p->x_start < z->w_list.a[i-1].x_start) is_srt = 0;
            }

            if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            error_rate = non_trim_error_rate(z, rref, uref, v_idx, dumy, g_read, e_rate, block_s);
            z->is_match = 0;
            // if(y_id == 4) {
            //     fprintf(stderr, "[M::%s::idx->%lld::] z::x_pos_s->%u, z::x_pos_e->%u, ovl->%lld, aln->%u, error_rate->%f, e_rate_final->%f, tt->%ld\n", 
            //     __func__, j, z->x_pos_s, z->x_pos_e, overlap_length, z->align_length, error_rate, e_rate_final, tt);
            // }

            if (error_rate <= e_rate_final/**asm_opt.max_ov_diff_final**/) {
                overlap_list->mapped_overlaps_length += overlap_length;
                z->is_match = 1; append_unmatched_wins(z, block_s);
                // if(j == 248) {
                //     fprintf(stderr, "3-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu, w_idx[0]->cidx::%u, w_idx[0]->clen::%u, w_idx[0]->cigar[0]:%u\n", __func__, 
                //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0], z->w_list.a[w_idx[0]].cidx, z->w_list.a[w_idx[0]].clen, z->w_list.c.a[z->w_list.a[w_idx[0]].cidx]);
                // }
                if(rref) {
                    calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
                } else {
                    calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
                }
                // if(j == 248) {
                //     fprintf(stderr, "4-[M::%s] j::%lld, nw::%ld, a_nw::%ld, z->x_pos_s::%u, z->x_pos_e::%u, z->y_pos_s::%u, z->y_pos_e::%u, w_idx[0]::%lu, w_idx[0]->cidx::%u, w_idx[0]->clen::%u, w_idx[0]->cigar[0]:%u\n", __func__, 
                //     j, nw, a_nw, z->x_pos_s, z->x_pos_e, z->y_pos_s, z->y_pos_e, w_idx[0], z->w_list.a[w_idx[0]].cidx, z->w_list.a[w_idx[0]].clen, z->w_list.c.a[z->w_list.a[w_idx[0]].cidx]);
                // }
                // if((int64_t)z->x_pos_s!=z->w_list.a[0].x_start || 
                //         (int64_t)z->x_pos_e!=z->w_list.a[z->w_list.n-1].x_end) {
                //     fprintf(stderr, "[M::%s] z::x_pos_s->%u, z::x_pos_e->%u, (0)::x_start->%d, (wn-1)x_end->%d, z->w_list.n->%ld\n", __func__,
                //      z->x_pos_s, z->x_pos_e, z->w_list.a[0].x_start, z->w_list.a[z->w_list.n-1].x_end, (int64_t)z->w_list.n);
                // }

                // assert(get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s)==(int64_t)z->w_list.n);
                // assert((int64_t)z->x_pos_s==z->w_list.a[0].x_start && 
                //                             (int64_t)z->x_pos_e==z->w_list.a[z->w_list.n-1].x_end);
            } else if (error_rate <= /**asm_opt.max_ov_diff_final**/e_rate_final * 1.5) {
                z->is_match = 3;
            }
            // fprintf(stderr, "[M::%s::idx->%lld::is_match->%u] z::y_id->%u, z::x_pos_s->%u, z::x_pos_e->%u, error_rate->%f, e_threshold->%f\n", 
            // __func__, j, z->is_match, z->y_id,  z->x_pos_s, z->x_pos_e, error_rate, e_rate);
        } else {///it impossible to be matched
            z->is_match = 0;
            // fprintf(stderr, "[M::%s::idx->%ld::is_match->%u] z::x_pos_s->%u, z::x_pos_e->%u, error_rate->-1, e_threshold->%f\n", 
            // __func__, j, z->is_match, z->x_pos_s, z->x_pos_e, e_rate);
        }

        
    }
    ///debug_window_cigar(overlap_list, g_read, dumy, rref, 1, 1);
}

uint32_t inline simi_pass(int64_t ol, int64_t aln_ol, uint32_t second_ck, double o_rate, double *e_rate)
{
    if(aln_ol == 0 || ol == 0) return 0;
    if((!second_ck) && (!e_rate)) {
        // if((ol*OVERLAP_THRESHOLD_FILTER) <= aln_ol) return 1;
        if((ol*o_rate) <= aln_ol) return 1;
    } else if(e_rate) {
        if((ol*((double)(((double)1.0)-(*e_rate)))) <= aln_ol) return 1;
    } else if(second_ck) {
        if(((ol*MIN_UL_ALIN_RATE) <= aln_ol) && (aln_ol >= MIN_UL_ALIN_LEN)) return 1;
    }

    // if(rref) {
    //     if((ol*OVERLAP_THRESHOLD_FILTER) <= aln_ol) return 1;
    // } else if(uref) {
    //     if(e_rate) {
    //         if((ol*((double)(((double)1.0)-(*e_rate)))) <= aln_ol) return 1;
    //     } else {
    //         if(((ol*MIN_UL_ALIN_RATE) <= aln_ol) && (aln_ol >= MIN_UL_ALIN_LEN)) return 1;
    //     }
    // }

    return 0;
}

inline uint32_t gen_backtrace(window_list *p, overlap_region *z, All_reads *rref, const ul_idx_t *uref, UC_Read* g_read, Correct_dumy* dumy,
int32_t y_strand, int32_t y_id)
{
    int64_t x_start, x_end, x_len, Window_Len, o_len; 
    int32_t threshold; long long y_start; 
    int real_y_start = 0, end_site, extra_begin, extra_end;
    char *x_string, *y_string; unsigned int error;
    ///there is no problem for x
    x_start = p->x_start; x_end = p->x_end; x_len = x_end - x_start + 1;
    threshold = p->error_threshold; Window_Len = x_len + (threshold << 1);


    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    y_start = p->y_start; extra_begin = p->extra_begin; extra_end = p->extra_end;
    o_len = Window_Len - extra_end - extra_begin;
    if(rref) {
        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
    } else {
        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
    }
    x_string = g_read->seq + x_start; y_string = dumy->overlap_region;


    ///note!!! need notification
    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
    &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - y_start);
    // assert(error != (unsigned int)-1);
    if(error != (unsigned int)-1) {
        ///this condition is always wrong
        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
        if (end_site == Window_Len - 1 || real_y_start == 0) {
            if(rref) {
                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                extra_begin, extra_end, y_id, Window_Len, rref, dumy, y_strand, error,
                &y_start, &real_y_start, &end_site, &extra_begin, &extra_end, &error)) {
                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                }
            } else {
                if(fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                &extra_end, &error)) {
                    p->error = error; p->extra_begin = extra_begin; p->extra_end = extra_end;
                }
            }
        }

        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);    

        ///note!!! need notification
        real_y_start = y_start + real_y_start - extra_begin;
        p->y_start = real_y_start;  
        p->y_end = y_start + end_site - extra_begin;
        p->error = error;         
        return 1;                     
    }
    p->error = -1;
    return 0;
}


inline uint32_t gen_backtrace_adv(window_list *p, overlap_region *z, All_reads *rref, hpc_t *hpc_g, const ul_idx_t *uref, 
char *qstr, char *tstr, char *tstr1, Correct_dumy* dumy, uint32_t rev, uint32_t id)
{
    int64_t qs, qe, ql, aln_l, t_pri_l, thres, ts; 
    int r_ts = 0, t_end; int64_t aux_beg, aux_end;
    char *q_string, *t_string; unsigned int error;
    ///there is no problem for x
    qs = p->x_start; qe = p->x_end; ql = qe + 1 - qs;
    thres = p->error_threshold; aln_l = ql + (thres<<1);

    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    ts = p->y_start; aux_beg = p->extra_begin; aux_end = p->extra_end;
    t_pri_l = aln_l - aux_beg - aux_end;

    q_string = qstr + qs; 
    if(rref) {
        fill_subregion(tstr, ts, t_pri_l, rev, rref, id, aux_beg, aux_end); t_string = tstr;
    } else {
        t_string = return_str_seq(tstr, ts, t_pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
    }
    
    t_end = Reserve_Banded_BPM_PATH(t_string, aln_l, q_string, ql, thres, &error, &r_ts,
        &(dumy->path_length), dumy->matrix_bit, dumy->path, p->error, p->y_end - ts);

    // assert(error != (unsigned int)-1);
    if(error != (unsigned int)-1) {
        /**
        bit_extz_t exz, exz64; init_bit_extz_t(&exz, thres); init_bit_extz_t(&exz64, thres);
        // ed_band_cal_extension_64_0_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // ed_band_cal_extension_64_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_extension_64_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));

        // ed_band_cal_extension_infi_0_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // ed_band_cal_extension_infi_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_infi_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);
        
        // ed_band_cal_extension_256_0_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // ed_band_cal_extension_256_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_256_0_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // char *qr, *tr;
        // gen_rev_str(q_string, &qr, ql); gen_rev_str(t_string+r_ts, &tr, t_end+1-r_ts); 

        // ed_band_cal_extension_64_1_w(tr, t_end+1-r_ts, qr, ql, thres, &exz64);
        // ed_band_cal_extension_64_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz64);
        // cigar_check(tr, qr, &exz64);
        // exz64.err = INT32_MAX;
        // ed_band_cal_extension_64_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz64);
        // cigar_check(tr, qr, &exz64);
        // assert(exz.err == exz64.err && exz.ps == (exz64.pl-exz64.pe-1) && exz.pe == (exz64.pl-exz64.ps-1) && exz.ts == (exz64.tl-exz64.te-1) && exz.te == (exz64.tl-exz64.ts-1));

        // ed_band_cal_extension_infi_1_w(tr, t_end+1-r_ts, qr, ql, thres, NULL, &exz);
        // ed_band_cal_extension_infi_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, NULL, &exz);
        // cigar_check(tr, qr, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_infi_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, NULL, &exz);
        // cigar_check(tr, qr, &exz);
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // ed_band_cal_extension_256_1_w(tr, t_end+1-r_ts, qr, ql, thres, &exz);
        // ed_band_cal_extension_256_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz);
        // cigar_check(tr, qr, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_extension_256_1_w_trace(tr, t_end+1-r_ts, qr, ql, thres, &exz);
        // cigar_check(tr, qr, &exz);
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // free(qr); free(tr);


        // ed_band_cal_global_64_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // ed_band_cal_global_64_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_global_64_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz64);
        // cigar_check(t_string+r_ts, q_string, &(exz64));

        // ed_band_cal_global_infi_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // ed_band_cal_global_infi_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_global_infi_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // ed_band_cal_global_256_w(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // ed_band_cal_global_256_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_global_256_w_trace(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres, &exz);
        // cigar_check(t_string+r_ts, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        
        // ed_band_cal_semi_64_w(t_string, aln_l, q_string, ql, thres, &exz64);
        // ed_band_cal_semi_64_w_trace(t_string, aln_l, q_string, ql, thres, &exz64);
        // cigar_check(t_string, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_semi_64_w_trace(t_string, aln_l, q_string, ql, thres, &exz64);
        // cigar_check(t_string, q_string, &(exz64));

        // ed_band_cal_semi_infi_w(t_string, aln_l, q_string, ql, thres, NULL, &exz);
        // ed_band_cal_semi_infi_w_trace(t_string, aln_l, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_semi_infi_w_trace(t_string, aln_l, q_string, ql, thres, NULL, &exz);
        // cigar_check(t_string, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);

        // ed_band_cal_semi_256_w(t_string, aln_l, q_string, ql, thres, &exz);
        // ed_band_cal_semi_256_w_trace(t_string, aln_l, q_string, ql, thres, &exz);
        // cigar_check(t_string, q_string, &exz);
        // exz.err = INT32_MAX;
        // ed_band_cal_semi_256_w_trace(t_string, aln_l, q_string, ql, thres, &exz);
        // cigar_check(t_string, q_string, &exz);
        // assert(exz.err <= (int32_t)error && exz.err >= 0); 
        // assert(exz.err == exz64.err && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);



        ed_band_cal_semi_64_w_absent_diag(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz64);
        // ed_band_cal_semi_64_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
        //                                                                         thres, aux_beg, &exz64);
        // cigar_check(t_string+aux_beg, q_string, &(exz64));
        // exz64.err = INT32_MAX;
        // ed_band_cal_semi_64_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
        //                                                                         thres, aux_beg, &exz64);
        // cigar_check(t_string+aux_beg, q_string, &(exz64));



        ed_band_cal_semi_infi_w_absent_diag(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, NULL, &exz);
        ed_band_cal_semi_infi_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, NULL, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        exz.err = INT32_MAX;
        ed_band_cal_semi_infi_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, NULL, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        assert(exz.err <= (int32_t)error && exz.err >= 0); 
        assert(exz.err == exz64.err && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);
        exz64.err = exz.err; exz64.ps = exz.ps; exz64.pe = exz.pe; exz64.ts = exz.ts; exz64.te = exz.te;


        ed_band_cal_semi_256_w_absent_diag(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz);
        ed_band_cal_semi_256_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        exz.err = INT32_MAX;
        ed_band_cal_semi_256_w_absent_diag_trace(t_string+aux_beg, aln_l-aux_beg-aux_end, q_string, ql, 
                                                                                thres, aux_beg, &exz);
        cigar_check(t_string+aux_beg, q_string, &(exz));
        assert(exz.err <= (int32_t)error && exz.err >= 0); 
        assert(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te);
        // if((!(exz.err == exz64.err && exz.ps == exz64.ps && exz.pe == exz64.pe && exz.ts == exz64.ts && exz.te == exz64.te)) || (exz.err != (int32_t)error)) {
        //     fprintf(stderr, "\n[M::%s::semi] error::%u, ql::%ld, thres::%ld, exz.err::%d, exz64.err::%d, exz.ps::%d, exz64.ps::%d, exz.pe::%d, exz64.pe::%d, exz.ts::%d, exz64.ts::%d, exz.te::%d, exz64.te::%d\n", __func__, 
        //     error, ql, thres, exz.err, exz64.err, exz.ps, exz64.ps, exz.pe, exz64.pe, exz.ts, exz64.ts, exz.te, exz64.te);
        //     fprintf(stderr, "[tstr] %.*s\n", (int32_t)aln_l, t_string);
        //     fprintf(stderr, "[qstr] %.*s\n", (int32_t)ql, q_string);
        // }
        destroy_bit_extz_t(&exz); destroy_bit_extz_t(&exz64);

        // if(exz.err > (int32_t)error && ql == 1) {
        //     fprintf(stderr, "[M::%s::] error::%u, ed_extension::%d, ql::%ld, thres::%ld\n", 
        //     __func__, error, exz.err, ql, thres);
        //     fprintf(stderr, "[tstr] %.*s\n", t_end+1-r_ts, t_string+r_ts);
        //     fprintf(stderr, "[qstr] %.*s\n", (int32_t)ql, q_string);
        // }
        
        // assert(ed_band_cal_global(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres) == 
        //                 ed_band_cal_global_128bit(t_string+r_ts, t_end+1-r_ts, q_string, ql, thres));
        **/
        
        ///this condition is always wrong
        ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
        if (((t_end+1) == aln_l) || (r_ts == 0)) {
            if(recal_boundary(q_string, tstr1, ql, thres, ts, r_ts, t_end,
            aux_beg, aux_end, error, id, aln_l, rev, dumy, rref, hpc_g, uref, 
            &ts, &r_ts, &t_end, &aux_beg, &aux_end, &error)) {
                p->error = error; p->extra_begin = aux_beg; p->extra_end = aux_end;
                t_string = update_des_str(tstr, ts, aln_l-aux_beg-aux_end, rev, rref, hpc_g, uref,
                                                    id, aux_beg, aux_end, hpc_g?NULL:tstr1);
            }
        }

        generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &r_ts, &t_end, &error, q_string, ql, t_string);    

        p->y_start = ts + r_ts - aux_beg; 
        p->y_end = ts + t_end - aux_beg;
        p->error = error;         
        return 1;                     
    }
    p->error = -1;
    return 0;
}

inline uint32_t aln_wlst_adv(overlap_region *z, All_reads *rref, hpc_t *hpc_g, 
const ul_idx_t *uref, char *qstr, char *tstr, char *tstr1, Correct_dumy* dumy,
uint32_t rev, uint32_t id, int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, 
double e_rate, uint32_t is_cigar)
{
    int64_t ql, aln_l, t_tot_l; window_list *p = NULL; int r_ts = 0, t_end;
    int64_t aux_beg, aux_end, t_pri_l;
    int64_t thres; char *q_string, *t_string; unsigned int error;
    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD), ql);
    } else {
        thres = double_ul_error_threshold(get_init_err_thres(ql, e_rate, block_s, THRESHOLD_MAX_SIZE), ql);
    }
    aln_l = ql + (thres << 1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; 
    if(rref) {
        fill_subregion(tstr, t_s, t_pri_l, rev, rref, id, aux_beg, aux_end); t_string = tstr;
    } else {
        t_string = return_str_seq(tstr, t_s, t_pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
    }

    if(is_cigar) {
        ///note!!! need notification
        t_end = Reserve_Banded_BPM_PATH(t_string, aln_l, q_string, ql, thres, &error, &r_ts,
        &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
    } else {
        ///note!!! need notification
        t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, ql, thres, &error);    
    }
    if(error!=(unsigned int)-1) {
        if(is_cigar) {
            ///this condition is always wrong
            ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
            if (((t_end+1) == aln_l) || (r_ts == 0)) {
                if(recal_boundary(q_string, tstr1, ql, thres, t_s, r_ts, t_end,
                aux_beg, aux_end, error, id, aln_l, rev, dumy, rref, hpc_g, uref, 
                &t_s, &r_ts, &t_end, &aux_beg, &aux_end, &error)) {
                    t_string = update_des_str(tstr, t_s, aln_l-aux_beg-aux_end, rev, rref, hpc_g, uref,
                                                    id, aux_beg, aux_end, hpc_g?NULL:tstr1);
                }
            }
        }

        kv_pushp(window_list, z->w_list, &p);
        p->x_start = qs; p->x_end = qe; ///must set x_start/x_end here
        if(is_cigar) {
            generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &r_ts, &t_end, &error, q_string, ql, t_string);  
        } else {
            p->cidx = p->clen = 0;
        }
        p->y_start = t_s + r_ts;///difference
        p->y_end = t_s + t_end;
        p->error = error;
        p->extra_begin = aux_beg;
        p->extra_end = aux_end;
        p->error_threshold = thres;
        z->align_length += ql; 
        
        return 1;
    }
    return 0;
}

void push_wcigar(window_list *idx, window_list_alloc *res, bit_extz_t *exz)
{
    idx->cidx = res->c.n; idx->clen = exz->cigar.n; res->c.n += exz->cigar.n;
    kv_resize(uint16_t, res->c, res->c.n);
    memcpy(res->c.a+idx->cidx, exz->cigar.a, exz->cigar.n*sizeof(*(res->c.a)));
}

inline uint32_t aln_wlst_adv_exz(overlap_region *z, All_reads *rref, hpc_t *hpc_g, 
const ul_idx_t *uref, char *qstr, char *tstr, bit_extz_t *exz, uint32_t max_err,
uint32_t rev, uint32_t id, int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, 
double e_rate, uint32_t is_cigar)
{
    int64_t ql, tl, aln_l, t_tot_l; window_list *p = NULL; ///int r_ts = 0, t_end;
    int64_t aux_beg, aux_end, t_pri_l; int64_t thres; char *q_string, *t_string; 
    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, max_err), ql);
    } else {
        thres = double_ul_error_threshold(get_init_err_thres(ql, e_rate, block_s, max_err), ql);
    }
    aln_l = ql + (thres << 1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; 
    if(rref) {
        recover_UC_Read_sub_region(tstr, t_s, t_pri_l, rev, rref, id); t_string = tstr;
    } else {
        t_string = return_str_seq_exz(tstr, t_s, t_pri_l, rev, hpc_g, uref, id);
    }
    tl = t_pri_l;
    if(is_cigar) {
        clear_align(*exz);
        ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thres, aux_beg, exz);
    } else {
        ed_band_cal_semi_64_w_absent_diag(t_string, tl, q_string, ql, thres, aux_beg, exz); exz->ps = 0;
    }

    // if(id == 40 && qs == 79670 && qe == 79824) {
    //     fprintf(stderr, "\n[M::%s::semi] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
    //     __func__, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
    // }

    if(is_align(*exz)) {
        kv_pushp(window_list, z->w_list, &p);
        p->x_start = qs; p->x_end = qe; ///must set x_start/x_end here
        p->y_start = t_s + exz->ps;///difference
        p->y_end = t_s + exz->pe;
        p->error = exz->err;
        p->cidx = p->clen = 0;
        if(is_cigar) {
            push_wcigar(p, &(z->w_list), exz);
            ///this condition is always wrong
            ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
            if ((((exz->pe+1) == tl) || (exz->ps == 0)) && (exz->err > 0)) {
                if(recal_boundary_exz(q_string, tstr, ql, tl, thres, t_s, exz->ps, exz->pe,
                exz->err, id, rev, exz, rref, hpc_g, uref, &t_s, &aux_beg, &aux_end)) {
                    //update cigar
                    z->w_list.c.n = p->cidx; push_wcigar(p, &(z->w_list), exz);

                    p->y_start = t_s + exz->ps;///difference
                    p->y_end = t_s + exz->pe;
                    p->error = exz->err;
                }
            }
        } 

        p->extra_begin = aux_beg;
        p->extra_end = aux_end;
        p->error_threshold = thres;
        z->align_length += ql; 
        return 1;
    }
    return 0;
}

inline uint32_t aln_wlst_adv_non_retrieve_exz(overlap_region *z, char *qstr, char *tstr, int64_t t_tot_l, bit_extz_t *exz, uint32_t max_err,
uint32_t rev, uint32_t id, int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, double e_rate, uint32_t is_cigar)
{
    int64_t ql, tl, aln_l; window_list *p = NULL; ///int r_ts = 0, t_end;
    int64_t aux_beg, aux_end, t_pri_l; int64_t thres; char *q_string, *t_string; 
    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, max_err), ql);

    aln_l = ql + (thres << 1);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; t_string = tstr + t_s;

    tl = t_pri_l;
    if(is_cigar) {
        clear_align(*exz);
        ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thres, aux_beg, exz);
    } else {
        ed_band_cal_semi_64_w_absent_diag(t_string, tl, q_string, ql, thres, aux_beg, exz); exz->ps = 0;
    }

    // if(id == 40 && qs == 79670 && qe == 79824) {
    //     fprintf(stderr, "\n[M::%s::semi] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
    //     __func__, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
    // }

    if(is_align(*exz)) {
        kv_pushp(window_list, z->w_list, &p);
        p->x_start = qs; p->x_end = qe; ///must set x_start/x_end here
        p->y_start = t_s + exz->ps;///difference
        p->y_end = t_s + exz->pe;
        p->error = exz->err;
        p->cidx = p->clen = 0;
        if(is_cigar) {
            push_wcigar(p, &(z->w_list), exz);
            ///this condition is always wrong
            ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
            if ((((exz->pe+1) == tl) || (exz->ps == 0)) && (exz->err > 0)) {
                if(recal_boundary_non_retrieve_exz(q_string, tstr, t_tot_l, ql, tl, thres, t_s, exz->ps, exz->pe,
                exz->err, id, rev, exz, &t_s, &aux_beg, &aux_end)) {
                    //update cigar
                    z->w_list.c.n = p->cidx; push_wcigar(p, &(z->w_list), exz);

                    p->y_start = t_s + exz->ps;///difference
                    p->y_end = t_s + exz->pe;
                    p->error = exz->err;
                }
            }
        } 

        p->extra_begin = aux_beg;
        p->extra_end = aux_end;
        p->error_threshold = thres;
        z->align_length += ql; 
        return 1;
    }
    return 0;
}

inline uint32_t aln_wlst(overlap_region *z, All_reads *rref, const ul_idx_t *uref, UC_Read* g_read, Correct_dumy* dumy,
int32_t y_strand, int32_t y_id, int64_t x_start, int64_t x_end, long long y_start, int64_t block_s, double e_rate, int32_t is_cigar)
{
    int64_t x_len, Window_Len; window_list *p = NULL; long long o_len;
    int32_t threshold; int real_y_start = 0, end_site, extra_begin, extra_end;
    char *x_string, *y_string; unsigned int error;
    x_len = x_end + 1 - x_start;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        threshold = double_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD), x_len);
    } else {
        threshold = double_ul_error_threshold(get_init_err_thres(x_len, e_rate, block_s, THRESHOLD_MAX_SIZE), x_len);
    }
    Window_Len = x_len + (threshold << 1);
    ///y_start might be less than 0
    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len)), 
    &extra_begin, &extra_end, &y_start, &o_len)) {
        return 0;
    }

    if(o_len + threshold < x_len) return 0;
    if(rref) {
        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, rref, y_id, extra_begin, extra_end);
    } else {
        fill_subregion_ul(dumy->overlap_region, y_start, o_len, y_strand, uref, y_id, extra_begin, extra_end);
    }
    x_string = g_read->seq + x_start;
    y_string = dumy->overlap_region;

    if(is_cigar) {
        ///note!!! need notification
        end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
        &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);
    } else {
        ///note!!! need notification
        end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);    
    }
    if(error!=(unsigned int)-1) {
        if(is_cigar) {
            ///this condition is always wrong
            ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
            if (end_site == Window_Len - 1 || real_y_start == 0) {
                if(rref) {
                    fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                    end_site, extra_begin, extra_end, y_id, Window_Len, rref, dumy, 
                    y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                    &extra_end, &error);
                } else {
                    fix_ul_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                    end_site, extra_begin, extra_end, y_id, Window_Len, uref, dumy, 
                    y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                    &extra_end, &error);
                }
            }
        }

        kv_pushp(window_list, z->w_list, &p);
        p->x_start = x_start; p->x_end = x_end; ///must set x_start/x_end here
        if(is_cigar) {
            generate_cigar(dumy->path, dumy->path_length, p, &(z->w_list), &real_y_start, &end_site, &error, x_string, x_len, y_string);  
        } else {
            p->cidx = p->clen = 0;
        }
        p->y_start = y_start + real_y_start;///difference
        p->y_end = y_start + end_site;
        p->error = error;
        p->extra_begin = extra_begin;
        p->extra_end = extra_end;
        p->error_threshold = threshold;
        z->align_length += x_len; 
        
        return 1;
    }
    return 0;
}
uint64_t realign_ed(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, 
double *e_rate_final, uint32_t sec_check, int64_t *is_sort);
inline void refine_ed_aln(overlap_region_alloc* overlap_list, All_reads *rref, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, double e_rate_final)
{
    int64_t j, k, i, on, y_id, y_readLen, x_start, x_end, x_len, total_y_start, total_y_end;
    int32_t y_strand, real_y_start;
    int64_t nw, a_nw, w_id, w_s, w_e, is_srt, mm_we, mm_ws, mm_aln, ovl;
    double error_rate; uint64_t *w_idx; overlap_region *z; window_list *p = NULL;

    overlap_list->mapped_overlaps_length = 0; on = overlap_list->length;
    for (j = 0; j < on; ++j) {
        // z = &(overlap_list->list[j]); ovl = z->x_pos_e+1-z->x_pos_s;
        // if(!realign_ed(z, uref, NULL, rref, g_read->seq, 
        //         dumy->overlap_region, dumy->overlap_region_fix, dumy, v_idx, block_s, e_rate, NULL, 1, &is_srt)) {
        //     continue;
        // }
        z = &(overlap_list->list[j]); z->is_match = 0; is_srt = 1;
        if(z->w_list.n == 0) continue;///no alignment
        nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s); a_nw = z->w_list.n;
        kv_resize(uint64_t, v_idx->a, (uint64_t)nw); 
        memset(v_idx->a.a, -1, sizeof((*v_idx->a.a))*nw); w_idx = v_idx->a.a;
        for (i = 0; i < a_nw; i++) { ///w_idx[] == (uint64_t) if unmatched
            assert(z->w_list.a[i].y_end != -1);
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, NULL);
            w_idx[w_id] = i;
        }

        y_id = z->y_id; y_strand = z->y_pos_strand; 
        ovl = z->x_pos_e+1-z->x_pos_s; mm_we = z->x_pos_s; mm_aln = 0; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        for (i = a_nw-1; i >= 0; i--) { //utilize the the end pos of pre-window in forward
            w_id = get_win_id_by_s(z, z->w_list.a[i].x_start, block_s, &w_e);
            assert(z->w_list.a[i].x_end == w_e);
            if(w_e > mm_we) mm_we = w_e;
            ///in most cases, extra_begin = 0
            total_y_start = z->w_list.a[i].y_end + 1 - z->w_list.a[i].extra_begin;
            for (k = w_id + 1; k < nw && total_y_start < y_readLen; k++) {
                if(w_idx[k] != (uint64_t)-1) break;
                w_s = w_e + 1;
                w_id = get_win_id_by_s(z, w_s, block_s, &w_e);
                assert(w_id == k);
                x_start = w_s; x_end = w_e;
                if(aln_wlst(z, rref, uref, g_read, dumy, y_strand, y_id, x_start, x_end, 
                                                            total_y_start, block_s, e_rate, 0)) {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    w_idx[k] = z->w_list.n - 1;
                    if(x_end > mm_we) mm_we = x_end;
                    if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                } else {
                    break;
                }
                total_y_start = p->y_end + 1 - p->extra_begin;
            }
        }
        mm_ws = z->x_pos_s; mm_aln = mm_we+1-mm_ws; 
        if(!simi_pass(ovl, mm_aln, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) continue;
        if(nw > 0 && w_idx[0] != (uint64_t)-1) mm_ws = z->w_list.a[w_idx[0]].x_end+1;
        
        for (i = 1; i < nw; i++) { //utilize the the start pos of next window in backward
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(w_idx[i] != (uint64_t)-1 && w_idx[i-1] == (uint64_t)-1) {
                w_s = z->w_list.a[w_idx[i]].x_start; mm_aln -= (w_s-mm_ws);
                ///check if the start pos of this matched window has been calculated
                if(z->w_list.a[w_idx[i]].clen == 0) {
                    p = &(z->w_list.a[w_idx[i]]);
                    gen_backtrace(p, z, rref, uref, g_read, dumy, y_strand, y_id);
                    assert(p->error != -1);
                    p->y_end += p->extra_begin;
                } 
                real_y_start = p->y_start;

                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && w_idx[k] == (uint64_t)-1 && total_y_end > 0; k--) {  
                    w_e = w_s - 1;
                    w_id = get_win_id_by_e(z, w_e, block_s, &w_s);
                    assert(w_id == k);
                    x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
                    if(aln_wlst(z, rref, uref, g_read, dumy, y_strand, y_id, x_start, x_end, total_y_end+1-x_len, block_s, e_rate, 1)) {
                        p = &(z->w_list.a[z->w_list.n-1]);
                        p->y_start -= p->extra_begin; ///y_start has no shift, but y_end has shift  
                        w_idx[k] = z->w_list.n - 1;
                        mm_aln += x_len;
                        if(is_srt && z->w_list.n > 1 && p->x_start < z->w_list.a[z->w_list.n-2].x_start) is_srt = 0;
                    } else {
                        break;
                    }
                    total_y_end = p->y_start - 1;
                }
                if(!simi_pass(ovl, mm_aln, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) break;
            }
            if(w_idx[i] != (uint64_t)-1) mm_ws = z->w_list.a[w_idx[i]].x_end+1;
        }

        if(i < nw) continue;
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) {
            z->is_match = 3; overlap_list->mapped_overlaps_length += z->align_length;
            ///sort for set_herror_win
            if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
        }
    }

    if(uref && overlap_list->mapped_overlaps_length > 0) {
        set_herror_win(overlap_list, dumy, v_idx, e_rate, g_read->length, block_s);
    }

    overlap_list->mapped_overlaps_length = 0;
    for (j = 0; j < (long long)overlap_list->length; j++) {
        z = &(overlap_list->list[j]);
        y_id = z->y_id; y_strand = z->y_pos_strand; 
        y_readLen = (rref?(Get_READ_LENGTH((*rref), y_id)):(uref->ug->u.a[y_id].len));
        ovl = z->x_pos_e + 1 - z->x_pos_s; //z->is_match = 0;
        // if(y_id == 4) {
        //     fprintf(stderr, "[M::%s::idx->%ld::] z::x_pos_s->%u, z::x_pos_e->%u, ovl->%ld, aln->%u\n", 
        //     __func__, j, z->x_pos_s, z->x_pos_e, ovl, z->align_length);
        // }
        ///debug_scan_cigar(&(overlap_list->list[j]));
        ///only calculate cigar for high quality overlaps
        // int64_t tt = 0;
        if(simi_pass(ovl, z->align_length, 0, OVERLAP_THRESHOLD_NOSI_FILTER, &e_rate)) {
            a_nw = z->w_list.n;
            for (i = 0, is_srt = 1; i < a_nw; i++) {
                p = &(z->w_list.a[i]);
                ///check if the cigar of this window has been got 
                if(p->clen == 0) {
                    gen_backtrace(p, z, rref, uref, g_read, dumy, y_strand, y_id);
                    assert(p->error != -1);
                    // if(y_id == 4) {
                    //     fprintf(stderr, "+[M::idx->%ld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                else {
                    p->y_end -= p->extra_begin;
                    // if(y_id == 4) {
                    //     fprintf(stderr, "-[M::idx->%ld::] y_start->%d, y_end->%d, error->%d\n", 
                    //     j, p->y_start, p->y_end, p->error);
                    // }
                }
                // tt += p->error;
                if(is_srt && i > 0 && p->x_start < z->w_list.a[i-1].x_start) is_srt = 0;
            }
            if(!is_srt) radix_sort_window_list_xs_srt(z->w_list.a, z->w_list.a + z->w_list.n);
            error_rate = non_trim_error_rate(z, rref, uref, v_idx, dumy, g_read, e_rate, block_s);
            z->is_match = 0;///must be here;
            // if(y_id == 4) {
            //     fprintf(stderr, "[M::%s::idx->%ld::] block_s->%ld, z::x_pos_s->%u, z::x_pos_e->%u, ovl->%ld, aln->%u, error_rate->%f, e_rate_final->%f\n", 
            //     __func__, j, block_s, z->x_pos_s, z->x_pos_e, ovl, z->align_length, error_rate, e_rate_final);
            //     exit(1);
            // }
            if (error_rate <= e_rate_final) {
                overlap_list->mapped_overlaps_length += ovl;
                z->is_match = 1; append_unmatched_wins(z, block_s);
                if(rref) {
                    calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
                } else {
                    calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
                }
                // assert(get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s)==(int64_t)z->w_list.n);
                // assert((int64_t)z->x_pos_s==z->w_list.a[0].x_start && 
                //                             (int64_t)z->x_pos_e==z->w_list.a[z->w_list.n-1].x_end);
            } else if (error_rate <= e_rate_final * 1.5) {
                z->is_match = 3;
            }
        } else {///it impossible to be matched
            z->is_match = 0;
            // fprintf(stderr, "[M::%s::idx->%ld::is_match->%u] z::x_pos_s->%u, z::x_pos_e->%u, error_rate->-1, e_threshold->%f\n", 
            // __func__, j, z->is_match, z->x_pos_s, z->x_pos_e, e_rate);
        }
    }
}


uint32_t align_ul_ed_post(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, char *tstr_1, 
Correct_dumy* dumy, double e_rate, int64_t w_l, double ovlp_cut, void *km);
double gen_extend_err(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, char *tstr_1, Correct_dumy* dumy, uint64_t *v_idx, int64_t block_s, double ovlp_cut, double e_rate, double e_max, int64_t *r_e);
inline void refine_ed_aln_test(overlap_region_alloc* overlap_list, All_reads *rref, const ul_idx_t *uref, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read, kvec_t_u64_warp* v_idx, int64_t block_s, double e_rate, double e_rate_final)
{
    int64_t j, on, ovl; uint64_t k; double rr; overlap_region *z; 

    overlap_list->mapped_overlaps_length = 0; on = overlap_list->length;
    for (j = 0; j < on; ++j) {
        z = &(overlap_list->list[j]); ovl = z->x_pos_e+1-z->x_pos_s;
        if(!align_ul_ed_post(z, uref, NULL, g_read->seq, dumy->overlap_region, dumy->overlap_region_fix, 
                        dumy, e_rate, block_s, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) {
            continue;
        }
        if(uref && simi_pass(ovl, z->align_length, uref?1:0, OVERLAP_THRESHOLD_NOSI_FILTER, NULL)) {
            z->is_match = 3; overlap_list->mapped_overlaps_length += z->align_length;
        }
    }

    if(uref && overlap_list->mapped_overlaps_length > 0) {
        set_herror_win(overlap_list, dumy, v_idx, e_rate, g_read->length, block_s);
    }

    double e_max = e_rate_final * 1.5;
    overlap_list->mapped_overlaps_length = 0; on = overlap_list->length;
    for (j = 0; j < on; j++) {
        z = &(overlap_list->list[j]); ovl = z->x_pos_e + 1 - z->x_pos_s;
        rr = gen_extend_err(z, uref, NULL, rref, g_read->seq, dumy->overlap_region, dumy->overlap_region_fix,
                    dumy, v_idx?v_idx->a.a:NULL, block_s, -1, e_rate, (e_max+0.000001), NULL);
        z->is_match = 0;///must be here;
        if (rr <= e_rate_final) {
            for (k = 0; k < z->w_list.n; k++) {
                if(z->w_list.a[k].clen) continue;
                gen_backtrace_adv(&(z->w_list.a[k]), z, rref, NULL, uref, g_read->seq, dumy->overlap_region, dumy->overlap_region_fix, 
                dumy, z->y_pos_strand, z->y_id);
            }
            
            overlap_list->mapped_overlaps_length += ovl;
            z->is_match = 1; append_unmatched_wins(z, block_s);
            if(rref) {
                calculate_boundary_cigars(z, rref, dumy, g_read, e_rate);
            } else {
                calculate_ul_boundary_cigars(z, uref, dumy, g_read, e_rate, block_s);
            }
        } else if (rr <= e_max) {
            z->is_match = 3;
        }
    }
}



inline void add_base_to_correct_read_directly(Correct_dumy* dumy, char base)
{
  
    if (dumy->corrected_read_length + 2 > dumy->corrected_read_size)
    {
        dumy->corrected_read_size = dumy->corrected_read_size * 2;
        dumy->corrected_read = (char*)realloc(dumy->corrected_read, dumy->corrected_read_size);
    }
    
    dumy->corrected_read[dumy->corrected_read_length] = base;
    dumy->corrected_read_length++;
    dumy->corrected_read[dumy->corrected_read_length] = '\0';

}

inline void add_base_to_correct_read(Correct_dumy* dumy, char base, int is_error)
{
    ///don't need to deal with deletion
    if (base != 'D')
    {
        if (dumy->corrected_read_length + 2 > dumy->corrected_read_size)
        {
            dumy->corrected_read_size = dumy->corrected_read_size * 2;
            dumy->corrected_read = (char*)realloc(dumy->corrected_read, dumy->corrected_read_size);
        }
        
        dumy->corrected_read[dumy->corrected_read_length] = base;
        dumy->corrected_read_length++;
        dumy->corrected_read[dumy->corrected_read_length] = '\0';
    }

    if (is_error)
    {
        dumy->corrected_base++;
    }
    
}


inline void add_segment_to_correct_read(Correct_dumy* dumy, char* segment, long long segment_length)
{

    if (dumy->corrected_read_length + segment_length + 2 > dumy->corrected_read_size)
    {
        dumy->corrected_read_size = dumy->corrected_read_length + segment_length + 2;
        dumy->corrected_read = (char*)realloc(dumy->corrected_read, dumy->corrected_read_size);
    }

    memcpy(dumy->corrected_read + dumy->corrected_read_length, segment, segment_length);
    dumy->corrected_read_length += segment_length;
    dumy->corrected_read[dumy->corrected_read_length] = '\0';
}




///return the ID of next node at backbone


///return the ID of next node at backbone









inline void generate_seq_from_path(Graph* DAGCon, Node* node, int direction)
{
    clear_Queue(&(DAGCon->node_q));
    RSet iter;
    Edge* e = NULL;
    uint64_t max;
    Node* max_node = NULL;
    

    if(direction == 0)
    {
        while (node->ID != DAGCon->s_end_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getOutputEdges(&iter, DAGCon, node, &e))
            {
                if(e->weight > max)
                {
                    max = e->weight;
                    max_node = &(G_Node(*DAGCon, e->out_node));
                }
            }
            node = max_node;
        }
    }
    else
    {
        while (node->ID != DAGCon->s_start_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getInputEdges(&iter, DAGCon, node, &e))
            {
                if(e->weight > max)
                {
                    max = e->weight;
                    max_node = &(G_Node(*DAGCon, e->in_node));
                }
            }
            node = max_node;
        }

        long long i, k;
        long long length = (DAGCon->node_q.end - DAGCon->node_q.beg);
        long long length_ex = length/2;
        long long* array = DAGCon->node_q.buffer + DAGCon->node_q.beg;
        for (i = 0; i < length_ex; i++)
        {
            k = array[i];
            array[i] = array[length - i - 1];
            array[length - i - 1] = k;
        }
    }
}




inline void generate_seq_from_node(Graph* DAGCon, Node* node, int direction)
{
    clear_Queue(&(DAGCon->node_q));
    RSet iter;
    uint64_t max;
    Node* max_node = NULL;
    Node* getNodes = NULL;
    

    if(direction == 0)
    {
        while (node->ID != DAGCon->s_end_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getOutputNodes(&iter, DAGCon, node, &getNodes))
            {
                if(getNodes->weight > max)
                {
                    max = getNodes->weight;
                    max_node = getNodes;
                }
            }
            node = max_node;
        }
    }
    else
    {
        while (node->ID != DAGCon->s_start_nodeID)
        {
            push_to_Queue(&(DAGCon->node_q), node->base);
            clear_RSet(&iter);
            max = 0;
            while(getInputNodes(&iter, DAGCon, node, &getNodes))
            {
                if(getNodes->weight > max)
                {
                    max = getNodes->weight;
                    max_node = getNodes;
                }
            }
            node = max_node;
        }

        long long i, k;
        long long length = (DAGCon->node_q.end - DAGCon->node_q.beg);
        long long length_ex = length/2;
        long long* array = DAGCon->node_q.buffer + DAGCon->node_q.beg;
        for (i = 0; i < length_ex; i++)
        {
            k = array[i];
            array[i] = array[length - i - 1];
            array[length - i - 1] = k;
        }
    }
}











///correct bases of current_dumy->corrected_read in [start_base, end_base]





inline int get_available_fully_covered_interval(long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, long long* real_length, long long* real_length_100)
{
    long long i, fud = 0;
    long long Len;
    long long overlap_length;

    if(window_start == 0) dumy->start_i = 0;
    for (i = dumy->start_i; i < (long long)overlap_list->length; i++)
    {
        if (window_end < (long long)overlap_list->list[i].x_pos_s)
        {
            dumy->start_i = 0;
            return 0;
        }
        else 
        {
            dumy->start_i = i;
            break;
        }
    }


    if (i >= (long long)overlap_list->length)
    {
        dumy->start_i = overlap_list->length;
        return -2;
    }
    



    
    long long fake_length = 0;
    overlap_length = window_end - window_start + 1;
    (*real_length) = 0; fud = 0;

    for (; i < (long long)overlap_list->length; i++)
    {
        if((Len = OVERLAP(window_start, window_end, (long long)overlap_list->list[i].x_pos_s, (long long)overlap_list->list[i].x_pos_e)) > 0)
        {
            fake_length++;

            if (overlap_length == Len && overlap_list->list[i].is_match == 1)
            {
                (*real_length)++;
            }

            if (overlap_length == Len && overlap_list->list[i].is_match == 100)
            {
                (*real_length_100)++;
            }
            if(fud == 0) fud = 1, dumy->start_i = i;
        }
        
        if((long long)overlap_list->list[i].x_pos_s > window_end)
        {
            break;
        }
    }

    if (fake_length == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


///mark SNPs at [xBeg, xEnd], note we need to deal with flag_offset carefully


///window_offset is still the x-based offset
///x_total_start and y_total_start are global positions, instead of local positions



/**
void addSNPtohaplotype(
long long window_offset, int overlapID,
char* x_string, long long x_total_start, long long x_length, 
char* y_string, long long y_total_start, long long y_length, 
CIGAR* cigar, haplotype_evdience_alloc* hap, int snp_threshold)
{
    
    int x_i, y_i, cigar_i;
    x_i = 0;
    y_i = 0;
    cigar_i = 0;
    int operation;
    int operationLen;
    int i;
    long long inner_offset = x_total_start - window_offset;
    haplotype_evdience ev;
    
    ///note that node 0 is the start node
    ///0 is match, 1 is mismatch, 2 is up, 3 is left
    ///2 represents thre are more bases at y
    ///3 represents thre are more bases at x
    while (cigar_i < cigar->length)
    {
        operation = cigar->C_C[cigar_i];
        operationLen = cigar->C_L[cigar_i];

        ///matches
        if (operation == 0)
        {
            for (i = 0; i < operationLen; i++)
            {
                ///should be at least 2 mismatches
                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 0;
                    addHaplotypeEvdience(hap, &ev, NULL);
                }


                inner_offset++;
                x_i++;
                y_i++;
            }

        }
        else if(operation == 1)
        {
            for (i = 0; i < operationLen; i++)
            {

                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  y_string[y_i];
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 1;
                    addHaplotypeEvdience(hap, &ev, NULL);
                }

                inner_offset++;
                x_i++;
                y_i++;
            }
        }///insertion
        else if (operation == 2)
        {
            y_i += operationLen;
        }
        else if (operation == 3)
        {
            //may have bugs
            for (i = 0; i < operationLen; i++)
            {
                if(hap->flag[inner_offset] > snp_threshold)
                {
                    ev.misBase =  'N';
                    ev.overlapID = overlapID;
                    ev.site = x_total_start + x_i;
                    ev.overlapSite = y_total_start + y_i;
                    ev.type = 2;
                    addHaplotypeEvdience(hap, &ev, NULL);
                }

                inner_offset++;
                x_i++;
            }
            //may have bugs
        }
        
        cigar_i++;
    }
}
**/


///mark SNPs at [xBeg, xEnd], note we need to deal with flag_offset carefully



/**
void cluster(char* r_string, long long window_start, long long window_end, 
overlap_region_alloc* overlap_list, Correct_dumy* dumy, All_reads* R_INF, haplotype_evdience_alloc* hap)
{
    ///window_start, window_end, and useful_length correspond to x, instead of y
    long long useful_length = window_end - window_start + 1;
    long long x_start;
    long long x_length; 
    char* x_string;
    char* y_string;
    long long i;
    long long y_start, y_length;
    long long overlapID, windowID;

    long long correct_x_pos_s;
    int snp_threshold;
    snp_threshold = 1;

    ///all overlaps related to the current window [window_start, window_end]
    ///first mark all snp pos
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap id, instead of the window id or the y id
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list[windowID].x_end 
                - overlap_list->list[overlapID].w_list[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list[windowID].y_end
                - overlap_list->list[overlapID].w_list[windowID].y_start + 1;

            
        markSNP(window_start, x_start, x_length, y_start, y_length, 
        &(overlap_list->list[overlapID].w_list[windowID].cigar), hap);
    }


    //may have bugs
    long long last_snp = -1;
    long long first_snp = -1;
    for (i = 0; i < useful_length; i++)
    {
        if(hap->flag[i] != 0)
        {
            last_snp = i;
            if(first_snp == -1)
            {
                first_snp = i;
            }
        }
        ///for a real snp, the coverage should be at least 2
        if(hap->flag[i] > snp_threshold)
        {
            // hap->snp++;
            hap->nn_snp++;
        }
    }
    ///if there are any >0 elements, both first_snp and last_snp should be != -1
    if(first_snp == -1 || last_snp == -1)
    {
        first_snp = 0;
        last_snp = -1;
    }
    //may have bugs



    ///add the information related to snp to haplotype_evdience_alloc
    for (i = 0; i < (long long)dumy->length; i++)
    {
        ///overlap ID, instead of the window ID
        overlapID = dumy->overlapID[i];

        ///overlap_list->list[overlapID].x_pos_s is the begining of the whole overlap
        correct_x_pos_s = (overlap_list->list[overlapID].x_pos_s / WINDOW) * WINDOW;
        ///window_start is the begining of this window in the whole x_read
        windowID = (window_start - correct_x_pos_s) / WINDOW;

        ///skip if this window is not matched
        if (overlap_list->list[overlapID].w_list[windowID].y_end == -1)
        {
            continue;
        }
        
        ///both x_start and y_start are the offsets of the whole x_read and y_read
        ///instead of the offsets of window
        x_start = overlap_list->list[overlapID].w_list[windowID].x_start;
        x_length = overlap_list->list[overlapID].w_list[windowID].x_end 
                - overlap_list->list[overlapID].w_list[windowID].x_start + 1;

        y_start = overlap_list->list[overlapID].w_list[windowID].y_start;
        y_length = overlap_list->list[overlapID].w_list[windowID].y_end
                - overlap_list->list[overlapID].w_list[windowID].y_start + 1;


        recover_UC_Read_sub_region(dumy->overlap_region, y_start, y_length, overlap_list->list[overlapID].y_pos_strand, 
                R_INF, overlap_list->list[overlapID].y_id);

        x_string = r_string + x_start;
        y_string = dumy->overlap_region;


        addSNPtohaplotype(window_start, overlapID, x_string, x_start, x_length, 
        y_string, y_start, y_length, &(overlap_list->list[overlapID].w_list[windowID].cigar), 
        hap, snp_threshold);        
    }

    RsetInitHaplotypeEvdienceFlag(hap, first_snp, last_snp + 1 - first_snp);
}
**/













































inline int check_informative_site(haplotype_evdience_alloc* hap, SnpStats* snp)
{
    long long vectorID = snp->id;
    int8_t *vector = Get_SNP_Vector((*hap), vectorID);
    snp->occ_0 = 0;
    snp->occ_1 = 0;
    snp->occ_2 = 0;
    long long i;
    for (i = 0; i < Get_SNP_Vector_Length((*hap)); i++)
    {
        if(vector[i] == 0)
        {
            snp->occ_0++;
        }
        else if(vector[i] == 1)
        {
            snp->occ_1++;
        }
        else if(vector[i] == 2)
        {
            snp->occ_2++;
        }
    }

    if(snp->occ_0 >= 2 || snp->occ_1 >= 2)
    {
        return 1;
    }

    return 0;
}







#define is_st_bs(s, rr, mm) (((mm) != ((uint64_t)-1)) && (((s).overlap_num + mm) >= ((s).occ_0)) && ((((s).occ_0*(rr) + (s).overlap_num)) >= ((s).occ_0)))





inline int64_t comput_sc_rphase(SnpStats *ai, uint64_t id, SnpStats *aj, uint64_t jd, haplotype_evdience *za, uint64_t occ0_cut)
{
    if(ai->site == aj->site) return INT64_MIN;
    // if(ai->occ_0 < occ0_cut || aj->occ_0 < occ0_cut) return INT64_MIN;
    haplotype_evdience *iz = NULL, *jz = NULL; int64_t in, jn, ik, jk, nn[2]; uint8_t fi, fj;
    iz = za + ai->non_homopolymer_num; in = ai->homopolymer_num - ai->non_homopolymer_num;
    jz = za + aj->non_homopolymer_num; jn = aj->homopolymer_num - aj->non_homopolymer_num;

    for (ik = jk = nn[0] = nn[1] = 0; (ik < in) && (jk < jn); ik++) {
        for (; (jk < jn) && (jz[jk].overlapID < iz[ik].overlapID); jk++);
        if((jk < jn) && (jz[jk].overlapID == iz[ik].overlapID)) {

            fi = 2;
            if(hh_tp(iz[ik]) == 0) {
                fi = 0;
            } else if(iz[ik].overlapSite == id){
                fi = 1;
            }

            fj = 2;
            if(hh_tp(jz[jk]) == 0) {
                fj = 0;
            } else if(jz[jk].overlapSite == jd){
                fj = 1;
            }

            if((fi == 2) && (fj == 2) && (iz[ik].overlapSite != ((uint32_t)-1)) && (jz[jk].overlapSite != ((uint32_t)-1))) {///for rare cases
                fi = fj = 0;
            }

            if(fi == 2 || fj == 2) return INT64_MIN;
            if(fi != fj) return INT64_MIN;
            nn[fi]++;
        }
    }

    if(nn[0] > 0 && nn[1] > 0) return 1;
    return INT64_MIN;
}

///idx->a:: [0, ch_n) -> tree;








inline void fill_incom(asg64_v *om, uint64_t oid, uint64_t pe, uint64_t *idx_a, int64_t idx_n, uint64_t qid)
{
    if(om->a[oid] == ((uint64_t)-1)) {
        om->a[oid] = pe;
        return;
    }

    uint64_t ps = om->a[oid]; int64_t k;
    for (k = idx_n - 1; idx_a[k] != ps; k--);
    // if(!(k >= 0)) {
    //     fprintf(stderr, "qid::%lu, oid::%lu, pe::%lu, ps::%lu, idx_n::%ld\n", qid, oid, pe, ps, idx_n);
    // }
    assert(k >= 0);
    for (k++; k < idx_n; k++) {
        kv_push(uint64_t, *om, ((oid<<32)|(idx_a[k])));
    }
    om->a[oid] = pe;
}




/**
void partition_overlaps(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, haplotype_evdience_alloc* hap,
                        int force_repeat)
{
    ResizeInitHaplotypeEvdience(hap);

    long long i;
    long long window_start, window_end;

    long long num_availiable_win = 0;
    
    Window_Pool w_inf;
    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;
    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;

        ///return overlaps that is overlaped with [window_start, window_end]
        flag = get_available_interval(window_start, window_end, overlap_list, dumy);
        switch (flag)
        {
            case 1:    ///found matched overlaps
                break;
            case 0:    ///do not find any matched overlaps
                break;
            case -2: ///do not find any matched overlaps, and the next window also cannot match
                break;
        }
        
        num_availiable_win = num_availiable_win + dumy->length;

        cluster(g_read->seq, window_start, window_end, overlap_list, dumy, R_INF, hap);
    }

    ///very time-consuming
    qsort(hap->list, hap->length, sizeof(haplotype_evdience), cmp_haplotype_evdience);

    
    

    ///debug_hap_information(overlap_list, R_INF, g_read, hap, dumy);

    SetSnpMatrix(hap, &(hap->nn_snp), &(overlap_list->length), 1, NULL);


    uint64_t pre_site = (uint64_t)-1;
    uint64_t num_of_snps = 0;
    long long pre_i = -1;
    long long sub_length;
    haplotype_evdience* sub_list;

    ////split reads
    for (i = 0; i < hap->length; i++)
    {
        if(pre_site != hap->list[i].site)
        {
            if(i != 0)
            {
                sub_list = hap->list + pre_i;
                sub_length = i - pre_i;
                split_sub_list(hap, sub_list, sub_length, overlap_list, R_INF, g_read);
            }
            num_of_snps++;
            pre_site = hap->list[i].site;
            pre_i = i;
        }
    }

    if(pre_i != -1)
    {
        sub_list = hap->list + pre_i;
        sub_length = i - pre_i;
        split_sub_list(hap, sub_list, sub_length, overlap_list, R_INF, g_read);
    }

    ///debug_snp_matrix(hap);
    generate_haplotypes_DP(hap, overlap_list, R_INF, g_read->length, force_repeat);
    ///generate_haplotypes_naive(hap, overlap_list, R_INF, g_read->length, force_repeat);

    lable_large_indels(overlap_list, g_read->length, dumy, asm_opt.max_ov_diff_ec);


    ///debug_snp_matrix(hap);
}
**/


inline void insert_snp_vv(haplotype_evdience_alloc* h, haplotype_evdience* a, uint64_t a_n, char misBase, UC_Read* g_read, void *km)
{
    if(a_n == 0) return;
    SnpStats *p = NULL; uint64_t /**nn = 0,**/ i;
    if(!km) kv_pushp(SnpStats, h->snp_stat, &p);
    else kv_pushp_km(km, SnpStats, h->snp_stat, &p);
    p->id = h->snp_stat.n-1;
    p->occ_0 = 1;
    p->occ_1 = 0;
    p->occ_2 = 0;
    p->overlap_num = 0;
    p->site = a[0].site;
    p->is_homopolymer = if_is_homopolymer_strict(p->site, g_read->seq, g_read->length);
    for (i = 0; i < a_n; i++) {
        if(a[i].type == 0) {
            a[i].overlapSite = p->id;
            h->snp_stat.a[p->id].occ_0 += a[i].cov;
        }
        else if(a[i].type == 1 && a[i].misBase == misBase) {
            a[i].overlapSite = p->id;
            h->snp_stat.a[p->id].occ_1 += a[i].cov;
        }
        else {
            h->snp_stat.a[p->id].occ_2 += a[i].cov;
        }
        h->snp_stat.a[p->id].overlap_num += a[i].cov;
    }
    h->snp_stat.a[p->id].score = -1;
}




















/**
void debug_phasing_status(overlap_region_alloc *olist, ma_ug_t *ug, uint64_t print_w_list,
haplotype_evdience_alloc* hap, UC_Read* g_read, int64_t flanking, uint64_t yid)
{
    uint64_t i, k, l, ii;
    int64_t t;
    SnpStats *s = NULL;
    
        for (i = 0; i < olist->length; i++) {
            if(olist->list[i].y_id != yid) continue;
            fprintf(stderr, "\n[M::utg%.6d%c::is_match->%u] rev->%u, x->[%u, %u), y->[%u, %u)\n", 
            (int)olist->list[i].y_id+1, "lc"[ug->u.a[olist->list[i].y_id].circ], olist->list[i].is_match, 
                olist->list[i].y_pos_strand, olist->list[i].x_pos_s, olist->list[i].x_pos_e+1, olist->list[i].y_pos_s, olist->list[i].y_pos_e+1);
            if(print_w_list) {
                for (k = 0; k < olist->list[i].w_list.n; k++) {
                    if(olist->list[i].w_list.a[k].y_end != -1) {
                        fprintf(stderr, "x->[%lu, %lu), y->[%d, %d), e->%d\n", 
                        olist->list[i].w_list.a[k].x_start, olist->list[i].w_list.a[k].x_end+1, 
                            olist->list[i].w_list.a[k].y_start, olist->list[i].w_list.a[k].y_end+1, 
                            olist->list[i].w_list.a[k].error);
                    } else {
                        fprintf(stderr, "x->[-1, -1), y->[-1, -1), e->-1\n");
                    }
                }
            }
        }
    

    for (k = 1, l = 0; k <= hap->length; ++k) {  
        if (k == hap->length || hap->list[k].overlapID != hap->list[l].overlapID) {
            ii = hap->list[l].overlapID;
            if(olist->list[ii].y_id != yid) {
                l = k;
                continue;
            }
            for (i = l; i < k; i++) {
                if(hap->list[i].type!=1) continue; 
                s = &(hap->snp_stat.a[hap->list[i].overlapSite]);
                if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2))) {
                    fprintf(stderr, "s->site:%u, s->occ_0:%u, s->occ_1:%u, s->occ_2:%u\n", s->site, s->occ_0, s->occ_1, s->occ_2);
                    for (t = s->site>=flanking?s->site-flanking:0; t<g_read->length && t<=s->site+flanking; t++){
                        if(t == s->site) fprintf(stderr,"[");
                        fprintf(stderr,"%c", g_read->seq[t]);
                        if(t == s->site) fprintf(stderr,"]");
                    }
                    fprintf(stderr,"\n");
                }
            }
            l = k;
        }
    }
}
**/

void align_ul_ed(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, char* qstr, char *tstr, double e_rate, int64_t w_l, void *km);



































/**
void recalcate_high_het_overlap(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    long long j, k, i;
    int threshold;
    long long y_id;
    int y_strand;
    long long y_readLen;
    long long x_start;
    long long x_end;
    long long x_len;
    long long total_y_start;
    long long total_y_end;
    long long y_start;
    long long Window_Len;
    char* x_string;
    char* y_string;
    int end_site;
    unsigned int error;
    int real_y_start;
    long long overlap_length;
    int extra_begin, extra_end;
    long long o_len;
    kvec_t(uint8_t) x_num;
    kvec_t(uint8_t) y_num;
    kv_init(x_num);
    kv_init(y_num);

    for (j = 0; j < (long long)overlap_list->length; j++)
    {
        
        if(overlap_list->list[j].w_list_length == 0) continue;
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);

        //i corresponding to each window of a overlap
        //utilize the the end pos of pre-window in backwards
        for (i = overlap_list->list[j].w_list_length - 1; i >= 0; i--)
        {
            ///the first matched window
            if(overlap_list->list[j].w_list[i].y_end != -1)
            {
                ///note!!! need notification
                ///this is the actual end postion in ystring
                total_y_start = overlap_list->list[j].w_list[i].y_end 
                            - overlap_list->list[j].w_list[i].extra_begin + 1;

                ///k corresponding to all unmatched windows at the right side of overlap_list->list[j].w_list[i]
                ///so k starts from i + 1, and end to the first matched window
                for (k = i + 1; k < (long long)overlap_list->list[j].w_list_length && overlap_list->list[j].w_list[k].y_end == -1; k++)
                {
                    extra_begin = extra_end = 0;

                    ///if y_start > y_readLen, direct terminate
                    if (total_y_start >= y_readLen)
                    {
                        break;
                    }
                    
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    ///there are two potiential reasons for unmatched window:
                    ///1. this window has a large number of differences
                    ///2. DP does not start from the right offset
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);
                    
                    y_start = total_y_start;
                    Window_Len = x_len + (threshold << 1);

                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, Get_READ_LENGTH((*R_INF), y_id), 
                    &extra_begin, &extra_end, &y_start, &o_len))
                    {
                        break;
                    }

                    if(o_len + threshold < x_len)
                    {
                        break;
                    }
                    
                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);

                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM(y_string, Window_Len, x_string, x_len, threshold, &error);
                    
                    ///if error==-1, unmatched
                    if (error!=(unsigned int)-1)
                    {
                        overlap_list->list[j].w_list[k].cigar.length = -1;
                        overlap_list->list[j].w_list[k].y_start = y_start;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = (int)error;
                        ///note!!! need notification
                        overlap_list->list[j].w_list[k].extra_begin = extra_begin;
                        overlap_list->list[j].w_list[k].extra_end = extra_end;
                        overlap_list->list[j].w_list[k].error_threshold = threshold;
                        
                        overlap_list->list[j].align_length += x_len;
                    }
                    else
                    {
                        break;
                    }

                    ///note!!! need notification
                    total_y_start = y_start + end_site - extra_begin + 1;
                }
                
            }
            
        }
        


        //i corresponding to each window of a overlap
        //utilize the the start pos of next window in forward
        for (i = 0; i < (long long)overlap_list->list[j].w_list_length; i++)
        {
            ///find the first matched window, which should not be the first window
            ///the pre-window of this matched window must be unmatched
            if(overlap_list->list[j].w_list[i].y_end != -1 && i != 0 && overlap_list->list[j].w_list[i - 1].y_end == -1)
            {
                ///check if the start pos of this matched window has been calculated
                if(overlap_list->list[j].w_list[i].cigar.length == -1)
                {
                    ///there is no problem for x
                    x_start = overlap_list->list[j].w_list[i].x_start;
                    x_end = overlap_list->list[j].w_list[i].x_end;
                    x_len = x_end - x_start + 1;
                    //may have bugs
                    threshold = overlap_list->list[j].w_list[i].error_threshold;
                    //may have bugs
                    //may have bugs
                    ///should not adjust threshold, since this window can be matched by the old threshold
                    ///threshold = Adjust_Threshold(threshold, x_len);
                    //may have bugs
                    Window_Len = x_len + (threshold << 1);


                    ///y_start is the real y_start
                    y_start = overlap_list->list[j].w_list[i].y_start;
                    extra_begin = overlap_list->list[j].w_list[i].extra_begin;
                    extra_end = overlap_list->list[j].w_list[i].extra_end;
                    o_len = Window_Len - extra_end - extra_begin;
                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, 
                    overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);


                    ///y_start has already been calculated
                    if (error != (unsigned int)-1)
                    {
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, 
                            end_site, extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, 
                            y_strand, error, &y_start, &real_y_start, &end_site, &extra_begin, 
                            &extra_end, &error))
                            {
                                overlap_list->list[j].w_list[i].error = error;
                                overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                overlap_list->list[j].w_list[i].extra_end = extra_end;
                            }
                        }
                                                 
                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);   

                        ///note!!! need notification
                        real_y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[i].y_start = real_y_start; 
                        ///I forget why don't reduce the extra_begin for y_end
                        ///it seems extra_begin will be reduced at the end of this function 
                        overlap_list->list[j].w_list[i].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[i].error = error;                         
                    }
                    else
                    {
                        fprintf(stderr, "error\n");
                    }
                }
                else
                {
                    real_y_start = overlap_list->list[j].w_list[i].y_start;
                }

                
                ///the end pos for pre window is real_y_start - 1
                total_y_end = real_y_start - 1;
                ///find the unmatched window on the left of current matched window
                ///k starts from i - 1
                for (k = i - 1; k >= 0 && overlap_list->list[j].w_list[k].y_end == -1; k--)
                {  
                    ///there is no problem in x
                    x_start = overlap_list->list[j].w_list[k].x_start;
                    x_end = overlap_list->list[j].w_list[k].x_end;
                    x_len = x_end - x_start + 1;
                    ///there are two potiential reasons for unmatched window:
                    ///1. this window has a large number of differences
                    ///2. DP does not start from the right offset
                    threshold = double_error_threshold(overlap_list->list[j].w_list[k].error_threshold, x_len);

                    Window_Len = x_len + (threshold << 1);

                    if(total_y_end <= 0)
                    {
                        break;
                    }

                    ///y_start might be less than 0
                    y_start = total_y_end - x_len + 1;
                    if(!determine_overlap_region(threshold, y_start, y_id, Window_Len, Get_READ_LENGTH((*R_INF), y_id),
                    &extra_begin, &extra_end, &y_start, &o_len))
                    {
                        break;
                    }

                    if(o_len + threshold < x_len)
                    {
                        break;
                    }

                    fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                    R_INF, y_id, extra_begin, extra_end);
                    x_string = g_read->seq + x_start;
                    y_string = dumy->overlap_region;

                    ///note!!! need notification
                    end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                    &(dumy->path_length), dumy->matrix_bit, dumy->path, -1, -1);

                    if (error!=(unsigned int)-1)
                    { 
                        ///this condition is always wrong
                        ///in best case, real_y_start = threshold, end_site = Window_Len - threshold - 1
                        if (end_site == Window_Len - 1 || real_y_start == 0)
                        {
                            fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                            extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                            &y_start, &real_y_start, &end_site,
                            &extra_begin, &extra_end, &error);
                        }

                        generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[k]),
                        &real_y_start, &end_site, &error, x_string, x_len, y_string);  

                        ///y_start has no shift, but y_end has shift           
                        overlap_list->list[j].w_list[k].y_start = y_start + real_y_start - extra_begin;
                        overlap_list->list[j].w_list[k].y_end = y_start + end_site;
                        overlap_list->list[j].w_list[k].error = error;
                        overlap_list->list[j].align_length += x_len;
                        overlap_list->list[j].w_list[k].extra_begin = extra_begin;
                        overlap_list->list[j].w_list[k].extra_end = extra_end;
                        overlap_list->list[j].w_list[k].error_threshold = threshold;
                    }
                    else
                    {
                        break;
                    }

                    total_y_end = y_start + real_y_start - 1 - extra_begin;
                }
            }
        }
    }



    overlap_list->mapped_overlaps_length = 0;
    
    double error_rate;
    int is_update = 0;
    for (j = 0; j < (long long)overlap_list->length; j++)
    {
        y_id = overlap_list->list[j].y_id;
        y_strand = overlap_list->list[j].y_pos_strand;
        y_readLen = Get_READ_LENGTH((*R_INF), y_id);
        overlap_length = overlap_list->list[j].x_pos_e - overlap_list->list[j].x_pos_s + 1;
        overlap_list->list[j].is_match = 0;
        is_update = 0;
        ///debug_scan_cigar(&(overlap_list->list[j]));
        if(overlap_list->list[j].w_list_length == 0 || overlap_length == 0 || overlap_list->list[j].align_length == 0) continue;
        ///only calculate cigar for high quality overlaps
        if (overlap_length * OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length)
        {
            
            for (i = 0; i < (long long)overlap_list->list[j].w_list_length; i++)
            {
                ///first we need to check if this window is matched
                if(overlap_list->list[j].w_list[i].y_end != -1)
                {
                    ///second check if the cigar of this window has been got 
                    if(overlap_list->list[j].w_list[i].cigar.length == -1)
                    {
                        ///there is no problem for x
                        x_start = overlap_list->list[j].w_list[i].x_start;
                        x_end = overlap_list->list[j].w_list[i].x_end;
                        x_len = x_end - x_start + 1;
                        //may have bugs
                        ///threshold = x_len * asm_opt.max_ov_diff_ec;
                        threshold = overlap_list->list[j].w_list[i].error_threshold;
                        //may have bugs
                        //may have bugs
                        ///should not adjust threshold, since this window can be matched by the old threshold
                        ///threshold = Adjust_Threshold(threshold, x_len);
                        //may have bugs
                        Window_Len = x_len + (threshold << 1);


                        ///y_start is the real y_start
                        ///for the window with cigar, y_start has already reduced extra_begin
                        y_start = overlap_list->list[j].w_list[i].y_start;
                        extra_begin = overlap_list->list[j].w_list[i].extra_begin;
                        extra_end = overlap_list->list[j].w_list[i].extra_end;
                        o_len = Window_Len - extra_end - extra_begin;
                        fill_subregion(dumy->overlap_region, y_start, o_len, y_strand, 
                        R_INF, y_id, extra_begin, extra_end);
                        x_string = g_read->seq + x_start;
                        y_string = dumy->overlap_region;


                        ///note!!! need notification
                        end_site = Reserve_Banded_BPM_PATH(y_string, Window_Len, x_string, x_len, threshold, &error, &real_y_start,
                        &(dumy->path_length), dumy->matrix_bit, dumy->path, 
                        overlap_list->list[j].w_list[i].error, overlap_list->list[j].w_list[i].y_end - y_start);

                        if (error != (unsigned int)-1)
                        {
                            if (end_site == Window_Len - 1 || real_y_start == 0)
                            {
                                
                                if(fix_boundary(x_string, x_len, threshold, y_start, real_y_start, end_site,
                                extra_begin, extra_end, y_id, Window_Len, R_INF, dumy, y_strand, error,
                                &y_start, &real_y_start, &end_site,
                                &extra_begin, &extra_end, &error))
                                {
                                    overlap_list->list[j].w_list[i].error = error;
                                    overlap_list->list[j].w_list[i].extra_begin = extra_begin;
                                    overlap_list->list[j].w_list[i].extra_end = extra_end;
                                }
                                
                            }

                            generate_cigar(dumy->path, dumy->path_length, &(overlap_list->list[j].w_list[i]),
                            &real_y_start, &end_site, &error, x_string, x_len, y_string);    

                            ///note!!! need notification
                            real_y_start = y_start + real_y_start - extra_begin;
                            overlap_list->list[j].w_list[i].y_start = real_y_start;  
                            overlap_list->list[j].w_list[i].y_end = y_start + end_site - extra_begin;
                            overlap_list->list[j].w_list[i].error = error;                              
                        }
                        else
                        {
                            fprintf(stderr, "error\n");
                        }
                    }
                    else
                    {
                        overlap_list->list[j].w_list[i].y_end -= overlap_list->list[j].w_list[i].extra_begin;
                    }


                }
            }

            error_rate = non_trim_error_rate(overlap_list, j, R_INF, dumy, g_read);
            
            
            if (error_rate <= HIGH_HET_ERROR_RATE)
            {
                is_update = 1;
            }
        }

        if((is_update == 0) && (overlap_list->list[j].align_length >= WINDOW) && 
           (overlap_length * HIGH_HET_OVERLAP_THRESHOLD_FILTER <=  overlap_list->list[j].align_length))
        {
            kv_resize(uint8_t, x_num, (uint64_t)(Get_READ_LENGTH((*R_INF), overlap_list->list[j].x_id)));
            kv_resize(uint8_t, y_num, (uint64_t)(Get_READ_LENGTH((*R_INF), overlap_list->list[j].y_id)));
            is_update = get_affine_gap_score(&(overlap_list->list[j]), g_read, overlap_read, x_num.a, y_num.a,
            overlap_list->list[j].x_pos_e + 1 - overlap_list->list[j].x_pos_s, 
            overlap_list->list[j].y_pos_e + 1 - overlap_list->list[j].y_pos_s);
        }

        if(is_update)
        {
            overlap_list->list[j].is_match = 2;
        }
    }

    kv_destroy(x_num);
    kv_destroy(y_num);
}



void correct_overlap_high_het(overlap_region_alloc* overlap_list, All_reads* R_INF, 
                        UC_Read* g_read, Correct_dumy* dumy, UC_Read* overlap_read)
{
    clear_Correct_dumy(dumy, overlap_list, NULL);

    long long window_start, window_end;

    Window_Pool w_inf;

    init_Window_Pool(&w_inf, g_read->length, WINDOW, (int)(1.0/asm_opt.max_ov_diff_ec));

    int flag = 0;

    while(get_Window(&w_inf, &window_start, &window_end) && flag != -2)
    {
        dumy->length = 0;
        dumy->lengthNT = 0;
        flag = get_interval(window_start, window_end, overlap_list, dumy, w_inf.window_length);

        switch (flag)
        {
            case 1:    ///no match here
                break;
            case 0:    ///no match here
                break;
            case -2: ///if flag == -2, loop would be terminated
                break;
        }

        ///dumy->lengthNT represent how many overlaps that the length of them is not equal to WINDOW; may larger or less than WINDOW
        ///dumy->length represent how many overlaps that the length of them is WINDOW
        //may improve
        ///now the windows which are larger than WINDOW are verified one-by-one, to improve it, we can do it group-bygroup
        verify_window(window_start, window_end, overlap_list, dumy, R_INF, g_read->seq);
    }

    recalcate_high_het_overlap(overlap_list, R_INF, g_read, dumy, overlap_read);
}
**/




///(char *qstr, kvec_t_u64_warp* q_idx) -> only used for hpc


///ts do not have aux_beg, while te has 




inline uint32_t gen_backtrace_adv_exz(window_list *p, overlap_region *z, All_reads *rref, hpc_t *hpc_g, const ul_idx_t *uref, 
char *qstr, char *tstr, bit_extz_t *exz, uint32_t rev, uint32_t id)
{
    if(p->error < 0 || p->y_end < 0) return 0;
    int64_t qs, qe, ql, tl, aln_l, t_pri_l, thres, ts, t_tot_l; 
    int64_t aux_beg, aux_end;
    char *q_string, *t_string; 
    ///there is no problem for x
    qs = p->x_start; qe = p->x_end; ql = qe + 1 - qs;
    thres = p->error_threshold; aln_l = ql + (thres<<1);

    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    ts = p->y_start; aux_beg = p->extra_begin; aux_end = p->extra_end;
    if(aux_end >= 0) {
        t_pri_l = aln_l - aux_beg - aux_end;
    } else {
        if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
        else if(uref) t_tot_l = uref->ug->u.a[id].len;
        else t_tot_l = Get_READ_LENGTH((*rref), id);

        t_pri_l = ts + aln_l - aux_beg; if(t_pri_l > t_tot_l) t_pri_l = t_tot_l;
        t_pri_l = t_pri_l - ts;
    }
    
    q_string = qstr + qs; tl = t_pri_l;
    if(rref) {
        recover_UC_Read_sub_region(tstr, ts, t_pri_l, rev, rref, id); t_string = tstr;
    } else {
        t_string = return_str_seq_exz(tstr, ts, t_pri_l, rev, hpc_g, uref, id);
    }

    exz->ts = 0; exz->te = p->x_end-p->x_start; exz->tl = ql;
    exz->ps = -1; exz->pe = p->y_end-p->y_start; exz->pl = tl;  
    exz->err = p->error; exz->thre = p->error_threshold; 
    // clear_align(*exz);
    ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thres, aux_beg, exz);
    // if(id == 178 && p->x_start == 86800 && p->x_end == 86807) {
    //     fprintf(stderr, "\n[M::%s::semi] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
    //     __func__, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
    //     fprintf(stderr, "[M::%s::semi] p->y_start::%d, p->y_end::%d, p->x_start::%d, p->x_end::%d, p->error::%d\n", 
    //     __func__, p->y_start, p->y_end, p->x_start, p->x_end, p->error);
    //     if(is_align(*exz)) {
    //         prt_cigar(exz->cigar.a, exz->cigar.n);
    //         fprintf(stderr, "[tstr] %.*s\n", exz->pe+1-exz->ps, t_string+exz->ps);
    //         fprintf(stderr, "[qstr] %.*s\n", exz->te+1-exz->ts, q_string+exz->ts);
    //     }
    // }
    // assert(is_align(*exz));
    // assert(cigar_check(t_string, q_string, exz));
    
            
    if(is_align(*exz)) {
        p->y_start = ts + exz->ps;///difference
        p->y_end = ts + exz->pe;
        p->error = exz->err;
        push_wcigar(p, &(z->w_list), exz);
        ///this condition is always wrong
        ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
        if ((((exz->pe+1) == tl) || (exz->ps == 0)) && (exz->err > 0)) {
            if(recal_boundary_exz(q_string, tstr, ql, tl, thres, ts, exz->ps, exz->pe,
            exz->err, id, rev, exz, rref, hpc_g, uref, &ts, &aux_beg, &aux_end)) {
                //update cigar
                z->w_list.c.n = p->cidx; push_wcigar(p, &(z->w_list), exz);

                p->y_start = ts + exz->ps;///difference
                p->y_end = ts + exz->pe;
                p->error = exz->err;
            }
        }
        p->extra_begin = aux_beg; 
        p->extra_end = aux_end;
        return 1;
    }
    p->error = -1;
    return 0;
}

inline uint32_t gen_backtrace_non_retrieve_adv_exz(window_list *p, overlap_region *z, char *qstr, char *tstr, int64_t t_tot_l, bit_extz_t *exz, uint32_t rev, uint32_t id)
{
    if(p->error < 0 || p->y_end < 0) return 0;
    int64_t qs, qe, ql, tl, aln_l, t_pri_l, thres, ts; 
    int64_t aux_beg, aux_end;
    char *q_string, *t_string; 
    ///there is no problem for x
    qs = p->x_start; qe = p->x_end; ql = qe + 1 - qs;
    thres = p->error_threshold; aln_l = ql + (thres<<1);

    ///y_start is the real y_start
    ///for the window with cigar, y_start has already reduced extra_begin
    ts = p->y_start; aux_beg = p->extra_begin; aux_end = p->extra_end;
    if(aux_end >= 0) {
        t_pri_l = aln_l - aux_beg - aux_end;
    } else {
        t_pri_l = ts + aln_l - aux_beg; if(t_pri_l > t_tot_l) t_pri_l = t_tot_l;
        t_pri_l = t_pri_l - ts;
    }
    
    q_string = qstr + qs; t_string = tstr + ts; tl = t_pri_l;

    exz->ts = 0; exz->te = p->x_end-p->x_start; exz->tl = ql;
    exz->ps = -1; exz->pe = p->y_end-p->y_start; exz->pl = tl;  
    exz->err = p->error; exz->thre = p->error_threshold; 
    // clear_align(*exz);
    ed_band_cal_semi_64_w_absent_diag_trace(t_string, tl, q_string, ql, thres, aux_beg, exz);
    // if(id == 178 && p->x_start == 86800 && p->x_end == 86807) {
    //     fprintf(stderr, "\n[M::%s::semi] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d\n", 
    //     __func__, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n);
    //     fprintf(stderr, "[M::%s::semi] p->y_start::%d, p->y_end::%d, p->x_start::%d, p->x_end::%d, p->error::%d\n", 
    //     __func__, p->y_start, p->y_end, p->x_start, p->x_end, p->error);
    //     if(is_align(*exz)) {
    //         prt_cigar(exz->cigar.a, exz->cigar.n);
    //         fprintf(stderr, "[tstr] %.*s\n", exz->pe+1-exz->ps, t_string+exz->ps);
    //         fprintf(stderr, "[qstr] %.*s\n", exz->te+1-exz->ts, q_string+exz->ts);
    //     }
    // }
    // assert(is_align(*exz));
    // assert(cigar_check(t_string, q_string, exz));
    
            
    if(is_align(*exz)) {
        p->y_start = ts + exz->ps;///difference
        p->y_end = ts + exz->pe;
        p->error = exz->err;
        push_wcigar(p, &(z->w_list), exz);
        ///this condition is always wrong
        ///in best case, r_ts = threshold, t_end = aln_l - thres - 1
        if ((((exz->pe+1) == tl) || (exz->ps == 0)) && (exz->err > 0)) {
            if(recal_boundary_non_retrieve_exz(q_string, tstr, t_tot_l, ql, tl, thres, ts, exz->ps, exz->pe,
            exz->err, id, rev, exz, &ts, &aux_beg, &aux_end)) {
                //update cigar
                z->w_list.c.n = p->cidx; push_wcigar(p, &(z->w_list), exz);

                p->y_start = ts + exz->ps;///difference
                p->y_end = ts + exz->pe;
                p->error = exz->err;
            }
        }
        p->extra_begin = aux_beg; 
        p->extra_end = aux_end;
        return 1;
    }
    p->error = -1;
    return 0;
}


///ts do not have aux_beg, while te has 

#define pass_qovlp(o, a, r) (((a)>0)&&((o)*(r)<=(a)))

///ts do not have aux_beg, while te has 
uint32_t push_hc_wlst_exz(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, overlap_region* ol, 
                        char* qstr, char *tstr, bit_extz_t *exz, uint32_t max_err,
                        int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t tl, 
                        int64_t aux_beg, int64_t aux_end, double e_rate, int64_t block_s, double ovlp_cut, int64_t force_aln, void *km)
{

    window_list p, t, *a; int64_t w_e, w_s, ce = qs - 1, cs = ol->x_pos_s, toff, ovl, ualn, aln, ys; 
    uint64_t a_n, k;
    p.x_start = qs; p.x_end = qe; p.y_start = ts; p.y_end = te; p.error = exz->err;
    p.extra_begin = aux_beg; p.extra_end = aux_end; p.error_threshold = exz->thre; p.cidx = p.clen = 0;
    if(ol->w_list.n > 0) { //utilize the the end pos of pre-window in forward
        w_e = ol->w_list.a[ol->w_list.n-1].x_end; 
        toff = ol->w_list.a[ol->w_list.n-1].y_end + 1;
        while ((w_e < ce) && (toff < tl)) {
            w_s = w_e + 1;
            get_win_id_by_s(ol, w_s, block_s, &w_e);
            // x_start = w_s; x_end = w_e;
            if(aln_wlst_adv_exz(ol, rref, hpc_g, uref, qstr, tstr, exz, max_err,
            ol->y_pos_strand, ol->y_id, w_s, w_e, toff, block_s, e_rate, 0)) {
                toff = ol->w_list.a[ol->w_list.n-1].y_end + 1;
            } else {
                break;
            }
        }
        cs = ol->w_list.a[ol->w_list.n-1].x_end + 1;
    }
    ///utilize the the start pos of next window in backward
    a_n = ol->w_list.n; w_s = qs;
    if(w_s > cs) {
        gen_backtrace_adv_exz(&p, ol, rref, hpc_g, uref, qstr, tstr, exz, ol->y_pos_strand, ol->y_id);
        toff = p.y_start - 1;
        while (w_s > cs) {
            w_e = w_s - 1;
            get_win_id_by_e(ol, w_e, block_s, &w_s); ys = toff+1-(w_e+1-w_s);
            // x_start = w_s; x_end = w_e; x_len = x_end + 1 - x_start;
            if((ys >= 0) && aln_wlst_adv_exz(ol, rref, hpc_g, uref, qstr, tstr, exz, max_err,
                ol->y_pos_strand, ol->y_id, w_s, w_e, ys, block_s, e_rate, 1)) {
                toff = ol->w_list.a[ol->w_list.n-1].y_start - 1;
            } else {
                break;
            }
        }
    }

    ol->align_length += qe + 1 - qs;
    ovl = ol->x_pos_e+1-ol->x_pos_s; ualn = (qe + 1 - ol->x_pos_s) - ol->align_length; aln = ovl-ualn;
    // if((!force_aln) && (!simi_pass(ovl, aln, 0, ovlp_cut, &e_rate)) && (!simi_pass(ovl, aln, sec_check, ovlp_cut, NULL))) {
    if((!force_aln) && (!pass_qovlp(ovl, aln, ovlp_cut))) {
        kv_push(window_list, ol->w_list, p); 
        return 0;
    }

    if(ol->w_list.n > a_n) {
        a = ol->w_list.a + a_n; a_n = ol->w_list.n - a_n; toff = a_n; a_n >>=1;
        for (k = 0; k < a_n; k++) {
            t = a[k]; a[k] = a[toff-1-k]; a[toff-1-k] = t;
        }
    }
    kv_push(window_list, ol->w_list, p); 
    return 1;
}


///ts do not have aux_beg, while te has 



uint32_t align_hc_ed_post_extz(overlap_region *z, All_reads *rref, char* qstr, char *tstr, bit_extz_t *exz, double e_rate, int64_t w_l, double ovlp_cut, int64_t force_aln, void *km)
{
    int64_t q_s, q_e, nw, k, q_l, t_l, t_tot_l, aux_beg, aux_end, t_s, thre, aln_l, t_pri_l;
    char *q_string, *t_string; 
    z->w_list.n = 0; z->is_match = 0; z->align_length = 0;
    nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, w_l);
    get_win_se_by_normalize_xs(z, (z->x_pos_s/w_l)*w_l, w_l, &q_s, &q_e);
    // if(z->x_id == 19350 && z->y_id == 19324) {
    //     fprintf(stderr, "-z-[M::%s] tid::%u(%c)\tq::[%u,%u)\tt::[%u,%u)\n", __func__, z->y_id, "+-"[z->y_pos_strand], z->x_pos_s, z->x_pos_e+1, z->y_pos_s, z->y_pos_e+1);
    // }
    for (k = 0; k < nw; k++) {
        aux_beg = aux_end = 0; q_l = 1 + q_e - q_s;
        thre = q_l*e_rate; thre = Adjust_Threshold(thre, q_l);
        if(thre > THRESHOLD_MAX_SIZE) thre = THRESHOLD_MAX_SIZE;
        ///offset of y
        t_s = (q_s - z->x_pos_s) + z->y_pos_s;
        t_s += y_start_offset(q_s, &(z->f_cigar));

        aln_l = q_l + (thre<<1); t_tot_l = Get_READ_LENGTH((*rref), z->y_id); 
        if(init_waln(thre, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) {
            q_string = qstr+q_s; 
            recover_UC_Read_sub_region(tstr, t_s, t_pri_l, z->y_pos_strand, rref, z->y_id); t_string = tstr;
            // t_string = return_str_seq_exz(tstr, t_s, t_pri_l, z->y_pos_strand, hpc_g, uref, z->y_id);
            
            t_l = t_pri_l;
            // t_end = Reserve_Banded_BPM(t_string, aln_l, q_string, q_l, thre, &error);
            ed_band_cal_semi_64_w_absent_diag(t_string, t_l, q_string, q_l, thre, aux_beg, exz);

            // if(z->x_id == 5569 && z->y_id == 5557 && q_s == 10075 && q_e == 10849) {
            //     fprintf(stderr, "\n[M::%s::semi::t_s->%ld::t_pri_l->%ld::aux_beg->%ld::aux_end->%ld::thre->%ld] exz->ps::%d, exz->pe::%d, exz->ts::%d, exz->te::%d, exz->err::%d, exz->cigar.n::%d, thre::%ld\n", 
            //     __func__, t_s, t_pri_l, aux_beg, aux_end, thre, exz->ps, exz->pe, exz->ts, exz->te, exz->err, (int32_t)exz->cigar.n, thre);
            //     fprintf(stderr, "[tstr::len->%ld] %.*s\n", t_l, (int32_t)t_l, t_string);
            //     fprintf(stderr, "[qstr::len->%ld] %.*s\n", q_l, (int32_t)q_l, q_string);
            // } 
            if (is_align(*exz)) {
                // ed_band_cal_semi_64_w(t_string, aln_l, q_string, q_l, thre, exz);
                // assert(exz->err <= exz->thre);
                // if(z->x_id == 19350 && z->y_id == 19324) {
                //     fprintf(stderr, "+[M::%s]\tq::[%ld,%ld)\tt::[%ld,%ld)\texz->err::%d\n", __func__, q_s, q_e + 1, t_s, t_s + exz->pe + 1, exz->err);
                // }
                ///t_s do not have aux_beg, while t_s + t_end (aka, te) has 
                if(!push_hc_wlst_exz(NULL, NULL, rref, z, qstr, tstr, exz, THRESHOLD_MAX_SIZE, q_s, q_e, t_s, t_s + exz->pe, 
                                        t_tot_l, aux_beg, aux_end, e_rate, w_l, ovlp_cut, force_aln, km)) {
                    return 0;
                }
                // append_window_list(z, q_s, q_e, t_s, t_s + t_end, error, aux_beg, aux_end, thre, w_l, km);
            } 
            // else {
                // if(z->x_id == 19350 && z->y_id == 19324) {
                //     fprintf(stderr, "-[M::%s]\tq::[%ld,%ld)\tt::[%ld,%ld)\texz->err::%d\n", __func__, q_s, q_e + 1, t_s, t_s + exz->pe + 1, exz->err);
                // }
            // }
        }
        q_s = q_e + 1; q_e = q_s + w_l - 1; 
        if(q_e >= (int64_t)z->x_pos_e) q_e = z->x_pos_e;
    }

    // if((!force_aln) && (!simi_pass(z->x_pos_e+1-z->x_pos_s, z->align_length, 0, ovlp_cut, &e_rate)) &&
    //     (!simi_pass(z->x_pos_e+1-z->x_pos_s, z->align_length, sec_check, ovlp_cut, NULL))) return 0;
    if((!force_aln) && (!pass_qovlp(z->x_pos_e+1-z->x_pos_s, z->align_length, ovlp_cut))) return 0;
    return 1;
}


inline uint32_t ed_cut(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char *qstr, char *tstr, uint32_t rev, uint32_t id, 
int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, double e_rate, int64_t max_err,
uint32_t aln_dir, int64_t* r_err, int64_t* qoff, int64_t* toff, int64_t* aln_qlen)
{
    (*aln_qlen) = 0; (*r_err) = INT32_MAX;
    if(qoff) (*qoff) = -1; if(toff) (*toff) = -1; 
    int64_t ql, aln_l, t_tot_l, aux_beg, aux_end, t_pri_l, thres;
    char *q_string, *t_string; unsigned int error; int t_end, q_end;

    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    if(rref) {
        thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, max_err), ql);
    } else {
        thres = double_ul_error_threshold(get_init_err_thres(ql, e_rate, block_s, max_err), ql);
    }
    
    aln_l = ql + (thres << 1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    // if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; 
    if(rref) {
        fill_subregion(tstr, t_s, t_pri_l, rev, rref, id, aux_beg, aux_end); t_string = tstr;
    } else {
        t_string = return_str_seq(tstr, t_s, t_pri_l, rev, hpc_g, uref, id, aux_beg, aux_end);
    }

    // if(id == 6) {
    //     fprintf(stderr, "-[M::%s::aln_dir->%u] qs->%ld, ts->%ld, thres->%ld, aux_beg->%ld, aux_end->%ld, t_pri_l->%ld\n", 
    //     __func__, aln_dir, qs, t_s, thres, aux_beg, aux_end, t_pri_l);
    // }
    if(aln_dir == 0) {
        Reserve_Banded_BPM_Extension(t_string, aln_l, q_string, ql, thres, &error, &t_end, &q_end);
    } else {
        Reserve_Banded_BPM_Extension_REV(t_string, aln_l, q_string, ql, thres, &error, &t_end, &q_end);
    }

    if(t_end != -1 && q_end != -1) (*aln_qlen) = (aln_dir?(ql-q_end):(q_end+1));
    if(qoff) (*qoff) = q_end; if(toff) (*toff) = t_end; (*r_err) = error;

    if((*aln_qlen) == 0) return 0;
    return 1;
}

inline uint32_t ed_non_retrieve_cut(char *qstr, char *tstr, int64_t t_tot_l, uint32_t rev, uint32_t id, 
int64_t qs, int64_t qe, int64_t t_s, int64_t block_s, double e_rate, int64_t max_err,
uint32_t aln_dir, int64_t* r_err, int64_t* qoff, int64_t* toff, int64_t* aln_qlen)
{
    (*aln_qlen) = 0; (*r_err) = INT32_MAX;
    if(qoff) (*qoff) = -1; if(toff) (*toff) = -1; 
    int64_t ql, aln_l, aux_beg, aux_end, t_pri_l, thres;
    char *q_string, *t_string; unsigned int error; int t_end, q_end;

    ql = qe + 1 - qs;
    ///there are two potiential reasons for unmatched window:
    ///1. this window has a large number of differences
    ///2. DP does not start from the right offset
    thres = double_error_threshold(get_init_err_thres(ql, e_rate, block_s, max_err), ql);
    
    aln_l = ql + (thres << 1);

    if(!init_waln(thres, t_s, t_tot_l, aln_l, &aux_beg, &aux_end, &t_s, &t_pri_l)) return 0;
    // if(t_pri_l + thres < ql) return 0;

    q_string = qstr + qs; t_string = tstr + t_s;

    // if(id == 6) {
    //     fprintf(stderr, "-[M::%s::aln_dir->%u] qs->%ld, ts->%ld, thres->%ld, aux_beg->%ld, aux_end->%ld, t_pri_l->%ld\n", 
    //     __func__, aln_dir, qs, t_s, thres, aux_beg, aux_end, t_pri_l);
    // }
    if(aln_dir == 0) {
        Reserve_Banded_BPM_Extension(t_string, aln_l, q_string, ql, thres, &error, &t_end, &q_end);
    } else {
        Reserve_Banded_BPM_Extension_REV(t_string, aln_l, q_string, ql, thres, &error, &t_end, &q_end);
    }

    if(t_end != -1 && q_end != -1) (*aln_qlen) = (aln_dir?(ql-q_end):(q_end+1));
    if(qoff) (*qoff) = q_end; if(toff) (*toff) = t_end; (*r_err) = error;

    if((*aln_qlen) == 0) return 0;
    return 1;
}



int64_t gen_extend_err_0_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, bit_extz_t *exz, uint64_t *v_idx, int64_t block_s, double e_rate, int64_t max_err,
int64_t qs, int64_t qe, int64_t pk)
{
    int64_t tot_e = 0, ts, di[2], al[2], tb[2], an = z->w_list.n; double rr;
    int64_t id = z->y_id, rev = z->y_pos_strand, ql = qe + 1 - qs;
    ///check if there are some windows that cannot be algined by any overlaps/unitigs
    ///if no, it is likely that the UL read itself has issues
    if(uref && v_idx && z->is_match == 4) {
        if(check_coverage_gap(v_idx, qs, qe, block_s)) {
            tot_e += THRESHOLD_MAX_SIZE; return tot_e;
        } 
    }
    ts = (qs - z->x_pos_s) + z->y_pos_s; ts += y_start_offset(qs, &(z->f_cigar));

    di[0] = di[1] = al[0] = al[1] = 0; tb[0] = tb[1] = -1;
    if((pk > 0) && (qs == (z->w_list.a[pk].x_end + 1))) {
        if(z->w_list.a[pk].clen == 0) {///do not have cigar
            gen_backtrace_adv_exz(&(z->w_list.a[pk]), z, rref, hpc_g, uref, qstr, tstr, exz, rev, id);
        }
        tb[0] = z->w_list.a[pk].y_end + 1;
    }

    if(((pk+1) < an) && ((qe+1) == (z->w_list.a[pk+1].x_start))) {
        if(z->w_list.a[pk+1].clen == 0) {///do not have cigar
            gen_backtrace_adv_exz(&(z->w_list.a[pk+1]), z, rref, hpc_g, uref, qstr, tstr, exz, rev, id);
        }
        tb[1] = z->w_list.a[pk+1].y_start-ql;
    }

    if(tb[0] == -1 && tb[1] == -1) tb[0] = tb[1] = ts;
    else if(tb[0] == -1 && tb[1] != -1) tb[0] = tb[1];
    else if(tb[1] == -1 && tb[0] != -1) tb[1] = tb[0];

    if(tb[0] != -1) {
        if(!ed_cut(uref, hpc_g, rref, qstr, tstr, rev, id, qs, qe, tb[0], block_s, e_rate, max_err,
                                                                    0, &(di[0]), NULL, NULL, &(al[0]))) {
            di[0] = ql; al[0] = 0;
        }
    }

    if(tb[1] != -1) {
        if(!ed_cut(uref, hpc_g, rref, qstr, tstr, rev, id, qs, qe, tb[1], block_s, e_rate, max_err,
                                                                    1, &(di[1]), NULL, NULL, &(al[1]))) {
            di[1] = ql; al[1] = 0;
        }
    }

    if(al[0] && al[1]) {///matched in both sides
        if((al[0] + al[1]) <= ql) {
            tot_e += di[0] + di[1] + ql - (al[0] + al[1]);
        } else {
            rr = ((double)ql)/((double)(al[0] + al[1]));
            tot_e += (di[0] + di[1])*rr;
        }
    } else if((!al[0]) && (!al[1])) {//failed
        tot_e += ql;
    } else if(al[0]) {
        tot_e += di[0] + (ql - al[0]);
    }else if(al[1]) {
        tot_e += di[1] + (ql - al[1]);
    }
    // if(z->y_id == 6) {
    //     fprintf(stderr, "-[M::%s] qs->%ld, ts->%ld, tb[0]->%ld, tb[1]->%ld, di[0]->%ld, di[1]->%ld, al[0]->%ld, al[1]->%ld, block_s->%ld, e_rate->%f\n", __func__, 
    //     qs, ts, tb[0], tb[1], di[0], di[1], al[0], al[1], block_s, e_rate);
    // }
    return tot_e;
}



double gen_extend_err_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, 
char *tstr, bit_extz_t *exz, uint64_t *v_idx, int64_t block_s, double ovlp_cut, double e_rate, double e_max, int64_t max_err, int64_t sec_check, int64_t *r_e)
{
    int64_t ovl, k, ce, an = z->w_list.n, tot_l, tot_e, ws, we, ql; 
    ovl = z->x_pos_e+1-z->x_pos_s; if(r_e) (*r_e) = INT64_MAX;
    if((sec_check) && (!simi_pass(ovl, z->align_length, 0, ovlp_cut, &e_rate))) return DBL_MAX;
    // nw = get_num_wins(z->x_pos_s, z->x_pos_e+1, block_s);
    // for (k = 0; k < an; k++) {
    //     if(z->w_list.a[k].clen) z->w_list.a[k].y_end -= z->w_list.a[k].extra_begin;
    // }

    tot_l = tot_e = 0;
    for (k = an-1, ce = z->x_pos_e; k >= 0; k--) {
        // assert(k == 0 || z->w_list.a[k].x_end > z->w_list.a[k-1].x_start);//sorted
        tot_l += z->w_list.a[k].x_end + 1 - z->w_list.a[k].x_start;
        tot_e += z->w_list.a[k].error;///matched window

        we = z->w_list.a[k].x_end;
        while (we < ce) {
            ws = we+1;
            get_win_id_by_s(z, ws, block_s, &we);
            ql = we+1-ws; tot_l += ql; 
            tot_e += gen_extend_err_0_exz(z, uref, hpc_g, rref, qstr, tstr, exz, v_idx, block_s, e_rate, max_err, ws, we, k);
            if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
        }
        ce = z->w_list.a[k].x_start-1;
        if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
    }

    if(ce >= ((int64_t)z->x_pos_s)) {
        we = ((int64_t)z->x_pos_s)-1;
        while (we < ce) {
            ws = we+1;
            get_win_id_by_s(z, ws, block_s, &we);
            ql = we+1-ws; tot_l += ql; 
            tot_e += gen_extend_err_0_exz(z, uref, hpc_g, rref, qstr, tstr, exz, v_idx, block_s, e_rate, max_err, ws, we, k);
            if((e_max > 0) && (tot_e > (ovl*e_max))) return DBL_MAX;
        }
    }

    assert(tot_l == ovl); if(r_e) (*r_e) = tot_e;
    return (double)(tot_e)/(double)(tot_l);
}




#define gen_hpc_max_len(x) ((x)+((x)>>1)+1)
///[off_s, off_e)


#define cl_pushp(type, v, p) do {									\
		if ((v).length == (v).size) {										\
			(v).size = (v).size? (v).size<<1 : 2;							\
			(v).list = (type*)realloc((v).list, sizeof(type) * (v).size);	\
		}															\
		*(p) = &(v).list[(v).length++]; \
	} while (0)

///ai is the suffix of aj










// void ul_phase(overlap_region *oa, int64_t on, uint64_t *idx, uint64_t *buf)
// {
//     int64_t i, oi; overlap_region *z;
//     for (i = 0; i < on; i++) {
//         oi = (uint32_t)idx[i]; z = &(oa[oi]);
//         if(z->is_match != 1) continue;
//     }
// }

inline uint32_t cigar_check_dbg(char *pstr, char *tstr, bit_extz_t *ez)
{
	int32_t pi = ez->ps, ti = ez->ts, err = 0; uint32_t ci = 0, cl, k; uint16_t c;
	while (ci < ez->cigar.n) {
		ci = pop_trace(&(ez->cigar), ci, &c, &cl);
		// fprintf(stderr, "# %u = %u, cigar_n::%u\n", c, cl, (uint32_t)ez->cigar.n); 
		if(c == 0) {
			for (k=0;(k<cl)&&(pstr[pi]==tstr[ti]);k++,pi++,ti++);
			if(k!=cl) {
				fprintf(stderr, "ERROR-d-0\n"); 
				return 0;
			}
		} else {
			err += cl;
			if(c == 1) {
				for (k=0;(k<cl)&&(pstr[pi]!=tstr[ti]);k++,pi++,ti++);
				if(k!=cl) {
					fprintf(stderr, "ERROR-d-1\n"); 
					return 0;
				}
			} else if(c == 2) {///more p
				pi+=cl;
			} else if(c == 3) {
				ti+=cl;
			}
		}
	}
	if(err != ez->err) {
		fprintf(stderr, "ERROR-err, err::%d, ez->err::%d\n", err, ez->err); 
		return 0;
	}
	return 1;
}





inline uint64_t scale_ed_thre(uint32_t err, uint32_t max_err)
{
    uint64_t bd = (err<<1)+1, w; 
    w = (bd>>bitw); w <<= bitw; if(w < bd) w += bitwbit;
    err = (w-1)>>1; if(err > max_err) err = max_err;
    return err;
}


///[qs, qe)
int64_t update_semi_coord(const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, overlap_region *z, 
int64_t qs, int64_t qe, int64_t thre, int64_t *ts, int64_t *te, int64_t *aux_beg)
{
    int64_t ql = qe - qs, aln_l, t_tot_l, id = z->y_id, aux_end, tl; 
    (*ts) = (qs - z->x_pos_s) + z->y_pos_s;
    (*ts) += y_start_offset(qs, &(z->f_cigar));
    aln_l = ql + (thre<<1);
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);
    if(!init_waln(thre, (*ts), t_tot_l, aln_l, aux_beg, &aux_end, ts, &tl)) {
        (*ts) = (*te) = (*aux_beg) = -1;
        return 0;
    }
    (*te) = (*ts) + tl; 
    return 1;
}


///[qs, qe)


void adjust_ext_offset(int64_t *qs, int64_t *qe, int64_t *ts, int64_t *te, int64_t ql, int64_t tl, int64_t thre, int64_t mode)
{
    int64_t qoff, toff;
    if(mode == 1) {///forward extension
        qoff = ql - (*qs); toff = tl - (*ts);
        if(qoff <= toff) {
            (*qe) = ql; (*te) = (*ts) + qoff + thre;
        } else {
            (*te) = tl; (*qe) = (*qs) + toff + thre;
        }
    } else if(mode == 2) {///backward extension
        qoff = (*qe); toff = (*te);
        if(qoff <= toff) {
            (*qs) = 0; (*ts) = (*te) - qoff - thre;
        } else {
            (*ts) = 0; (*qs) = (*qe) - toff - thre;
        }
    }
    if((*qs) < 0) (*qs) = 0;
    if((*ts) < 0) (*ts) = 0;
    if((*qe) > ql) (*qe) = ql;
    if((*te) > tl) (*te) = tl;
}
















///[ys, ye)
uint64_t inline get_win_aln(overlap_region *z, uint64_t wid, int64_t *ys, int64_t *ye, int64_t *err)
{
	(*err) = -2;
	if((wid > 0) && (z->w_list.a[wid].y_end != -1) && (z->w_list.a[wid-1].y_end != -1) && (z->w_list.a[wid].y_end > z->w_list.a[wid-1].y_end)) {
		(*ys) = z->w_list.a[wid-1].y_end+1; 
		(*ye) = z->w_list.a[wid].y_end+1; 
		(*err) = z->w_list.a[wid].error;
		return 1;
	}
	return 0;
}



inline int64_t khit_long_gap(k_mer_hit *a, k_mer_hit *b, double small_bw_rate, int64_t min_small_bw)
{
    int64_t dq, dr, dd, dm;
    dq = b->self_offset-a->self_offset;
    dr = b->offset-a->offset;
	dd = dq>=dr? ((dq)-(dr)): ((dr)-(dq));

    dm = dq>=dr?dr:dq;
    if((dd > (dm*small_bw_rate)) && (dd > min_small_bw)) return 0;
    return 1;
}















int64_t cal_estimate_err_hc(overlap_region *z, int64_t wl, int64_t qs, int64_t qe, int64_t ts, int64_t te, double e_rate, int64_t *exact)
{
    int64_t k, ws, we, wid, os, oe, ovlp, tot, cov_l, est = (qe-qs)*e_rate, wn = z->w_list.n, exa = 1, ots, ote, q[2], t[2]; 
    if(exact) (*exact) = 0;
    if(!wn) return est;
    if(qs < z->x_pos_s) qs = z->x_pos_s;
    if(qe > z->x_pos_e+1) qe = z->x_pos_e+1;
    ws = qs/wl; ws *= wl; wid = get_win_id_by_s(z, ws, wl, NULL);
    if(wid>=wn) wid = wn-1;
    for (k = wid;k < wn && qs > z->w_list.a[k].x_end; k++); if(k == wn) return est;
    ///qs <= z->w_list.a[k].x_end
    for (;k>=0 && qs < z->w_list.a[k].x_start; k--); if(k < 0) k = 0;
    ///qs >= z->w_list.a[k].x_start

    for (tot = cov_l = 0, ots = ote = -1; k < wn && z->w_list.a[k].x_start < qe; k++) {
        if(z->w_list.a[k].y_end == -1) continue;
        ws = z->w_list.a[k].x_start; we = z->w_list.a[k].x_end+1;
        os = MAX(qs, ws); oe = MIN(qe, we);
	    ovlp = ((oe>os)? (oe-os):0);
        if(!ovlp) continue;
        cov_l += ovlp;
        
        if(ovlp == (we-ws)) {
            tot += z->w_list.a[k].error;
        } else {
            tot += ((double)z->w_list.a[k].error)*((double)ovlp)/((double)(we-ws));
        }
        if(z->w_list.a[k].error > 0) exa = 0;
        if(exa) {
            q[0] = os - ws; q[1] = we - oe; 
            we = z->w_list.a[k].y_end+1; ws = we - (z->w_list.a[k].x_end+1-z->w_list.a[k].x_start);
            os = MAX(ts, ws); oe = MIN(te, we);
	        ovlp = ((oe>os)? (oe-os):0);
            t[0] = os - ws; t[1] = we - oe; 
            if((ovlp) && (q[0] == t[0]) && (t[0] == t[0])) {
                if(ote == -1) {
                    ots = os; ote = oe; 
                } else if(ote == os) {
                    ote = oe; 
                } else {
                    exa = 0;
                }
            } else {
                exa = 0;
            }
        }
    }
    tot += ((qe-qs)-cov_l)*e_rate; 
    if((exact) && exa && (((qe-qs) == cov_l))) {
        if(((qe-qs) == (ote - ots)) && (ots == ts) && (ote == te)) (*exact) = 1;
    }
    return tot;
}


///[s, e); [ps, pe)
inline char *retrieve_str_seq_exz(UC_Read *tu, int64_t s, int64_t l, 
int64_t ps, int64_t pl, uint8_t rev, const ul_idx_t *uref, hpc_t *hpc_g, 
All_reads *rref, int64_t id)
{
    if(!hpc_g) {
        char *str; int64_t ss = s, sl = l; tu->length = l;
        UC_Read_resize(*tu, sl); str = tu->seq; 
        if(s == ps) {
            if(l <= pl) return tu->seq;
            str = tu->seq + pl; ss = ps + pl; sl = l - pl;
        }
        if(uref) {
            retrieve_u_seq(NULL, str, &(uref->ug->u.a[id]), rev, ss, sl, NULL);
        } else if(rref) {
            recover_UC_Read_sub_region(str, ss, sl, rev, rref, id);
        }
        return tu->seq;
    } else {
        return hpc_str(*hpc_g, id, rev) + s;
    }
}

void cal_exz_global(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_global_64_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_global_128_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_global_infi_w_trace(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}


void cal_exz_extension_0(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_extension_64_0_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_extension_128_0_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_extension_infi_0_w_trace(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}


void cal_exz_extension_1(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_extension_64_1_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else if(nword == 2) {
        ed_band_cal_extension_128_1_w_trace(pstr, pn, tstr, tn, thre, ez);
    } else {
        ed_band_cal_extension_infi_1_w_trace(pstr, pn, tstr, tn, thre, &nword, ez);
    }
}


void cal_exz_semi(char *pstr, int32_t pn, char *tstr, int32_t tn, int32_t thre, int32_t aux_beg, bit_extz_t *ez)
{
    int32_t bd, nword;
    bd = (((thre)<<1)+1); nword = ((bd>>bitw)+(!!(bd&bitz)));

    if(nword <= 1) {
        ed_band_cal_semi_64_w_absent_diag_trace(pstr, pn, tstr, tn, thre, aux_beg, ez);
    } else if(nword == 2) {
        ed_band_cal_semi_128_w_absent_diag_trace(pstr, pn, tstr, tn, thre, aux_beg, ez);
    } else {
        ed_band_cal_semi_infi_w_absent_diag_trace(pstr, pn, tstr, tn, thre, aux_beg, &nword, ez);
    }
}



int64_t cal_exz_infi_adv(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t *pts, int64_t *pte, int64_t thre, int64_t *pthre, int64_t q_tot_l, int64_t mode)
{
    clear_align(*exz);
    int64_t aux_beg = 0, ql, tl, t_tot_l = -1, dd; 
    char *q_string, *t_string; int32_t rev = z->y_pos_strand, id = z->y_id; 
    ql = qe - qs; tl = te - ts; dd = MAX(ql, tl); 
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);
   
    if(mode == 3) {
        update_semi_coord(uref, hpc_g, rref, z, qs, qe, ((thre>dd)?dd:thre), &ts, &te, &aux_beg);
    } else if(mode == 1 || mode == 2) {
        adjust_ext_offset(&qs, &qe, &ts, &te, q_tot_l, t_tot_l, ((thre>dd)?dd:thre), mode);
    }
    
    if((qe > qs) && (te > ts) && (ts != -1) && (te != -1)) {
        ql = qe - qs; tl = te - ts; 
        dd = MAX(ql, tl); 
        if(thre > dd) thre = dd; 
        if(thre <= (*pthre)) return 0;
        (*pthre) = thre;
        
        q_string = qstr + qs; 
        t_string = retrieve_str_seq_exz(tu, ts, tl, (*pts), (*pte)-(*pts), rev, uref, hpc_g, rref, id);
        (*pts) = ts; (*pte) = te;

        if(mode == 0) { //global
            cal_exz_global(t_string, tl, q_string, ql, thre, exz);
        } else if(mode == 1) {///forward extension
            cal_exz_extension_0(t_string, tl, q_string, ql, thre, exz);
        } else if(mode == 2) {///backward extension
            cal_exz_extension_1(t_string, tl, q_string, ql, thre, exz);
        } else if(mode == 3) {//semi-global
            cal_exz_semi(t_string, tl, q_string, ql, thre, aux_beg, exz);
        }

        if(is_align(*exz)) {
            // cigar_check(t_string, q_string, exz); 
            // if(mode == 1) {
            //     fprintf(stderr, "\n[M::%s::ql::%ld] qs::%ld, qe::%ld, ts::%ld, te::%ld, mode::%ld, err::%d, thre::%d, exz_q[%d, %d], exz_t[%d, %d]\n", 
            //                      __func__, ql, qs, qe, ts, te, mode, exz->err, exz->thre, exz->ts, exz->te, exz->ps, exz->pe);
            //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tl, t_string);
            //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)ql, q_string);
            //     fprintf(stderr, "[M::%s::] exz->cigar.n::%d\n", __func__, (int32_t)exz->cigar.n);
            // }
            exz->ps += ts; exz->pe += ts;
            exz->ts += qs; exz->te += qs;
            return 1;
        }
        return 0;
    }
    return 0;
}




int64_t cal_exact_exz(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
bit_extz_t *exz, char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t *pts, int64_t *pte, int64_t q_tot_l, int64_t mode)
{
    clear_align(*exz); exz->thre = 0; exz->cigar.n = 0;
    int64_t ql, tl, t_tot_l = -1; 
    char *q_string, *t_string; int32_t rev = z->y_pos_strand, id = z->y_id; ql = qe - qs;
    if(hpc_g) t_tot_l = hpc_len(*hpc_g, id);
    else if(uref) t_tot_l = uref->ug->u.a[id].len;
    else t_tot_l = Get_READ_LENGTH((*rref), id);

    if(mode == 3) {//semi
        ts = (qs - z->x_pos_s) + z->y_pos_s; ts += y_start_offset(qs, &(z->f_cigar));
        te = ts + ql; 
    } else if(mode == 1) {///forward extension
        te = ts + ql; 
    } else if(mode == 2) {///backward extension
        ts = te - ql; 
    }
    if(ts < 0) ts = 0;
    if(ts > t_tot_l) ts = t_tot_l;
    if(te > t_tot_l) te = t_tot_l;
    ql = qe - qs; tl = te - ts; 
    if(ql != tl) return 0;

    q_string = qstr + qs; 
    t_string = retrieve_str_seq_exz(tu, ts, tl, (*pts), (*pte)-(*pts), rev, uref, hpc_g, rref, id);
    (*pts) = ts; (*pte) = te;

    if(memcmp(q_string, t_string, ql)) return 0;
    exz->err = 0; push_trace(&(exz->cigar), 0, ql);
    exz->pl = tl; exz->ps = 0; exz->pe = tl-1;
    exz->tl = ql; exz->ts = 0; exz->te = ql-1;
    // cigar_check(t_string, q_string, exz);
    // if(!cigar_check(t_string, q_string, exz)) {
    //     fprintf(stderr, "[M::%s::] cigar_n::%d\n", __func__, (int32_t)exz->cigar.n);
    //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tl, t_string);
    //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)ql, q_string);
    // }
    exz->ps += ts; exz->pe += ts;
    exz->ts += qs; exz->te += qs;
    return 1;
}


//[qmin, qmax) && [tmin, tmax)






void append_wcigar(window_list *idx, window_list_alloc *res, bit_extz_t *exz)
{
    // fprintf(stderr, "[M::%s::] idx->cidx::%u, idx->clen::%u, res->c.n_0::%u, exz->cigar.n::%u, ", 
    //         __func__, idx->cidx, idx->clen, (uint32_t)res->c.n, (uint32_t)exz->cigar.n);
    if(idx->clen > 0) {
        // if(!(idx->cidx+idx->clen == res->c.n)) {
        //     fprintf(stderr, "[M::%s::] idx->cidx::%u, idx->clen::%u, res->c.n::%u\n", 
        //     __func__, idx->cidx, idx->clen, (uint32_t)res->c.n);
        // }
        assert(idx->cidx+idx->clen == res->c.n);
        if(exz->cigar.n > 0) {
            uint16_t c0, c; uint32_t l0, l, ci = 0, cn; 
            ///last item of old cigar
            c0 = (res->c.a[res->c.n-1]>>14); l0 = (res->c.a[res->c.n-1]&(0x3fff));
            ///first item of new cigar
            ci = pop_trace(&(exz->cigar), ci, &c, &l);
            if(c0 == c) {l += l0; res->c.n--; idx->clen--;}
            push_trace(((asg16_v *)(&(res->c))), c, l); 
            idx->clen = res->c.n-idx->cidx;

            cn = exz->cigar.n-ci;
            if(cn > 0) {
                kv_resize(uint16_t, res->c, (res->c.n+cn));
                memcpy(res->c.a+res->c.n, exz->cigar.a+ci, cn*sizeof(*(res->c.a)));
                idx->clen += cn; res->c.n += cn;
            }
        }
    } else {
        push_wcigar(idx, res, exz);///if exz is the first item
    }
    // fprintf(stderr, "[M::%s::] res->c.n::%u\n", __func__, (uint32_t)res->c.n);
}

void push_alnw(overlap_region *aux_o, bit_extz_t *exz)
{
    window_list *p = NULL; int64_t t;
    if(aux_o->w_list.n > 0) {
        p = &(aux_o->w_list.a[aux_o->w_list.n-1]);
        // fprintf(stderr, "+[M::%s::wn->%d] px::[%d, %d], py::[%d, %d], pe::%d, exz->t::[%u, %u], exz->p::[%u, %u], exz->e::%d, clen::%u\n", 
        //     __func__, (int32_t)(aux_o->w_list.n), p->x_start, p->x_end, p->y_start, p->y_end, p->error, 
        //     exz->ts, exz->te, exz->ps, exz->pe, exz->err, p->clen);
        // assert((p->x_end<exz->ts)&&(p->y_end<exz->ps));
        
        if(p->clen > 0) {
            t = ((int64_t)p->error) + ((int64_t)exz->err);
            ///note: t cannot be equal to INT16_MAX; otherwise it is unable to distiguish unaligned regions
            if(((p->x_end+1) == exz->ts) && ((p->y_end+1) == exz->ps) && (t < INT16_MAX)) {
                p->x_end = exz->te; p->y_end = exz->pe; p->error += exz->err;
                append_wcigar(p, &(aux_o->w_list), exz);
            //     fprintf(stderr, "-[M::%s::wn->%d] px::[%d, %d], py::[%d, %d], pe::%d, exz->t::[%u, %u], exz->p::[%u, %u], exz->e::%d, clen::%u\n", 
            // __func__, (int32_t)(aux_o->w_list.n), p->x_start, p->x_end, p->y_start, p->y_end, p->error, 
            // exz->ts, exz->te, exz->ps, exz->pe, exz->err, p->clen);
                return;
            }
        }
    }
    // fprintf(stderr, "[M::%s::wn->%d] exz->t::[%u, %u], exz->p::[%u, %u], exz->e::%d\n", 
    //         __func__, (int32_t)(aux_o->w_list.n), exz->ts, exz->te, exz->ps, exz->pe, exz->err);
    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = exz->ts; p->x_end = exz->te;
    p->y_start = exz->ps; p->y_end = exz->pe;
    p->error_threshold = 0; p->error = exz->err;///single round of alignment cannot have INT16_MAX errors
    push_wcigar(p, &(aux_o->w_list), exz);
}

///[qs, qe] && [ts, te]
void push_unmap_alnw(overlap_region *aux_o, int32_t qs, int64_t qe, int64_t ts, int64_t te, int64_t mode)
{
    window_list *p = NULL;
    kv_pushp(window_list, aux_o->w_list, &p);
    p->x_start = qs; p->x_end = qe;
    p->y_start = ts; p->y_end = te;
    p->error_threshold = mode; p->error = INT16_MAX;
    p->extra_begin = p->extra_end = -1;
    p->cidx = p->clen = 0;
}

///[qs, qe] && [ts, te]

///[qs, qe] && [ts, te]




void set_exact_exz(bit_extz_t *exz, int64_t qs, int64_t qe, int64_t ts, int64_t te)
{
    clear_align(*exz); exz->thre = 0; exz->cigar.n = 0;
    exz->err = 0; push_trace(&(exz->cigar), 0, qe - qs);
    exz->pl = te - ts; exz->ps = 0; exz->pe = exz->pl-1;
    exz->tl = qe - qs; exz->ts = 0; exz->te = exz->tl-1;
    exz->ps += ts; exz->pe += ts;
    exz->ts += qs; exz->te += qs;
}


int64_t hc_aln_exz_adv_hc(overlap_region *z, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, 
char* qstr, UC_Read *tu, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t mode, int64_t wl, 
bit_extz_t *exz, int64_t q_tot, double e_rate, int64_t maxl, int64_t maxe, int64_t force_l, 
int64_t estimate_err, overlap_region *aux_o)
{
    clear_align(*exz); exz->thre = 0; 
    if(((ts == -1) && (te == -1))) mode = 3;///set to semi-global
    int64_t thre, ql = qe - qs, thre0, pts = -1, pte = -1, pthre = -1, full = 0;
    if(ql == 0 && (te-ts) == 0) return 1;
    if((ql <= 0) || (te-ts) <= 0) return 0;
    if(estimate_err < 0) estimate_err = cal_estimate_err_hc(z, wl, qs, qe, ts, te, e_rate, &full);
    

    
    // if(ql <= 16) {
    if(estimate_err == 0) {
        if(full) {
            // if(!cal_exact_exz(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, q_tot, mode)) {
            //     fprintf(stderr, "[M::%s::ql::%ld::%c] xid::%d, yid::%d, qs::[%ld, %ld), ts::[%ld, %ld), mode::%ld, est_err::%ld, e_rate::%f, maxe::%ld\n", 
            //                                             __func__, ql, "+-"[z->y_pos_strand], z->x_id, z->y_id, qs, qe, ts, te, mode, estimate_err, e_rate, maxe);
            //     exit(1);
            // }
            set_exact_exz(exz, qs, qe, ts, te); push_alnw(aux_o, exz);
            return 1;
        } else if(cal_exact_exz(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, q_tot, mode)) {
            // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
            // fprintf(stderr, ", err::%d, thre::%d, scale::0(+)\n", exz->err, exz->thre);
            push_alnw(aux_o, exz);
            return 1;
        }
    }

    if(ql <= maxl && (estimate_err>>1) <= maxe) {
        thre = scale_ed_thre(estimate_err, maxe); 
        if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
            // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
            // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(+)\n", exz->err, exz->thre, thre);
            push_alnw(aux_o, exz);
            return 1;
        }

        thre0 = thre; thre = ql*e_rate; thre = scale_ed_thre(thre, maxe); 
        if(thre > thre0) {
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }

        thre0 = thre; thre <<= 1; thre = scale_ed_thre(thre, maxe); 
        if(thre > thre0) {
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(-)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }

        thre0 = thre; thre = ql*0.51; thre = scale_ed_thre(thre, maxe); 
        if(thre > thre0) {
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }

        if(ql <= force_l) {
            thre = maxe; 
            if(cal_exz_infi_adv(z, uref, hpc_g, rref, exz, qstr, tu, qs, qe, ts, te, &pts, &pte, thre, &pthre, q_tot, mode)) {
                // ref_cigar_check(qstr, tu, uref, hpc_g, rref, z->y_id, z->y_pos_strand, exz);
                // fprintf(stderr, ", err::%d, thre::%d, scale::%ld(*)\n", exz->err, exz->thre, thre);
                push_alnw(aux_o, exz);
                return 1;
            }
        }
    }
    // fprintf(stderr, ", err::%d, thre::%d\n", INT32_MAX, exz->thre);
    // if(mode == 0) {
    //     fprintf(stderr, "[M::%s::] pstr::%.*s\n", __func__, (int32_t)tu->length, tu->seq);
    //     fprintf(stderr, "[M::%s::] tstr::%.*s\n", __func__, (int32_t)(qe-qs), qstr+qs);
    // }
    return 0;
    
}









inline int64_t translate_double_mode(uint64_t double_mode, uint64_t is_backward)
{
    if(double_mode == 0) return 0;
    if(double_mode == 4) return 3;
    if(double_mode == 1 || double_mode == 2) return double_mode;
}


// void chain_win_aln(overlap_region *z, Chain_Data *dp, Candidates_list *cl, int64_t qs, int64_t qe, 
// int64_t ts, int64_t te, int64_t ql, int64_t tl, int64_t wl, int64_t mode, bit_extz_t *exz)
// {

// }



inline void push_khit(Candidates_list *res, int32_t xs, int32_t ys, uint32_t len, uint32_t h_khit, uint32_t *ic)
{
    uint32_t p, c; k_mer_hit *z;
    c = ((len >= h_khit)?1:2); if(ic) c = *ic; c <<= 8;
    if(len > 0) {
        while (len >= (0xffu)) {
            p = (c + (0xffu)); 
            kv_pushp_cl(k_mer_hit, (*res), &z); 
            memset(z, 0, sizeof((*z)));
            z->self_offset = xs; z->offset = ys; z->cnt = p;
            len -= (0xffu);
        }
        if(len) {
            p = (c + len); 
            kv_pushp_cl(k_mer_hit, (*res), &z); 
            memset(z, 0, sizeof((*z)));
            z->self_offset = xs; z->offset = ys; z->cnt = p;
        }
    } else {
        p = (c + len); 
        kv_pushp_cl(k_mer_hit, (*res), &z); 
        memset(z, 0, sizeof((*z)));
        z->self_offset = xs; z->offset = ys; z->cnt = p;
    }
}

uint32_t extract_exact_cigar(asg16_v *ez, int32_t ps, int32_t ts, int32_t pmin, int32_t pmax, 
int32_t tmin, int32_t tmax, Candidates_list *res, int32_t minl, int64_t min_w_l, int64_t h_khit)
{
    uint32_t ci = 0, cl, occ = 0; uint16_t c; 
    int32_t pi = ps, ti = ts, p[2], t[2], poff, toff, maxl;
    int32_t pos, poe, tos, toe, l; poff = toff = maxl = -1;
    while (ci < ez->n && pi < pmax && ti < tmax) {
        ci = pop_trace(ez, ci, &c, &cl);
        if(c == 0) {
            p[0] = pi; p[1] = pi + cl;
            t[0] = ti; t[1] = ti + cl;
            pos = MAX(p[0], pmin); poe = MIN(p[1], pmax);
            tos = MAX(t[0], tmin); toe = MIN(t[1], tmax);
            if((poe > pos) && (toe > tos)) {
                l = poe - pos;
                if(l == (toe - tos)) {
                    poe--; toe--;
                    if(l > maxl) {
                        poff = poe; toff = toe; maxl = l;
                    }
                    if(l >= minl) {
                        push_khit(res, toe, poe, l, h_khit, NULL); occ++;
                    }
                }
            }
            pi+=cl; ti+=cl;
        } else if(c == 1) {
            pi+=cl; ti+=cl;
        } else if(c == 2) {///more p
            pi+=cl;
        } else if(c == 3) {
            ti+=cl;
        }
    }
    ///(ts >= tmin) && (ti >= (min_w_l + ts)): here is a whole window
    if(maxl > 0 && maxl < minl && (ts >= tmin) && (ti >= (min_w_l + ts))) {
        uint32_t w = 3;
        push_khit(res, toff, poff, maxl, h_khit, &w); occ++;
    }
    return occ;
}


int64_t gen_single_khit(Candidates_list *cl, int64_t ch_n, int64_t h_khit, int64_t mode, int64_t qs, int64_t qe, int64_t ts, int64_t te, int64_t max_skip, int64_t max_iter, int64_t rid)
{
    // if(ch_n != 3 || mode != 0 || qs != 171728 || qe != 172258) return 0;
    k_mer_hit *ch_a = cl->list + cl->length; int64_t k, i, j, occ, m, ncn, prefix, suffix, srt = 1;
    prefix = suffix = 0;
    if(mode == 0 || mode == 2) suffix = 1;
    if(mode == 0 || mode == 1) prefix = 1;
    // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "[M::%s::mode->%ld] ch_n::%ld, q::[%ld, %ld), t::[%ld, %ld)\n", 
        //     __func__, mode, ch_n, qs, qe, ts, te);
    // }
    
    for (k = occ = m = 0; k < ch_n; k++) {
        // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "+i::%ld[M::%s::] x::[%u, %u), y::[%u, %u), cnt::%u\n", k, __func__, 
        //     ch_a[k].self_offset+1-(ch_a[k].cnt&((uint32_t)(0xffu))), ch_a[k].self_offset+1,
        //     ch_a[k].offset+1-(ch_a[k].cnt&((uint32_t)(0xffu))), ch_a[k].offset+1, (ch_a[k].cnt&(0xffu)));
        // }
        if(!(ch_a[k].cnt&(0xffu))) continue;
        occ++;
        if((ch_a[k].cnt&(0xffu)) > 1) occ++;
        ch_a[m++] = ch_a[k];
    }
    ch_n = m; if(!ch_n) return ch_n;
    occ += prefix + suffix;
    // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "+[M::%s::] occ::%ld\n", __func__, occ);
    // }

    ncn = occ + cl->length;
    if(cl->size < ncn) {
        cl->size = ncn;
        REALLOC(cl->list, cl->size);
        // cl->list = (k_mer_hit*)realloc(cl->list, (sizeof((*(cl->list)))*cl->length));
    }
    ch_a = cl->list + cl->length; assert((cl->length+occ)<= cl->size);

    k_mer_hit cht;
    ///global or backward
    if(suffix) {
        cht.self_offset = qe;
        cht.offset = te;
        cht.cnt = 1; cht.readID = 1;//make it as primary chain
        cht.strand = 0;
        ch_a[--occ] = cht;
        // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
        //     cht.self_offset, cht.offset, cht.cnt, cht.readID);
    }
    for (k = ch_n-1; k >= 0; k--) {
        if(!(ch_a[k].cnt&(0xffu))) continue;
        ///end
        cht.self_offset = ch_a[k].self_offset+1;
        cht.offset = ch_a[k].offset+1;
        cht.strand = 0;
        cht.cnt = cht.readID = (ch_a[k].cnt&(0xffu));
        //make it as non-primary chain
        if((ch_a[k].cnt&(0xffu)) < h_khit) cht.readID = cht.cnt + 1;
        ch_a[--occ] = cht;
        // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
        //     cht.self_offset, cht.offset, cht.cnt, cht.readID);

        if((ch_a[k].cnt&(0xffu)) > 1) {
            ///start
            cht.self_offset = ch_a[k].self_offset+1-(ch_a[k].cnt&(0xffu));
            cht.offset = ch_a[k].offset+1-(ch_a[k].cnt&(0xffu));
            cht.strand = 0;
            cht.cnt = cht.readID = (ch_a[k].cnt&(0xffu));
            //make it as non-primary chain
            if((ch_a[k].cnt&(0xffu)) < h_khit) cht.readID = cht.cnt + 1;
            ch_a[--occ] = cht;
            // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
            // cht.self_offset, cht.offset, cht.cnt, cht.readID);
        }
    }
    
    if(prefix) { ///global or forward
        cht.self_offset = qs;
        cht.offset = ts;
        cht.cnt = 1; cht.readID = 1;//make it as primary chain
        cht.strand = 0;
        ch_a[--occ] = cht;
        // fprintf(stderr, "occ::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", occ, __func__, 
        //     cht.self_offset, cht.offset, cht.cnt, cht.readID);
    }
    // if(ch_n == 2 && mode == 2 && qe - qs == 2419 && te - ts == 2419) {
        // fprintf(stderr, "-[M::%s::] occ::%ld\n", __func__, occ);
    // }
    // if(!(occ == 0)) {
    //     fprintf(stderr, "[M::%s] rid::%ld, name::%.*s\n", __func__, rid,
    //      (int32_t)UL_INF.nid.a[rid].n, UL_INF.nid.a[rid].a);
    // }
    assert(occ == 0);
    ch_n = occ = ncn - cl->length;
    uint64_t q[2], t[2]; 
    q[0] = q[1] = t[0] = t[1] = (uint64_t)-1;
    if(prefix) {
        q[0] = qs; t[0] = ts;
    }
    if(suffix) {
        q[1] = qe; t[1] = te;
    }

    // for (k = 0; k < ch_n; k++) {
    //     fprintf(stderr, "<i::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", k, __func__, 
    //         ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    // }
    for (k = m = occ = 0; k < ch_n; k++) {
        if((k>0) && (ch_a[k].self_offset==q[0]) && (ch_a[k].offset=t[0])) continue;
        if(((k+1)<ch_n) && (ch_a[k].self_offset==q[1]) && (ch_a[k].offset=t[1])) continue;
        if(m>0) {
            if((ch_a[k].self_offset>ch_a[m-1].self_offset) && (ch_a[k].offset>ch_a[m-1].offset)) {
                occ++;
            }
            if(ch_a[k].self_offset<=ch_a[m-1].self_offset) srt = 0;
        } else {
            occ++;
        }
        ch_a[m++] = ch_a[k];
    }
    ch_n = m;
    if(occ == ch_n) return ch_n;///already colinear
    if(!srt) {
        radix_sort_k_mer_hit_self(ch_a, ch_a + ch_n);
        for (i = 1, j = 0; i <= ch_n; i++) {
            if (i == ch_n || ch_a[i].self_offset != ch_a[j].self_offset) {
                if(i - j > 1) radix_sort_k_mer_hit_off(ch_a+j, ch_a+i);
                j = i;
            }
        }
    }

    // for (k = 0; k < ch_n; k++) {
    //     fprintf(stderr, ">i::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", k, __func__, 
    //         ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    // }
    occ = ch_n;
    ch_n = lchain_simple(ch_a+prefix, ch_n-prefix-suffix, ch_a+prefix, &(cl->chainDP), max_skip, max_iter);
    ch_n += prefix + suffix; if(suffix) ch_a[ch_n-1] = ch_a[occ-1];
    // for (k = 0; k < ch_n; k++) {
    //     fprintf(stderr, "-i::%ld[M::%s::] x::%u, y::%u, cnt::%u, cov::%u\n", k, __func__, 
    //         ch_a[k].self_offset, ch_a[k].offset, ch_a[k].cnt, ch_a[k].readID);
    // }
    
    return ch_n;
}

///[qs, qe) && [ts, te)
int64_t gen_win_chain(overlap_region *z, Candidates_list *cl, int64_t qs, int64_t qe, int64_t ts, int64_t te, 
int64_t wl, const ul_idx_t *uref, hpc_t *hpc_g, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, 
int64_t ql, int64_t tl, double e_rate, int64_t h_khit, int64_t mode, int64_t rid, int64_t is_accurate)
{
    assert(mode < 3);
    int64_t k, ws, we, os, oe, wsk, rcn = cl->length, ncn, occ = 0, ovlp, wn = z->w_list.n; asg16_v ez; uint32_t w = 1;
    ws = qs; if(ws < z->x_pos_s) ws = z->x_pos_s; 
    we = qe-1; if(we > z->x_pos_e) we = z->x_pos_e; 
    wsk = get_win_id_by_s(z, ((ws/wl)*wl), wl, NULL); 
    // if(rid == 7) {
    //     fprintf(stderr, "[M::%s::]\tutg%.6u%c\t%u\t%u\t%c\tutg%.6u%c\t%u\t%u\tq::[%ld,%ld)\tw::[%ld,%ld]\twn::%ld\twsk::%ld\n", __func__, 
	// 			z->x_id+1, "lc"[uref->ug->u.a[z->x_id].circ], 
	// 			z->x_pos_s, z->x_pos_e+1, "+-"[z->y_pos_strand],
	// 			z->y_id+1, "lc"[uref->ug->u.a[z->y_id].circ], 
	// 			z->y_pos_s, z->y_pos_e+1, qs, qe, ws, we, wn, wsk);
    // }
    // assert((ws>=z->w_list.a[wsk].x_start) && (ws<=z->w_list.a[wsk].x_end));
    // wek = get_win_id_by_e(z, ((we/wl)*wl), wl, NULL);
    // assert((we>=z->w_list.a[wek].x_start) && (we<=z->w_list.a[wek].x_end));
    for(wsk=((wsk<wn)?(wsk):(wn-1)); wsk<wn && qs>z->w_list.a[wsk].x_end; wsk++);
    for(wsk=((wsk<wn)?(wsk):(wn-1)); wsk>=0 && qs<z->w_list.a[wsk].x_start; wsk--);
    if(wsk < 0) wsk = 0; ///qs >= z->w_list.a[wsk].x_start && qs <= z->w_list.a[wsk].x_end  
    ///global or forward
    if(mode == 0 || mode == 1) push_khit(cl, qs, ts, 0, 0, &w);
    //[ws, we] && [wsk, wek]; [qs, qe) && [ts, te)
    for (k = wsk; k<wn && z->w_list.a[k].x_start<qe; k++) {
        if(z->w_list.a[k].y_end == -1) continue;
        ws = z->w_list.a[k].x_start; we = z->w_list.a[k].x_end+1;
        os = MAX(qs, ws); oe = MIN(qe, we); ovlp = ((oe>os)? (oe-os):0);
        if(!ovlp) continue;
        if(!(z->w_list.a[k].clen)) {
            gen_backtrace_adv_exz(&(z->w_list.a[k]), z, rref, hpc_g, uref, qstr, tu->seq, exz, z->y_pos_strand, z->y_id);
        }
        ez.a = z->w_list.c.a + z->w_list.a[k].cidx; 
        ez.n = ez.m = z->w_list.a[k].clen; 
        occ += extract_exact_cigar(&ez, z->w_list.a[k].y_start, z->w_list.a[k].x_start, ts, te, qs, qe, cl, 10, wl, h_khit);
    }
    ///global or backward
    if(mode == 0 || mode == 2) push_khit(cl, qe-1, te-1, 0, 0, &w);
    // fprintf(stderr, "[M::%s::] rcn::%ld, cl->length::%lld\n", __func__, rcn, cl->length);
    if(!occ) {
        cl->length = rcn; return 0;
    }
    ncn = cl->length; cl->length = rcn;
    k_mer_hit *ch_a = cl->list + rcn; int64_t ch_n0 = ncn - rcn, ch_n;
    int64_t max_skip, max_iter, max_dis, quick_check; double chn_pen_gap, chn_pen_skip;
	set_lchain_dp_op(is_accurate, h_khit, &max_skip, &max_iter, &max_dis, &chn_pen_gap, &chn_pen_skip, &quick_check);
    max_dis = MAX_SIN_L>>1;
    // for (k = 0; k < ch_n0; k++) {
    //     assert(debug_k_mer_hit_retrive(&(ch_a[k]), hpc_g, rref, uref, qstr, tu, z->y_id, z->y_pos_strand));
    // }
    ch_n = lchain_qdp_fix(ch_a, ch_n0, &(cl->chainDP), max_skip, max_iter, max_dis, chn_pen_gap, chn_pen_skip, 
                e_rate, ql, tl, 1, ((mode==0)||(mode==1))?1:0,  ((mode==0)||(mode==2))?1:0);
    for (k = occ = 0; k < ch_n; k++) {
        ch_a[k] = ch_a[cl->chainDP.tmp[k]];
        if((ch_a[k].cnt&(0xffu))) occ++;
        // assert(debug_k_mer_hit_retrive(&(ch_a[k]), hpc_g, rref, uref, qstr, tu, z->y_id, z->y_pos_strand));
    }
    // fprintf(stderr, "[M::%s::] ch_n0::%ld, ch_n::%ld, mode::%ld, ql::%ld, tl::%ld, occ::%ld\n", 
    // __func__, ch_n0, ch_n, mode, qe-qs, te-ts, occ);
    if(occ <= 0) return 0;
    ch_n = gen_single_khit(cl, ch_n, h_khit, mode, qs, qe, ts, te, max_skip, max_iter, rid);
    return ch_n;
}


///[qs, qe) && [ts, te)





void update_overlap_region(overlap_region *des, overlap_region *src, int64_t xl, int64_t yl)
{
    kv_resize(uint16_t, des->w_list.c, src->w_list.c.n);
    des->w_list.c.n = src->w_list.c.n;
    memcpy(des->w_list.c.a, src->w_list.c.a, src->w_list.c.n*(sizeof((*(src->w_list.c.a)))));

    kv_resize(window_list, des->w_list, src->w_list.n);
    des->w_list.n = src->w_list.n;
    memcpy(des->w_list.a, src->w_list.a, src->w_list.n*(sizeof((*(src->w_list.a)))));

    if(src->w_list.n) {
        des->x_pos_s = src->w_list.a[0].x_start; des->x_pos_e = src->w_list.a[src->w_list.n-1].x_end; 
        des->y_pos_s = src->w_list.a[0].y_start; des->y_pos_e = src->w_list.a[src->w_list.n-1].y_end;
    }

    int64_t xr, yr; 
    if(des->x_pos_s <= des->y_pos_s) {
        des->y_pos_s -= des->x_pos_s; des->x_pos_s = 0;
    } else {
        des->x_pos_s -= des->y_pos_s; des->y_pos_s = 0;
    }

    xr = xl-des->x_pos_e-1; yr = yl-des->y_pos_e-1;
    if(xr <= yr) {
        des->x_pos_e = xl-1; des->y_pos_e += xr;        
    } else {
        des->y_pos_e = yl-1; des->x_pos_e += yr; 
    }
}




void hc_ovlp_base_direct(overlap_region *z, k_mer_hit *ch_a, int64_t ch_n, int64_t wl, All_reads *rref, char* qstr, UC_Read *tu, 
bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, int64_t tl, uint64_t rid, int64_t pre_mode)
{
    int64_t i, l, mode, q[2], t[2], qr, tr, is_done, zn, si, ei;

    if((pre_mode < 0) && (z->non_homopolymer_errors == 0) && (z->w_list.n)) {
        zn = z->w_list.n;
        for (i = 1; i < zn; i++) {
            if((z->w_list.a[i].error == 0 && z->w_list.a[i-1].error == 0) && (z->w_list.a[i].x_start == z->w_list.a[i-1].x_end + 1) && 
                (z->w_list.a[i].y_end == (z->w_list.a[i-1].y_end + (z->w_list.a[i].x_end-z->w_list.a[i-1].x_end)))) {
                    continue;
            }
            break;
        }
        if(i >= zn) {
            q[0] = z->w_list.a[0].x_start; q[1] = z->w_list.a[z->w_list.n-1].x_end;  
            t[1] = z->w_list.a[z->w_list.n-1].y_end; t[0] = z->w_list.a[0].y_end - (z->w_list.a[0].x_end-z->w_list.a[0].x_start); 

            if(q[0] <= t[0]) {
                t[0] -= q[0]; q[0] = 0;
            } else {
                q[0] -= t[0]; t[0] = 0;
            }

            qr = ql-q[1]-1; tr = tl-t[1]-1;
            if(qr <= tr) {
                q[1] = ql-1; t[1] += qr;        
            } else {
                t[1] = tl-1; q[1] += tr; 
            }

            if(q[0] == z->w_list.a[0].x_start && q[1] == z->w_list.a[z->w_list.n-1].x_end) {
                // fprintf(stderr, "[M::%s::%u->%u::%c] ovlp::%u, w_list.n::%u\n", __func__, z->x_id, z->y_id+1, "+-"[z->y_pos_strand], z->x_pos_e+1-z->x_pos_s, (uint32_t)z->w_list.n);
                set_exact_exz(exz, q[0], q[1] + 1, t[0], t[1] + 1); push_alnw(aux_o, exz);
                return;
            }
        }
    }

    si = 0; ei = ch_n;
    if(pre_mode == 0) {
        si = 1; ei = ch_n - 1;
    } else if(pre_mode == 1) {
        si = 1;
    } else if(pre_mode == 2) {
        ei = ch_n - 1;
    }

    for (l = si - 1, i = si; i <= ei; i++) {
        q[0] = q[1] = t[0] = t[1] = mode = -1; is_done = 0;
        if(l >= 0) {
            q[0] = ch_a[l].self_offset; t[0] = ch_a[l].offset;
        } else {
            q[0] = 0;
        }

        if(i < ch_n) {
            q[1] = ch_a[i].self_offset; t[1] = ch_a[i].offset;
        } else {
            q[1] = ql;
        }

        if((t[0] != -1) && (t[1] != -1)) {
            mode = 0;//global
        } else if((t[0] != -1) && (t[1] == -1)) {
            mode = 1;///forward extension
        } else if((t[0] == -1) && (t[1] != -1)) {
            mode = 2;///backward extension
        } else {
            mode = 3;///no primary hit within [ibeg, iend] 
        }

        // if(z->x_id == 57 && z->y_id == 2175) {
        if(mode == 1 || mode == 2) adjust_ext_offset(&(q[0]), &(q[1]), &(t[0]), &(t[1]), ql, tl, 0, mode);
        //     fprintf(stderr, "#[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld\n", 
        //             __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode);
        // }
        is_done = hc_aln_exz_adv_hc(z, NULL, NULL, rref, qstr, tu, q[0], q[1], t[0], t[1], mode, wl, exz, ql, e_rate, 
        MAX_SIN_L, MAX_SIN_E, FORCE_SIN_L, -1, aux_o);

        // if(z->x_id == 57 && z->y_id == 2175) {
            //     fprintf(stderr, "-is_done::%ld[M::%s::] utg%.6dl(%c), q::[%ld, %ld), t::[%ld, %ld), mode::%ld, ch_n::%ld\n", 
            //         is_done, __func__, (int32_t)z->y_id+1, "+-"[z->y_pos_strand], q[0], q[1], t[0], t[1], mode, ch_n);
        // }
                
        if(!is_done) {///postprocess
            push_unmap_alnw(aux_o, q[0], q[1]-1, t[0], t[1]-1, mode);
        }
        l = i;
    }
}



void rechain_aln_hc(overlap_region *z, Candidates_list *cl, overlap_region *aux_o, int64_t aux_i, int64_t wl, 
All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, double e_rate, int64_t ql, int64_t tl, int64_t h_khit, int64_t rid)
{
    int64_t rcn = cl->length, ch_n, qs, qe, ts, te, mode, an0, an, todo; 
    k_mer_hit *ch_a; uint8_t q[2], t[2]; ///ul_ov_t idx; 
    ///[qs, qe) && [ts, te)
    qs = aux_o->w_list.a[aux_i].x_start; qe = aux_o->w_list.a[aux_i].x_end+1;
    ts = aux_o->w_list.a[aux_i].y_start; te = aux_o->w_list.a[aux_i].y_end+1;
    if(qe - qs < FORCE_SIN_L || te - ts < FORCE_SIN_L) return;
    mode = aux_o->w_list.a[aux_i].error_threshold;
    ch_n = gen_win_chain(z, cl, qs, qe, ts, te, wl, NULL, NULL, rref, qstr, tu, exz, ql, tl, e_rate, h_khit, mode, rid, 1);
    ch_a = cl->list + rcn;
    if(ch_n) {
        todo = 1; ///idx.ts = idx.te = (uint32_t)-1; idx.qs = 0; idx.qe = ql; 
        if(mode == 0) {//global
            // idx.qn = 0; idx.tn = ch_n - 1;
            // idx.qs = ch_a[idx.qn].self_offset;
            // idx.ts = ch_a[idx.qn].offset;
            // idx.qe = ch_a[idx.tn].self_offset;
            // idx.te = ch_a[idx.tn].offset;
            assert(ch_a[0].self_offset == qs && ch_a[0].offset == ts);
            assert(ch_a[ch_n-1].self_offset == qe && ch_a[ch_n-1].offset == te);
            if(ch_n <= 2) todo = 0;
        } else if(mode == 1) {//forward ext
            // idx.qn = 0; idx.tn = ch_n;
            // idx.qs = ch_a[idx.qn].self_offset;
            // idx.ts = ch_a[idx.qn].offset;
            // idx.qe = ql; 
            assert(ch_a[0].self_offset == qs && ch_a[0].offset == ts);
            if(ch_n <= 1) todo = 0;
        } else if(mode == 2) {///backward ext
            // idx.qn = (uint32_t)-1; idx.tn = ch_n-1;
            // idx.qs = 0; 
            // idx.qe = ch_a[idx.tn].self_offset;
            // idx.te = ch_a[idx.tn].offset;
            assert(ch_a[ch_n-1].self_offset == qe && ch_a[ch_n-1].offset == te);
            if(ch_n <= 1) todo = 0;
        }
        if(todo) {
            an0 = aux_o->w_list.n;
            // if(z->x_id == 29033 && z->y_id == 21307) {
            //     fprintf(stderr, "[M::%s]\tan0::%ld\tq::[%u,\t%u)\tt::[%u,\t%u)\tlw::%u\trw::%u\n", __func__, an0,
            //     idx.qs, idx.qe, idx.ts, idx.te, idx.qn, idx.tn);
            // }
            // ovlp_base_aln(z, ch_a, ch_n, &idx, wl, uref, hpc_g, rref, qstr, tu, exz, aux_o, e_rate, ql, tl, (uint64_t)-1);
            hc_ovlp_base_direct(z, ch_a, ch_n, wl, rref, qstr, tu, exz, aux_o, e_rate, ql, tl, (uint64_t)-1, mode);
            an = aux_o->w_list.n; q[0] = q[1] = t[0] = t[1] = 0; todo = 0;
            // if(z->x_id == 29033 && z->y_id == 21307) {
            //     fprintf(stderr, "[M::%s]\tan::%ld\n", __func__, an);
            // }
            // fprintf(stderr, "[M::%s::] awn0::%ld, awn::%lu\n", __func__, an0, an);
            ///old unaligned window could be replaced by the new aligned window
            if((an == (an0 + 1)) && (!(is_ualn_win(aux_o->w_list.a[an-1])))) {
                if(aux_o->w_list.a[aux_i].x_start == aux_o->w_list.a[an-1].x_start) q[0] = 1;
                if(aux_o->w_list.a[aux_i].x_end == aux_o->w_list.a[an-1].x_end) q[1] = 1;
                if(aux_o->w_list.a[aux_i].y_start == aux_o->w_list.a[an-1].y_start) t[0] = 1;
                if(aux_o->w_list.a[aux_i].y_end == aux_o->w_list.a[an-1].y_end) t[1] = 1;
                if((mode == 0) && q[0] && q[1] && t[0] && t[1]) todo = 1;
                if((mode == 1) && q[0] && t[0]) todo = 1;
                if((mode == 2) && q[1] && t[1]) todo = 1;
                if(todo) {
                    aux_o->w_list.a[aux_i] = aux_o->w_list.a[an-1]; aux_o->w_list.n--;
                }
            }
            // if(an > an0) {///should always > 0 as there are unmapped windows
            // }
            // aux_o->w_list.n = an0;
        }
    }
    cl->length = rcn;///must reset!!!!
}


uint64_t gen_hc_fast_cigar0(overlap_region *z, Candidates_list *cl, uint64_t wl, All_reads *rref, char* qstr, UC_Read *tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, uint64_t rid, int64_t h_khit, int64_t *re)
{
    int64_t ch_idx = z->shared_seed, ch_n;
    int64_t i, tl, id = z->y_id, m, tot_e, aln, xe, ye; 
    k_mer_hit *ch_a = cl->list + ch_idx;
    tl = Get_READ_LENGTH((*rref), id);
    for (i = ch_idx; i < cl->length && cl->list[i].readID == cl->list[ch_idx].readID; i++); ch_n = i-ch_idx;
    if(ch_n <= 0) return 0;

    ///debug for memory
    // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);
    
    // fprintf(stderr, "[M::%s::rid->%ld] utg%.6dl(%c), z::[%u, %u)\n", 
    // __func__, rid, (int32_t)z->y_id+1, "+-"[z->y_pos_strand],  z->x_pos_s, z->x_pos_e+1);
    aux_o->w_list.n = aux_o->w_list.c.n = 0; 
    aux_o->y_id = z->y_id; aux_o->y_pos_strand = z->y_pos_strand;
    aux_o->x_pos_s = z->x_pos_s; aux_o->x_pos_e = z->x_pos_e;
    aux_o->y_pos_s = z->y_pos_s; aux_o->y_pos_e = z->y_pos_e;

    hc_ovlp_base_direct(z, ch_a, ch_n, wl, rref, qstr, tu, exz, aux_o, e_rate, ql, tl, rid, -1);

    int64_t aux_n = aux_o->w_list.n;
    for (i = 0; i < aux_n; i++) {
        // if(z->y_id == 30129) {
        //     fprintf(stderr, "[aln::-i->%ld::ql->%d] q::[%d, %d), t::[%d, %d), err::%d, clen::%u, mode::%d\n", i, 
        //         aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start,
        //         aux_o->w_list.a[i].x_start, aux_o->w_list.a[i].x_end+1, 
        //         aux_o->w_list.a[i].y_start, aux_o->w_list.a[i].y_end+1, 
        //         aux_o->w_list.a[i].error, aux_o->w_list.a[i].clen, aux_o->w_list.a[i].error_threshold);
        // }
        if(!(is_ualn_win(aux_o->w_list.a[i]))) continue;
        // if((aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start) <= FORCE_CNS_L) {
            // fprintf(stderr, "[aln::-i->%ld::ql->%d] q::[%d, %d), t::[%d, %d), err::%d, clen::%u, mode::%d\n", i, 
            //             aux_o->w_list.a[i].x_end+1-aux_o->w_list.a[i].x_start,
            //             aux_o->w_list.a[i].x_start, aux_o->w_list.a[i].x_end+1, 
            //             aux_o->w_list.a[i].y_start, aux_o->w_list.a[i].y_end+1, 
            //             aux_o->w_list.a[i].error, aux_o->w_list.a[i].clen, aux_o->w_list.a[i].error_threshold);
        // }
        //will overwrite ch_a; does not matter
        rechain_aln_hc(z, cl, aux_o, i, wl, rref, qstr, tu, exz, e_rate, ql, tl, h_khit, rid);
    }

    if(((int64_t)aux_o->w_list.n) > aux_n) {
        for (i = m = 0; i < ((int64_t)aux_o->w_list.n); i++) {
            if((i < aux_n) && (is_ualn_win(aux_o->w_list.a[i]))) continue;
            aux_o->w_list.a[m++] = aux_o->w_list.a[i];
        }
        aux_o->w_list.n = m;
        radix_sort_window_list_xs_srt(aux_o->w_list.a, aux_o->w_list.a+aux_o->w_list.n);
    }

    ///update z by aux_o
    update_overlap_region(z, aux_o, ql, tl);

    aux_n = z->w_list.n;
    for (i = tot_e = aln = 0; i < aux_n; i++) {
        if(is_ualn_win(z->w_list.a[i])) {
            xe = z->w_list.a[i].x_end + 1 - z->w_list.a[i].x_start; 
            ye = z->w_list.a[i].y_end + 1 - z->w_list.a[i].y_start;
            tot_e += ((xe >= ye)?(xe):(ye));
        } else {
            tot_e += z->w_list.a[i].error; aln += z->w_list.a[i].x_end + 1 - z->w_list.a[i].x_start;
        }
    }
    *re = tot_e;
    // fprintf(stderr, "[M::%s::%u->%u::%c] ovlp::%u, aln::%ld, tot_e::%ld, w_list.n::%u, ch_n::%ld\n", 
    // __func__, z->x_id, z->y_id+1, "+-"[z->y_pos_strand], z->x_pos_e+1-z->x_pos_s, aln, tot_e, (uint32_t)z->w_list.n, ch_n);

    // debug_overlap_region(aux_o, qstr, tu, NULL, NULL, rref);


    // ch_a = cl->list + ch_idx; //update
    // for (i = 0; i < wn; i++) z->w_list.a[i].clen = 0;///clean cigar
    // if(on > 1) {
    //     fprintf(stderr, "[M::%s::] rid::%lu, on::%ld\n", __func__, rid, on);
    // }
    // if(z->y_id == 126) prt_k_mer_hit(ch_a, ch_n);
    // for (i = ch_i = 0; i < on; i++) {
    //     assert((i<=0)||(ov[i].qs > ov[i-1].qe));
    //     ov[i].sec = 16;///do not know the aln type
    //     ch_i = sub_base_aln(z, dp, ch_a, ch_n, pe, ov[i].qs, ov[i].qe, wl, uref, hpc_g, rref, qstr, tu, exz, e_rate, ql, tl, ch_i, rid);
    //     pe = ov[i].qe;
    // }
    
   return 1;
}




#define gen_err_unaligned(xl, yl) (((xl)<=FORCE_SIN_L)?(MAX((xl), (yl))):MAX((MIN((xl), (yl))), ((xl*0.51)+1)))


///[s, e)

///[s, e)

#define bst_ov(x) ((x).misBase)
#define ov_dif(x) ((x).cov)
#define ov_id(x) ((x).overlapID)
#define ov_xoff(x) ((x).site)
#define var_id(x) ((x).overlapSite)

#define var_s(x) ((x).site)
#define var_l(x) ((x).overlap_num)
#define var_occ(x) ((x).occ_0)
#define var_min_dif(x) ((x).score)
#define var_min_ovid(x) ((x).id)
#define var_h_idx(x) ((x).occ_1)





// #define id_mm ((uint64_t)0x7fffffffffffffff)
#define id_set ((uint64_t)0x8000000000000000)
#define id_get(a) ((uint32_t)(a))
#define err_get(a) (((a)&((uint64_t)0x7fffffffffffffff))>>32)

///[s, e)


///[s, e)

///[s, e)

///[s, e)
inline int64_t detect_near_cc_tlen(bit_extz_t *ez, int64_t ck0, int64_t xk0, int64_t yk0, uint8_t rev)
{
    int64_t ck = ck0, xk = xk0, yk = yk0, cn = ez->cigar.n, op;
    if(!rev) {
        if(ck >= cn) return yk;
        while (ck < cn) {
            op = ez->cigar.a[ck]>>14;
            if(op == 0) return yk;
            if(op!=2) xk += (ez->cigar.a[ck]&(0x3fff));
            if(op!=3) yk += (ez->cigar.a[ck]&(0x3fff)); 
            ck++; 
        }
    } else {
        if(ck <= 0) return yk;
        while (ck > 0) {
            --ck;
            op = ez->cigar.a[ck]>>14;
            if(op == 0) return yk;
            if(op!=2) xk -= (ez->cigar.a[ck]&(0x3fff));
            if(op!=3) yk -= (ez->cigar.a[ck]&(0x3fff)); 
        }
    }

    return yk;
}

///[s, e)

///[s, e)





///[s, e)





















inline uint64_t set_cgid(ul_ov_t *z, uint64_t *ia, uint64_t *iak, uint64_t ian, uint64_t *ca, uint64_t ci, overlap_region *oa, int64_t *gni)
{
    assert(ca[z->tn] <= 1);
    if(z->ts != ((uint32_t)-1)) return 0;
    if(ca[z->tn] != 0) {
        // fprintf(stderr, "+[M::%s] sec::%u\trid::%u\t%.*s\tq::[%u,%u)\terr::%u\tca[z->tn]::%lu\n", __func__, z->sec, oa[z->tn].y_id, (int)Get_NAME_LENGTH(R_INF, oa[z->tn].y_id), Get_NAME(R_INF, oa[z->tn].y_id), z->qs, z->qe, z->qn, ca[z->tn]);
        return 0;
    }

    z->ts = ci; ca[z->tn]++; 
    ia[(*iak)++] = z->tn;
    assert((*iak) <= ian);
    (*gni)++;
    // fprintf(stderr, "-[M::%s] sec::%u\trid::%u\t%.*s\tq::[%u,%u)\terr::%u\tca[z->tn]::%lu\n", __func__, z->sec, oa[z->tn].y_id, (int)Get_NAME_LENGTH(R_INF, oa[z->tn].y_id), Get_NAME(R_INF, oa[z->tn].y_id), z->qs, z->qe, z->qn, ca[z->tn]);
    
    return 1;
}

///need to print some examples for double check











#define is_ul_ov_pe(z, zn, i, mm) ((((i)+1)<(zn))&&((z)[(i)+1].qn==(z)[(i)].qn)&&((z)[(i)+1].tn==(z)[(i)].tn)&&((z)[(i)+1].rev==(z)[(i)].rev)\
            &&((z)[(i)+1].sec==(z)[(i)].sec+1)&&(((z)[(i)].sec)<(mm).a[(z)[(i)].qn].n)&&(((z)[(i)+1].sec)<(mm).a[(z)[(i)+1].qn].n)\
            &&((mm).a[(z)[(i)].qn].a[((z)[(i)].sec)].pe)&&((mm).a[(z)[(i)+1].qn].a[((z)[(i)+1].sec)].pe)\
            &&((mm).a[(z)[(i)].qn].a[((z)[(i)].sec)].dir==0)&&((mm).a[(z)[(i)+1].qn].a[((z)[(i)+1].sec)].dir==1))

#define is_exact_ov(z, mm)  (((((z).sec)<(mm).a[(z).qn].n))&&(!((mm).a[(z).qn].a[((z).sec)].pe))\
            &&((mm).a[(z).qn].a[((z).sec)].el==(uint32_t)-1))


///ref_n <= 2



#define aln2ov(in, ou, tl) \
        {(ou).qn=(in).x_id;(ou).tn=(in).y_id;(ou).rev=(in).y_pos_strand;\
        (ou).qs=(in).x_pos_s;(ou).qe=(in).x_pos_e+1;\
        if((in).y_pos_strand){(ou).ts=(tl)-(in).y_pos_e-1;(ou).te=(tl)-(in).y_pos_s;}\
        else {(ou).ts=(in).y_pos_s;(ou).te=(in).y_pos_e+1;}}





typedef struct {
	int64_t qi, qs, qe;
	int64_t ti;
    int64_t qis, qie, tis, tie;
    int64_t wi, wn;
	int64_t tot, lc, cql, ctl, ci, lerr;
    overlap_region *z;
} pe_cigar_iter_t;

#define citer_end(m)  ((m).qi>=(m).qe)

/**
inline void pop_citer(pe_cigar_iter_t *it)
{
    it->lc = it->lerr = -1; it->cql = it->ctl = 0; 
    it->qis = it->qie; it->tis = it->tie;
    if(it->wi >= it->wn) return;
    int64_t ws, we, os, oe, ovlp, ql, tl, werr; 
    window_list *m; bit_extz_t ez;

    while (it->wi < it->wn) {
        m = &(it->z->w_list.a[it->wi]); 
        ws = m->x_start; we = m->x_end+1;
        os = MAX(it->qs, ws); oe = MIN(it->qe, we);
        ovlp = ((oe>os)? (oe-os):0);

        if(ovlp) {
            ql = m->x_end+1-m->x_start;
            tl = m->y_end+1-m->y_start;
            if((is_ualn_win((*m))) || (is_est_aln((*m)))) {
                it->lc = 4;//not an ordinary cigar
                if(is_ualn_win((*m))) { //unmapped
                    werr = gen_err_unaligned(ql, tl);
                } else {
                    werr = m->error;//shared window
                    if(!werr) it->lc = 0;///treat it as match
                }
                it->lerr = werr;
                if(ovlp < ql) {
                    werr = (((double)ovlp)/((double)ql))*((double)werr);
                }
                //skip the whole window
                it->tot += werr; 
                it->qis = it->qi; it->qi = it->qie = m->x_end+1; it->cql = it->qie - it->qis;
                it->tis = it->ti; it->ti = it->tie = m->x_start; it->ctl = it->tie - it->tis;
                it->ci = 0; it->wi++;
                return;
            } else {
                if(it->qi < m->x_start || it->ti < m->y_start) {
                    it->qi = m->x_start; it->ti = m->y_start; it->ci = 0;
                }
                set_bit_extz_t(ez, (*(it->z)), it->wi);
                // if(ovlp == ql) {
                //     //skip the whole window
                //     err += m->error; xk = m->x_end+1; ck = m->clen;
                // } else {
                //     set_bit_extz_t(ez, (*z), wk); 
                //     err += retrieve_cigar_err(&ez, os, oe, &xk, &ck);
                // }                
            }
        }
    }
}
**/













inline void push_rphase(overlap_region *z, uint64_t zs, uint64_t ze, uint64_t zerr, uint64_t is_pri)
{
    uint64_t pe, p_pri; window_list *p; 
    pe = z->x_pos_s; p_pri = 0;
    if(z->w_list.n) {
        pe = z->w_list.a[z->w_list.n-1].x_end;
        p_pri = z->w_list.a[z->w_list.n-1].cidx;
    }

    // fprintf(stderr, "+[M::%s] zs::%lu, ze::%lu, zerr::%lu, pe::%lu, is_pri::%lu, p_pri::%lu, wn::%u\n", 
    // __func__, zs, ze, zerr, pe, is_pri, p_pri, (uint32_t)z->w_list.n);

    if(pe <= zs) {///not overlap with [zs, ze)
        if(pe < zs) {
            kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
            p->clen = 0; p->x_start = pe; p->x_end = zs; p->cidx = 0;
        }
        kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
        p->clen = zerr; p->x_start = zs; p->x_end = ze; p->cidx = is_pri;
    } else {
        if(pe <= ze) {///prefix-suffix overlap
            if(is_pri) {
                if(p_pri) {
                    kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
                    p->clen = zerr; p->cidx = is_pri;
                    p->x_start = zs; p->x_end = ze;
                    z->w_list.a[z->w_list.n-2].x_end = zs;///trim the previous primary window
                } else {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    p->clen = zerr; p->cidx = is_pri;
                    p->x_end = ze;
                }
            } else {///if current is not primary
                p = &(z->w_list.a[z->w_list.n-1]);
                p->x_end = ze;
            }
        } else {//pe > ze; [zs, ze) is contained within [ps, pe)
            if(is_pri) {///is_pri == 0, no need to do anything
                if(!p_pri) {
                    p = &(z->w_list.a[z->w_list.n-1]);
                    p->clen = zerr; p->cidx = is_pri;
                } else {
                    kv_pushp(window_list, z->w_list, &p); memset(p, 0, sizeof((*p)));
                    p->clen = zerr; p->cidx = is_pri; 
                    p->x_start = zs; p->x_end = z->w_list.a[z->w_list.n-2].x_end;
                    z->w_list.a[z->w_list.n-2].x_end = zs;///trim the previous primary window
                }
            }
        }
    }
    // fprintf(stderr, "-[M::%s] zs::%lu, ze::%lu, zerr::%lu, pe::%lu, is_pri::%lu, p_pri::%lu, wn::%u, x::[%d, %d), sec::%u\n", 
    // __func__, zs, ze, zerr, pe, is_pri, p_pri, (uint32_t)z->w_list.n, 
    // z->w_list.a[z->w_list.n-1].x_start, z->w_list.a[z->w_list.n-1].x_end, z->w_list.a[z->w_list.n-1].clen);
}








inline uint32_t ovlp_win_check(overlap_region *z, uint32_t id0, uint32_t id1, int64_t max_lgap, double small_bw_rate, int64_t min_small_bw)
{
    if(id0 == (uint32_t)-1 || id1 == (uint32_t)-1) return 1;
    int64_t qs0, qe0, ts0, te0, qs1, qe1, ts1, te1, err, dd, dm, dq, dr;
    if(!get_win_aln(z, id0, &ts0, &te0, &err)) return 1; 
    qs0 = z->w_list.a[id0].x_start; qe0 = z->w_list.a[id0].x_end+1;
    if(!get_win_aln(z, id1, &ts1, &te1, &err)) return 1; 
    qs1 = z->w_list.a[id1].x_start; qe1 = z->w_list.a[id1].x_end+1;

    if(qs1 < qs0 || qe1 < qe0) return 0; 
    if(ts1 < ts0 || te1 < te0) return 0; 
    dq = qe1 - qs0; dr = te1 - ts0; dd = dq>=dr? ((dq)-(dr)): ((dr)-(dq));
    if((ts1 < te0) && (qe0 == qs1)) {//has overlap in y
        dm = dq>=dr?dr:dq;
        if((dd > (dm*small_bw_rate)) && (dd > min_small_bw)) return 0;
    } else {
        if(dd > max_lgap) return 0;
    }
    return 1;
}
















///q[2], t[2]






 /**
void cal_simi_ul_ov_t(overlap_region *z, ul_ov_t *o, kv_ul_ov_t *res, int64_t ql)
{
    int64_t beg_q[2], beg_t[2], end_q[2], end_t[2], wk[2], wq[2], wt[2], kbeg, kend; 
    int64_t k, tot_l, qs, qe, ts, te, os, oe, ovlp, aln_l, wts, wte, wtl, wn; ul_ov_t *p;
    qs = 0; qe = ql-1; ts = o->ts; te = o->te-1; //[qs, qe] && [ps, pe]
    k = (o->qs==(uint32_t)-1)?(-1):(o->qs);
    k = get_win_yoff(z, ts, k, qs, qe, ts, te, beg_q, beg_t, &(wq[0]), &(wt[0])); wk[0] = k; 
    k = o->qe;
    k = get_win_yoff(z, te, k, qs, qe, ts, te, end_q, end_t, &(wq[1]), &(wt[1])); wk[1] = k;
    o->sec = 0;
    kbeg = ((wk[0]>=0)?(wk[0]>>1):(0)); kend = ((wk[1]>=0)?(wk[1]>>1):(0)); aln_l = 0;
    for (k = kbeg; k <= kend; k++) {
        wts = z->w_list.a[k].y_start; wte = z->w_list.a[k].y_end; wtl = wte + 1 - wts;
        os = MAX(wts, ts); oe = MIN(wte, te) + 1;
        ovlp = ((oe>os)? (oe-os):0); aln_l += ovlp;
        assert(ovlp > 0);
        if(ovlp == wtl) {

        } else {
            
        }
    }
    if(aln_l < te+1-ts) {
        wn = z->w_list.n;
        assert(wk[0] == -1 || wk[1] == wn-1);
    }

    if(wk[0] == wk[1]) {///one window cover the whole [ts, te]
        kv_pushp(ul_ov_t, *res, &p); p->ts = ts; p->te = te;
    }
    if(wk[0] < 0) {
        
    }
}
**/

///[ts, te) ->  this is the reverse coordinates of t, not the original coordinates of t



int64_t return_t_chain(overlap_region *z, Candidates_list *cl)
{
    int64_t i, cn = cl->length, scn; uint64_t pid; k_mer_hit *ca;

    // if(z->y_id == 30129) {
    //     fprintf(stderr, "\n-0-[M::%s]\tutg%.6ul\tx::[%u,\t%u)\t%c\tutg%.6ul\ty::[%u,\t%u)\n", 
    //         __func__, z->x_id+1, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand], z->y_id+1, z->y_pos_s, z->y_pos_e+1);
    //     i = z->shared_seed; pid = cl->list[i].readID;
    //     for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++) {
    //         fprintf(stderr, "i::%ld[M::%s]\treadID::%u\tself_offset::%u\toffset::%u\t%c\n", 
    //             i, __func__, cl->list[i].readID, cl->list[i].self_offset, cl->list[i].offset, 
    //             "+-"[cl->list[i].strand]);
    //     }
    // }



    i = z->shared_seed; pid = cl->list[i].readID;
    for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++);
    scn = i - z->shared_seed; ca = cl->list+z->shared_seed;
    i = lchain_refine(ca, scn, ca, &(cl->chainDP), 50, 5000, 512, 16); cn = i;
    for (; i < scn; i++) ca[i].readID = ((uint32_t)(0x7fffffff));
    // if(z->y_id == 30129) fprintf(stderr, "\n-a-[M::%s]\tcn::%ld\n", __func__, cn);


    // if(z->y_id == 30129) {
    //     fprintf(stderr, "\n-1-[M::%s]\tutg%.6ul\tx::[%u,\t%u)\t%c\tutg%.6ul\ty::[%u,\t%u)\n", 
    //         __func__, z->x_id+1, z->x_pos_s, z->x_pos_e+1,
	// 		"+-"[z->y_pos_strand], z->y_id+1, z->y_pos_s, z->y_pos_e+1);
    //     i = z->shared_seed; pid = cl->list[i].readID; cn = cl->length;
    //     for (; i < cn && cl->list[i].readID == pid && cl->list[i].readID != ((uint32_t)(0x7fffffff)); i++) {
    //         fprintf(stderr, "i::%ld[M::%s]\treadID::%u\tself_offset::%u\toffset::%u\t%c\n", 
    //             i, __func__, cl->list[i].readID, cl->list[i].self_offset, cl->list[i].offset, 
    //             "+-"[cl->list[i].strand]);
    //     }
    // }
    return cn;
}

// int64_t gen_mix_tchain(Candidates_list *cl, int64_t kidx, int64_t kn, ul_ov_t *oa, int64_t on)
// {
//     int64_t rcn = cl->length, kk, ok;
//     kv_resize_cl(k_mer_hit, *cl, rcn+kn);
//     kk = ok = 0;
//     while(kk < kn && ok < on) {
//         if(cl->list[kidx+kk].offset)
//     }
// }



//return [rq, rt]


#define update_ul_ov_t_coor(z, qbeg, qend, tbeg, tend) do {\
        if(((int64_t)(z).qs) > (qbeg)) (z).qs = (qbeg);\
        if(((int64_t)(z).qe) < (qend)) (z).qe = (qend);\
        if(((int64_t)(z).ts) > (tbeg)) (z).ts = (tbeg);\
        if(((int64_t)(z).te) < (tend)) (z).te = (tend);\
	} while (0)

///[ps, pe) && [ts, te])


inline void update_trace_idx(rtrace_t *tc, int64_t wid, int64_t wid_s, int64_t wid_e, 
int64_t qs, int64_t qe, int64_t ts, int64_t te)
{
    if(qs < tc->c_qs || ts < tc->c_ts) {
        tc->c_qs = qs; 
        tc->c_ts = ts;
        tc->c_wsid = wid;
        tc->c_wsii = wid_s;
    }

    if(qe > tc->c_qe || te > tc->c_te) {
        tc->c_qe = qe; 
        tc->c_te = te;
        tc->c_weid = wid;
        tc->c_weii = wid_e;
    }
}













//get error within [qe, ql)


//get error within [qe, ql)


///[wsid, weid) && [ts, te)


///[ts, te) ->  this is the reverse coordinates of t, not the original coordinates of t












uint64_t gen_hc_fast_cigar(overlap_region *z, Candidates_list *cl, All_reads *rref, int64_t wl, char *qstr, UC_Read* tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t ql, int64_t rid, int64_t khit, int64_t *re)
{
    return_t_chain(z, cl);
    gen_hc_fast_cigar0(z, cl, wl, rref, qstr, tu, exz, aux_o, e_rate, ql, rid, khit, re);
    return 1;
}


void append_cigar(window_list *idx, window_list_alloc *res, uint16_t c, uint32_t l)
{
    if(l <= 0) return;
    uint16_t c0 = (uint16_t)-1; uint32_t l0 = 0; 
    ///last item of old cigar
    if(idx->clen > 0) {
        c0 = (res->c.a[res->c.n-1]>>14); 
        l0 = (res->c.a[res->c.n-1]&(0x3fff));
    }
    ///first item of new cigar
    if(c0 == c) {l += l0; res->c.n--; idx->clen--;}

    push_trace(((asg16_v *)(&(res->c))), c, l); 
    idx->clen = res->c.n-idx->cidx;
}

uint16_t adjust_gap(window_list *idx, window_list_alloc *res, char *pstr, char *tstr, int64_t pi, int64_t ti, uint16_t op0, asg16_v* buf, int64_t *rd_err)
{
    (*rd_err) = 0;

    if(idx->clen == 0) {
        append_cigar(idx, res, op0, 1);
        return 0;///no move
    }
    if(op0 != 2 && op0 != 3) return 0;///no move
    if(op0 == 2) {///more p -> y
        ti--;
    } else if(op0 == 3) {///more t -> x
        pi--;
    }

    // if(z->y_id == 3199 && z->x_id == 3196) {
    //     fprintf(stderr, "\n[M::%s]\tqi::%ld\tti::%ld\n", __func__, ti, pi);
    // }

    uint16_t *ca = res->c.a+idx->cidx; 
    int64_t ci = idx->clen; int64_t op, cl, k, l[2]; uint16_t p, ff;
    for (ci--, buf->n = ff = 0; ci >= 0; ci--) {
        op = ca[ci]>>14; cl = (ca[ci]&(0x3fff)); 
        // if(z->y_id == 3199 && z->x_id == 3196) {
        //     fprintf(stderr, "[M::%s]\tqi::%ld\tti::%ld\tci::%ld\tcl::%ld\top::%ld\n", __func__, ti, pi, ci, cl, op);
        // }
        if(op == 2 || op == 3) {
            p = op0; p <<= 14; p += 1; kv_push(uint16_t, *buf, p);
            l[0] = cl;
            p = op; p <<= 14; p += l[0]; kv_push(uint16_t, *buf, p);
            break;
        } else if(op == 0) {        
            for (k = cl-1, l[0] = l[1] = 0; k >= 0; k--, pi--, ti--) {
                if(pstr[pi] != (tstr[ti])) break;
                // if(op == 1) {
                //     if(!(pstr[pi] != (tstr[ti]))) {
                //         fprintf(stderr, "[M::%s] xid::%u, yid::%u, qi::%ld, ti::%ld\n", __func__, z->x_id, z->y_id, ti, pi);
                //     }
                //     assert(pstr[pi] != (tstr[ti]));
                // }
            }
            l[1] = cl - k - 1; l[0] = k + 1;
            if(l[1] > 0) {
                p = op; p <<= 14; p += l[1]; kv_push(uint16_t, *buf, p); ff = 1;
            }

            if(l[0] > 0) {
                p = op0; p <<= 14; p += 1; kv_push(uint16_t, *buf, p);
            }

            if(l[0] > 0) {
                p = op; p <<= 14; p += l[0]; kv_push(uint16_t, *buf, p);
            }
        } else {///op == 1; it is possible since cigar is not optimal
            for (k = cl-1, l[0] = cl, l[1] = 0; k >= 0; k--, pi--, ti--) {
                if(pstr[pi] == (tstr[ti])) {
                    l[1] = l[0] - k - 1;
                    l[0] = k;

                    if(l[1] > 0) {
                        p = op; p <<= 14; p += l[1]; kv_push(uint16_t, *buf, p); ff = 1;///push unmatch
                    }
                    p = 0; p <<= 14; p += 1; kv_push(uint16_t, *buf, p);///push match
                    (*rd_err)++;
                }
            }

            if(l[0] > 0) {
                p = op; p <<= 14; p += l[0]; kv_push(uint16_t, *buf, p);
            }
            l[0] = 0;
        }

        // fprintf(stderr, "[M::%s] ci::%ld, l[0]::%ld, l[1]::%ld\n", __func__, ci, l[0], l[1]);

        if(l[0] > 0) break;
    }

    // fprintf(stderr, "[M::%s] ff::%u, ci::%ld, buf->n::%lu\n", __func__, ff, ci, (uint64_t)buf->n);

    if(!ff) {///no move
        append_cigar(idx, res, op0, 1);
        return 0;///no move
    } else if(ci >= 0) {
        // if(!(idx->cidx + idx->clen == res->c.n)) {
        //     fprintf(stderr, "[M::%s] idx->cidx::%lu, idx->clen::%lu, res->c.n::%lu\n", __func__, (uint64_t)idx->cidx, (uint64_t)idx->clen, (uint64_t)res->c.n);
        // }
        assert(idx->cidx + idx->clen == res->c.n);
        idx->clen = ci; res->c.n = idx->cidx + idx->clen;
        for (k = ((int64_t)buf->n)-1; k >= 0; k--) {
            op = buf->a[k]>>14; cl = (buf->a[k]&(0x3fff)); 
            append_cigar(idx, res, op, cl);
        }
    } else {
        idx->clen = 0; res->c.n = idx->cidx + idx->clen;
        append_cigar(idx, res, op0, 1);
        for (k = ((int64_t)buf->n)-1; k >= 0; k--) {
            op = buf->a[k]>>14; cl = (buf->a[k]&(0x3fff)); 
            append_cigar(idx, res, op, cl);
        }
    }
    return 1;///no move
}

uint16_t ajust_end_cigar(window_list *idx, window_list_alloc *res)
{
    uint16_t *ca = res->c.a+idx->cidx, p, rr = 0; int64_t ci, cn = idx->clen, op, cl;
    if(cn <= 0) return rr;
    for (ci = 0; ci < cn; ci++) {
        op = ca[ci]>>14; cl = (ca[ci]&(0x3fff)); 
        if(op != 1) break;
        p = 3; p <<= 14; p += cl; ca[ci] = p; idx->y_start += cl; rr = 1;
    }

    for (ci = cn-1; ci >= 0; ci--) {
        op = ca[ci]>>14; cl = (ca[ci]&(0x3fff)); 
        if(op != 1) break;
        p = 3; p <<= 14; p += cl; ca[ci] = p; idx->y_end -= cl; rr = 1;
    }

    // if(rr) fprintf(stderr, "[M::%s] rr::%u\n", __func__, rr);

    return rr;
}

uint16_t move_wins(overlap_region *z, uint32_t wid, overlap_region *aux, All_reads *rref, char *tstr, char *pstr, UC_Read *pu, asg16_v* buf, int64_t *tot_re)
{
    // if(z->y_id == 3199 && z->x_id == 3196) {
    //     fprintf(stderr, "\n[M::%s] x_id::%u, y_id::%u, x::[%u, %u), y::[%u, %u)\n", __func__, z->x_id, z->y_id, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1);
    //     fprintf(stderr, "[M::%s] wid::%u, x::[%d, %d), y::[%d, %d), err::%d\n", __func__, wid, z->w_list.a[wid].x_start, z->w_list.a[wid].x_end + 1, z->w_list.a[wid].y_start, z->w_list.a[wid].y_end + 1, z->w_list.a[wid].error);
    // }
    if(is_ualn_win((z->w_list.a[wid]))) {
        kv_push(window_list, aux->w_list, (z->w_list.a[wid]));
        return 0;///no move
    }

    window_list *p = NULL; char *tseq = tstr, *pseq = pstr;
    bit_extz_t ez; set_bit_extz_t(ez, (*z), wid);
    kv_pushp(window_list, aux->w_list, &p);
    p->x_start = ez.ts; p->x_end = ez.te;
    p->y_start = ez.ps; p->y_end = ez.pe;
    p->error_threshold = 0; p->error = ez.err;///single round of alignment cannot have INT16_MAX errors
    p->cidx = aux->w_list.c.n; p->clen = 0; ///aux->w_list.c.n += ez.cigar.n;

    if(ez.err == 0) {
        push_wcigar(p, &(aux->w_list), &ez);
        return 0;///no move
    }

    int64_t pi = ez.ps, ti = ez.ts/**, err = 0**/, cl, op, k, rr = 0, re; uint64_t ci = 0, mm = 0; 
    // char cm[4]; cm[0] = 'M'; cm[1] = 'S'; cm[2] = 'I'; cm[3] = 'D'; 
    // if(z->y_id == 3199 && z->x_id == 3196) {
    //     fprintf(stderr, "******\n");
    //     for (ci = 0; ci < ez.cigar.n; ci++) {
    //         op = ez.cigar.a[ci]>>14; cl = (ez.cigar.a[ci]&(0x3fff)); 
    //         fprintf(stderr, "%ld%c(q::[%ld,%ld))(t::[%ld,%ld))(ck::%ld)\n", cl, cm[op], ti, ti + (((op<2)||(op==3))?(cl):(0)), pi, pi + (((op<2)||(op==2))?(cl):(0)), ci);

    //         if(op < 2) {
    //             pi+=cl; ti+=cl;
    //         } else {
    //             if(op == 2) {///more p -> y
    //                 pi+=cl;
    //             } else if(op == 3) {///more t -> x
    //                 ti+=cl;
    //             }
    //         }
    //     }
    //     pi = ez.ps; ti = ez.ts;
    //     fprintf(stderr, "******\n");
    // }

    for (ci = mm = 0; ci < ez.cigar.n; ci++) {
        op = ez.cigar.a[ci]>>14; ///cl = (ez.cigar.a[ci]&(0x3fff)); 
        if(op < 2) {
            mm = 1;
        } else if(mm) {
            break;
        }
        // if(op == 2 || op == 3) break;
    }
    if(ci >= ez.cigar.n) {
        push_wcigar(p, &(aux->w_list), &ez);
        if(ajust_end_cigar(p, &(aux->w_list))) rr = 1;
        return rr;
    }

    if(!pseq) {
        if(z->y_pos_strand) {
            recover_UC_Read_RC(pu, rref, aux->y_id);
        } else {
            recover_UC_Read(pu, rref, aux->y_id);
        }
        pseq = pu->seq;
    }
    
    kv_resize(uint16_t, aux->w_list.c, (aux->w_list.c.n + ez.cigar.n));

    for (ci = mm = 0; ci < ez.cigar.n; ci++) {
        op = ez.cigar.a[ci]>>14; cl = (ez.cigar.a[ci]&(0x3fff)); 
        // if(z->y_id == 24 && z->x_id == 25) {
        //     fprintf(stderr, "[M::%s] op::%ld, cl::%ld, pi::%ld, ti::%ld\n", __func__, op, cl, pi, ti);
        //     fprintf(stderr, "[M::%s] p->cidx::%lu, p->clen::%lu, cc->n::%lu\n", __func__, (uint64_t)p->cidx, (uint64_t)p->clen, (uint64_t)aux->w_list.c.n);
        // }
        // if(z->y_id == 3199 && z->x_id == 3196) {
        //     fprintf(stderr, "%ld%c(q::[%ld,%ld))(t::[%ld,%ld))(ck::%ld)\tmm::%lu\n", cl, cm[op], ti, ti + (((op<2)||(op==3))?(cl):(0)), pi, pi + (((op<2)||(op==2))?(cl):(0)), ci, mm);
        // }

        if(op < 2) {
            append_cigar(p, &(aux->w_list), op, cl);
            pi+=cl; ti+=cl; mm = 1;
        } else {
            if(mm == 0) {
                append_cigar(p, &(aux->w_list), op, cl);
                if(op == 2) {///more p -> y
                    pi+=cl;
                } else if(op == 3) {///more t -> x
                    ti+=cl;
                }
            } else {
                for (k = 0; k < cl; k++) {
                    // if(z->y_id == 3199 && z->x_id == 3196) {
                    //     fprintf(stderr, "+++k::%ld\tqi::%ld\tti::%ld\n", k, ti, pi);
                    // }
                    if(adjust_gap(p, &(aux->w_list), pseq, tseq, pi, ti, op, buf, &re)) {
                        rr = 1; p->error -= re; (*tot_re) += re;
                    }
                    if(op == 2) {///more p -> y
                        pi++;
                    } else if(op == 3) {///more t -> x
                        ti++;
                    }
                }
            }
        }
    }

    // if(rr) fprintf(stderr, "[M::%s] rr::%ld\n", __func__, rr);

    if(ajust_end_cigar(p, &(aux->w_list))) rr = 1;
    return rr;
}

void reassign_gaps(overlap_region *z, overlap_region *aux, char* qstr, int64_t ql, char* tstr, int64_t tl, All_reads *rref, UC_Read* tu, asg16_v* buf)
{
    // if(z->y_id != 24 || z->x_id != 25) return;
    if(z->non_homopolymer_errors == 0) return;
    aux->w_list.n = aux->w_list.c.n = 0; 
    aux->y_id = z->y_id; aux->y_pos_strand = z->y_pos_strand;
    aux->x_pos_s = z->x_pos_s; aux->x_pos_e = z->x_pos_e;
    aux->y_pos_s = z->y_pos_s; aux->y_pos_e = z->y_pos_e;

    int64_t k, z_n = z->w_list.n, rr = 0, re = 0;
    for (k = 0; k < z_n; k++) {
        if(move_wins(z, k, aux, rref, qstr, tstr, tu, buf, &re)) rr = 1;
    }
    ///update z by aux_o
    if(rr) update_overlap_region(z, aux, ql, ((rref)?(Get_READ_LENGTH((*rref), z->y_id)):(tl)));
    z->non_homopolymer_errors -= re;

    // if(z->y_id == 1) {
    //     fprintf(stderr, "[M::%s] rr::%ld\tx_id::%u\ty_id::%u\tx::[%u, %u)\ty::[%u, %u)\tz_n::%ld\n", 
    //             __func__, rr, z->x_id, z->y_id, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1, z_n);
    //     fprintf(stderr, "qstr(%lld)::%.*s\n", qu->length, (int32_t)(qu->length), qu->seq);
    //     fprintf(stderr, "tstr(%lld)::%.*s\n", tu->length, (int32_t)(tu->length), tu->seq);
    //     ///debug
        
        // if(z->y_pos_strand) {
        //     recover_UC_Read_RC(tu, rref, z->y_id);
        // } else {
        //     recover_UC_Read(tu, rref, z->y_id);
        // }
        // if(tstr) {
        //     bit_extz_t ez; 
        //     for (k = 0; k < z_n; k++) {
        //         if(is_ualn_win((z->w_list.a[k]))) continue;
        //         set_bit_extz_t(ez, (*z), k);
        //         if(!cigar_check(tstr, qstr, &ez)) {
        //             fprintf(stderr, "\n[M::%s] x_id::%u, y_id::%u, x::[%u, %u), y::[%u, %u)\n", __func__, z->x_id, z->y_id, z->x_pos_s, z->x_pos_e + 1, z->y_pos_s, z->y_pos_e + 1);
        //             exit(1);
        //         }
        //     }
        // }
        
    // }
}


uint32_t inline ff_tend(overlap_region *z, int64_t wn, int64_t dn, double dr, double er, int64_t min_err)
{
    int64_t k, zwn = z->w_list.n, err, mm, qi, ci, cn, ql, ws, we, zs, ze, s, e, os, oe, ovlp; bit_extz_t ez; uint32_t cl; uint16_t c;
    zs = z->x_pos_s; ze = z->x_pos_e + 1; ql = ze - zs; 
    if(ql < wn) return 0; if(dn > (ql*dr)) dn = ql*dr; if(dn < wn) return 0; if(ql < dn) return 0;
    
    s = zs; e = zs + dn;
    // if(z->y_id == 27) fprintf(stderr, "-a-[M::%s] tid::%u\t%.*s\tz::[%ld,%ld)\ti::[%ld,%ld)\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), zs, ze, s, e);
    for (k = err = mm = 0, qi = zs; (k < zwn) && (z->w_list.a[k].x_start < e); k++) {
        // if(z->y_id == 27) fprintf(stderr, "-a-[M::%s] tid::%u\t%.*s\tw::[%d,%d)\terr::%d\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), z->w_list.a[k].x_start, z->w_list.a[k].x_end + 1, z->w_list.a[k].error);
        if(!(is_ualn_win(z->w_list.a[k]))) {
            set_bit_extz_t(ez, (*z), k); 
            ci = 0; cn = ez.cigar.n; qi = ez.ts; //ti = ez.ps;
            while (ci < cn && qi < e) {
                ws = qi;
                ci = pop_trace(&(ez.cigar), ci, &c, &cl);
                if(c!=2) qi += cl;
                // if(c!=3) ti += cl; 
                we = qi; 

                if(c == 0) {
                    os = MAX(s, ws); oe = MIN(e, we);
                    ovlp = ((oe>os)? (oe-os):0); mm += ovlp;
                } else {
                    err += cl;
                }
                // assert(is_ovlp_debug(s, e, ws, we, c));

                if((err > min_err) && ((mm + err) > wn) && (err > ((mm + err)*er))) {
                    // fprintf(stderr, "-0-[M::%s] tid::%u\t%.*s\tmm::%ld\terr::%ld\tdn::%ld\twn::%ld\tdif::%f\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), mm, err, dn, wn, er);
                    return 1;
                }
            }
        } else {
            ws = z->w_list.a[k].x_start; we = z->w_list.a[k].x_end + 1;
            err += we - ws;

            // assert(is_ovlp_debug(s, e, ws, we, -1));

            if((err > min_err) && ((mm + err) > wn) && (err > ((mm + err)*er))) {
                // fprintf(stderr, "-1-[M::%s] tid::%u\t%.*s\tmm::%ld\terr::%ld\tdn::%ld\twn::%ld\tdif::%f\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), mm, err, dn, wn, er);
                return 1;
            }
        } 
    }

    s = ze - dn; e = ze;
    for (k = zwn - 1, err = mm = 0, qi = ze; (k >= 0) && ((z->w_list.a[k].x_end + 1) > s); k--) {
        if(!(is_ualn_win(z->w_list.a[k]))) {
            set_bit_extz_t(ez, (*z), k); 
            ci = ((int64_t)ez.cigar.n) - 1; qi = ez.te + 1; //ti = ez.pe + 1;
            while (ci >= 0 && qi > s) {
                we = qi;
                ci = pop_trace_back(&(ez.cigar), ci, &c, &cl);
                if(c!=2) qi -= cl;
                // if(c!=3) ti += cl; 
                ws = qi; 

                if(c == 0) {
                    os = MAX(s, ws); oe = MIN(e, we);
                    ovlp = ((oe>os)? (oe-os):0); mm += ovlp;
                } else {
                    err += cl;
                }

                // assert(is_ovlp_debug(s, e, ws, we, c));

                if((err > min_err) && ((mm + err) > wn) && (err > ((mm + err)*er))) {
                    // fprintf(stderr, "-2-[M::%s] tid::%u\t%.*s\tmm::%ld\terr::%ld\tdn::%ld\twn::%ld\tdif::%f\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), mm, err, dn, wn, er);
                    return 1;
                }
            }
        } else {
            ws = z->w_list.a[k].x_start; we = z->w_list.a[k].x_end + 1;
            err += we - ws;

            // assert(is_ovlp_debug(s, e, ws, we, -1));

            if((err > min_err) && ((mm + err) > wn) && (err > ((mm + err)*er))) {
                // fprintf(stderr, "-3-[M::%s] tid::%u\t%.*s\tmm::%ld\terr::%ld\tdn::%ld\twn::%ld\tdif::%f\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), mm, err, dn, wn, er);
                return 1;
            }
        } 
    }
    
    return 0;
}

uint32_t inline ff_lunalign(overlap_region *z, double erate, double gap_rate, int64_t max_gap)
{
    if((z->w_list.n == 1) && (!(is_ualn_win(z->w_list.a[0]))) 
        && (z->w_list.a[0].x_start == ((int64_t)z->x_pos_s)) && (z->w_list.a[0].x_end == ((int64_t)z->x_pos_e)) 
            && (z->w_list.a[0].y_start == ((int64_t)z->y_pos_s)) && (z->w_list.a[0].y_end == ((int64_t)z->y_pos_e))) {
                return 1;
    }

    int64_t k, zwn = z->w_list.n, zq, zt, wq, wt, tot_e, tot_g, ql, tl; 

    // fprintf(stderr, "[M::%s]\n", __func__);
    
    zq = z->x_pos_s; zt = z->y_pos_s; tot_e = tot_g = 0; 
    for (k = 0; k < zwn; k++) {
        // fprintf(stderr, "[M::%s]\twk::%ld\tq::[%d, %d)\tt::[%d, %d)\terr::%d\n", __func__, 
        //     k, z->w_list.a[k].x_start, z->w_list.a[k].x_end + 1, z->w_list.a[k].y_start, z->w_list.a[k].y_end + 1, z->w_list.a[k].error);
        if(is_ualn_win(z->w_list.a[k])) continue;
        wq = z->w_list.a[k].x_start; 
        wt = z->w_list.a[k].y_start;

        if(wq != zq) tot_g += ((wq>=zq)?(wq-zq):(zq-wq));
        if(wt != zt) tot_g += ((wt>=zt)?(wt-zt):(zt-wt));

        zq = z->w_list.a[k].x_end + 1; 
        zt = z->w_list.a[k].y_end + 1;
        tot_e += z->w_list.a[k].error;

        // fprintf(stderr, "[M::%s]\twk::%ld\tq::[%d, %d)\tt::[%d, %d)\terr::%d\n", __func__, 
        //     k, z->w_list.a[k].x_start, z->w_list.a[k].x_end + 1, z->w_list.a[k].y_start, z->w_list.a[k].y_end + 1, z->w_list.a[k].error);
    }

    wq = z->x_pos_e + 1;
    wt = z->y_pos_e + 1;
    if(wq != zq) tot_g += ((wq>=zq)?(wq-zq):(zq-wq));
    if(wt != zt) tot_g += ((wt>=zt)?(wt-zt):(zt-wt));

    // fprintf(stderr, "[M::%s]\t%.*s(id::%u)\tq::[%u, %u)\t%.*s(id::%u)\tt::[%u, %u)\tre::%ld\trg::%ld\n", __func__, (int)Get_NAME_LENGTH(R_INF, z->x_id), Get_NAME(R_INF, z->x_id), z->x_id, z->x_pos_s, z->x_pos_e + 1, 
    // (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), z->y_id, z->y_pos_s, z->y_pos_e + 1, tot_e, tot_g);

    if(!tot_g) return 1;

    // fprintf(stderr, "-0-[M::%s]\n", __func__);

    if(tot_g > max_gap) return 0;

    // fprintf(stderr, "-1-[M::%s]\n", __func__);

    ql = z->x_pos_e + 1 - z->x_pos_s;
    tl = z->y_pos_e + 1 - z->y_pos_s;

    if((tot_g > (ql*gap_rate)) || (tot_g > (tl*gap_rate))) return 0;

    // fprintf(stderr, "-2-[M::%s]\n", __func__);

    tot_e += tot_g;

    if((tot_e > (ql*erate)) || (tot_e > (tl*erate))) return 0;

    // fprintf(stderr, "-3-[M::%s]\n", __func__);
    
    return 1;
}

void gen_hc_r_alin(overlap_region_alloc* ol, Candidates_list *cl, All_reads *rref, UC_Read* qu, UC_Read* tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t wl, int64_t rid, int64_t khit, int64_t move_gap, asg16_v* buf, uint8_t chem_drop, double align_gap_rate, int64_t align_gap_max)
{
    uint64_t i, bs, k, ql = qu->length; Window_Pool w; double err, e_max, rr; int64_t re;
    overlap_region t; overlap_region *z; //asg64_v iidx, buf, buf1;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;
    // if(ol->length && ol->list[0].x_id == 19350) e_rate = 0.1;

    ///base alignment
    err = e_rate; e_max = err * 1.5;
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1); 
    // fprintf(stderr, "[M::%s] window_length::%lld, err::%f\n", __func__, w.window_length, err);

    for (i = k = 0; i < ol->length; i++) {
        z = &(ol->list[i]); z->shared_seed = z->non_homopolymer_errors;///for index

        // if(z->x_id == 19350 && z->y_id == 19324) fprintf(stderr, "-z-[M::%s] tid::%u\t%.*s\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id));
        
        if(!align_hc_ed_post_extz(z, rref, qu->seq, tu->seq, exz, err, w.window_length, OVERLAP_THRESHOLD_HIFI_FILTER, 0, NULL)) continue;

        // if(z->x_id == 19350 && z->y_id == 19324) fprintf(stderr, "-m-[M::%s] tid::%u\t%.*s\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id));

        rr = gen_extend_err_exz(z, NULL, NULL, rref, qu->seq, tu->seq, exz, NULL, w.window_length, -1, err, (e_max+0.000001), THRESHOLD_MAX_SIZE, 0, &re);
        z->is_match = 0; 

        // if(z->x_id == 19350 && z->y_id == 19324) fprintf(stderr, "-0-[M::%s] tid::%u\t%.*s\trr::%f\tre::%ld\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, re);

        if (rr > err) continue;
        z->non_homopolymer_errors = re;

        if(!gen_hc_fast_cigar(z, cl, rref, w.window_length, qu->seq, tu, exz, aux_o, e_rate, ql, rid, khit, &re)) continue;

        if((align_gap_max >= 0) && (!ff_lunalign(z, err, align_gap_rate, align_gap_max))) continue;

        if(chem_drop && ff_tend(z, 384, 2000, 0.1, (((e_rate*10)<0.36)?(e_rate*10):(0.36)), 128)) continue;


        // if(z->x_id == 3196 && z->y_id == 3199) fprintf(stderr, "-1-[M::%s] tid::%u\t%.*s\trr::%f\tre::%ld\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, re);

        reassign_gaps(z, aux_o, qu->seq, ql, NULL, -1, rref, tu, buf);

        // fprintf(stderr, "-2-[M::%s] tid::%u\t%.*s\trr::%f\tre::%u\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, z->non_homopolymer_errors);

        // if(z->x_id == 19350) fprintf(stderr, "-1-[M::%s] tid::%u\t%.*s\trr::%f\terr::%u\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, z->non_homopolymer_errors);

        if(k != i) {
            t = ol->list[k];
            ol->list[k] = ol->list[i];
            ol->list[i] = t;
        }
        z = &(ol->list[k++]); z->is_match = 1; ///z->non_homopolymer_errors = re;
        z->strong = z->without_large_indel = 0;
    }
    ol->length = k;
    if(ol->length <= 0) return;
}


void gen_hc_r_alin_nec(overlap_region_alloc* ol, Candidates_list *cl, All_reads *rref, UC_Read* qu, UC_Read* tu, bit_extz_t *exz, overlap_region *aux_o, double e_rate, int64_t wl, int64_t rid, int64_t khit, int64_t move_gap, asg16_v* buf, uint8_t chem_drop, double align_gap_rate, int64_t align_gap_max)
{
    uint64_t i, bs, k, ql = qu->length; Window_Pool w; double err, e_max, rr; int64_t re;
    overlap_region t; overlap_region *z; //asg64_v iidx, buf, buf1;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;
    // if(ol->length && ol->list[0].x_id == 19350) e_rate = 0.1;

    ///base alignment
    err = e_rate; e_max = err * 1.5;
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1); 
    // fprintf(stderr, "[M::%s] window_length::%lld\n", __func__, w.window_length);

    ///debug for memory
    // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

    for (i = k = 0; i < ol->length; i++) {
        z = &(ol->list[i]); 
        if(z->is_match != 1) {
            z->shared_seed = z->non_homopolymer_errors;///for index

            // if(z->x_id == 19350 && z->y_id == 19324) fprintf(stderr, "-z-[M::%s] tid::%u\t%.*s\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id));
            
            if(!align_hc_ed_post_extz(z, rref, qu->seq, tu->seq, exz, err, w.window_length, OVERLAP_THRESHOLD_HIFI_FILTER, 0, NULL)) continue;

            // if(z->x_id == 19350 && z->y_id == 19324) fprintf(stderr, "-m-[M::%s] tid::%u\t%.*s\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id));

            rr = gen_extend_err_exz(z, NULL, NULL, rref, qu->seq, tu->seq, exz, NULL, w.window_length, -1, err, (e_max+0.000001), THRESHOLD_MAX_SIZE, 0, &re);
            z->is_match = 0; 

            // if(z->x_id == 19350 && z->y_id == 19324) fprintf(stderr, "-0-[M::%s] tid::%u\t%.*s\trr::%f\tre::%ld\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, re);

            if (rr > err) continue;
            z->non_homopolymer_errors = re;

            ///debug for memory
            // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

            if(!gen_hc_fast_cigar(z, cl, rref, w.window_length, qu->seq, tu, exz, aux_o, e_rate, ql, rid, khit, &re)) continue;

            if((align_gap_max >= 0) && (!ff_lunalign(z, err, align_gap_rate, align_gap_max))) continue;

            if(chem_drop && ff_tend(z, 384, 2000, 0.1, (((e_rate*10)<0.36)?(e_rate*10):(0.36)), 128)) continue;

            // if(z->x_id == 3196 && z->y_id == 3199) fprintf(stderr, "-1-[M::%s] tid::%u\t%.*s\trr::%f\tre::%ld\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, re);
            ///debug for memory
            // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);

            reassign_gaps(z, aux_o, qu->seq, ql, NULL, -1, rref, tu, buf);

            ///debug for memory
            // snprintf(NULL, 0, "dwn::%u\tdcn::%u", (uint32_t)aux_o->w_list.n, (uint32_t)aux_o->w_list.c.n);
        }

        // if(z->x_id == 3196 && z->y_id == 3199) fprintf(stderr, "-2-[M::%s] tid::%u\t%.*s\trr::%f\tre::%u\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, z->non_homopolymer_errors);

        // if(z->x_id == 19350) fprintf(stderr, "-1-[M::%s] tid::%u\t%.*s\trr::%f\terr::%u\n", __func__, z->y_id, (int)Get_NAME_LENGTH(R_INF, z->y_id), Get_NAME(R_INF, z->y_id), rr, z->non_homopolymer_errors);

        if(k != i) {
            t = ol->list[k];
            ol->list[k] = ol->list[i];
            ol->list[i] = t;
        }
        z = &(ol->list[k++]); z->is_match = 1; ///z->non_homopolymer_errors = re;
        z->strong = z->without_large_indel = 0;
    }
    ol->length = k;
    if(ol->length <= 0) return;
}



/**
void ul_raw_lalign_adv(overlap_region_alloc* ol, Candidates_list *cl, const ul_idx_t *uref, All_reads *rdb, const ug_opt_t *uopt, 
        char *qstr, uint64_t ql, UC_Read* qu, UC_Read* tu, Correct_dumy* dumy, bit_extz_t *exz, haplotype_evdience_alloc* hap, 
        kvec_t_u64_warp* v_idx, overlap_region *aux_o, double e_rate, int64_t wl, kv_ul_ov_t *aln, kv_ul_ov_t *aln1,
        int64_t sid, uint64_t khit, st_mt_t *stb, void *km)
{
    uint64_t i, bs, k, aln_occ; Window_Pool w; double err; 
    overlap_region t; overlap_region *z; asg64_v iidx, buf, buf1;
    ol->mapped_overlaps_length = 0;
    if(ol->length <= 0) return;

    ///base alignment
    clear_Correct_dumy(dumy, ol, km); err = e_rate; 
    init_Window_Pool(&w, ql, wl, (int)(1.0/err));
    bs = (w.window_length)+(THRESHOLD_MAX_SIZE<<1)+1;
    resize_UC_Read(tu, bs<<1); 

    if(!aux_o) {
        resize_UC_Read(qu, ql); qu->length = ql; memcpy(qu->seq, qstr, ql); 
        copy_asg_arr(iidx, hap->snp_srt);
        for (i = k = 0, aln->n = 0; i < ol->length; i++) {
            z = &(ol->list[i]); z->shared_seed = z->non_homopolymer_errors;///for index
            align_ul_ed_post_extz(z, uref, NULL, qu->seq, tu->seq, exz, err, w.window_length, -1, 1, km);
            aln_occ = gen_r_aln(uref, z, k, &iidx, OVERLAP_THRESHOLD_FILTER, aln, 1000);
            if(!aln_occ) continue;
            if(k != i) {
                t = ol->list[k];
                ol->list[k] = ol->list[i];
                ol->list[i] = t;
            }
            z = &(ol->list[k++]); z->is_match = 1; 
        }
        copy_asg_arr(hap->snp_srt, iidx);
        ol->length = k;
        if(ol->length <= 0) return;
    } else {
        copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
        ul_gap_filling_local(ol, cl, aln, wl, uref, NULL, NULL, qu->seq, tu, exz, aux_o, &buf, &iidx, err, ql, sid, khit, 1, MAX_LGAP(ql));
        copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);

        copy_asg_arr(iidx, hap->snp_srt);
        gen_aln_local(ol, aln, aln1, &iidx, uref, qu->seq, tu, exz, aux_o, err, ql, OVERLAP_THRESHOLD_FILTER);
        copy_asg_arr(hap->snp_srt, iidx);
    }
    // } else {
    //     // fprintf(stderr, "-[M::%s] on::%lu\n", __func__, ol->length);
    //     if(ol->length <= 1) return;
    //     ///coordinates for all intervals with cov > 1
    //     copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
    //     // fprintf(stderr, "\n[M::%s] iidx_n::%ld\n", __func__, (int64_t)iidx.n);
    //     ul_gap_filling_adv(ol, cl, aln, wl, uref, NULL, NULL, qu->seq, tu, exz, aux_o, &buf, &iidx, err, ql, sid, khit, 1, MAX_LGAP(ql));
    //     copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);

    //     copy_asg_arr(iidx, hap->snp_srt); copy_asg_arr(buf, v_idx->a); copy_asg_arr(buf1, (*stb));
    //     region_phase(ol, uref, uopt, aln, &iidx, &buf, &buf1);
    //     copy_asg_arr(hap->snp_srt, iidx); copy_asg_arr(v_idx->a, buf); copy_asg_arr((*stb), buf1);
    // }
}
**/


