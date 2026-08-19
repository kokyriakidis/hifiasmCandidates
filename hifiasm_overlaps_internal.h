/*
 * hifiasm_overlaps_internal.h
 *
 * Internal (NOT part of the public bridge API) declarations shared between the
 * bridge implementation (hifiasm_overlaps.cpp) and the two overlap emitters
 * (candidates.cpp for the raw pre-alignment set, ecovlp.cpp for the aligned
 * set). It exposes an in-memory "overlap sink": while active, the emitters push
 * overlaps into a process-global buffer instead of writing a PAF file, so
 * hifiasm_detect_overlaps_mem() can return them without the text round-trip.
 *
 * The sink and all its functions use hifiasm's process-global read/option
 * stores and are NOT thread-safe with respect to begin/end/capture; push() is
 * internally serialized so the per-read worker threads may call it
 * concurrently. This matches the existing single-call-at-a-time contract of the
 * bridge entry points.
 */

#ifndef HIFIASM_OVERLAPS_INTERNAL_H
#define HIFIASM_OVERLAPS_INTERNAL_H

#include <stdint.h>
#include "hifiasm_overlaps.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Returns non-zero while an in-memory collection is in progress (i.e. between
 * hifiasm_ovlp_sink_begin() and hifiasm_ovlp_sink_end()). The emitters check
 * this to decide between pushing into the sink and writing PAF text. */
int hifiasm_ovlp_sink_active(void);

/* Append one overlap to the sink. Serialized internally; safe to call from the
 * kt_for/kt_pipeline worker threads. No-op if the sink is not active. Fields
 * mirror hifiasm_overlap_t (see hifiasm_overlaps.h). */
/* cigar/cigar_len: hifiasm packed uint16_t CIGAR tokens for this overlap (may be
 * NULL/0 for the raw candidate path). cigar_t_start: the target anchor in the
 * alignment frame the tokens run in (see hifiasm_overlap_t::cigar_t_start).
 * chain/chain_len: hifiasm's dense native chain anchors for this overlap, each
 * packed (q_start<<32)|t_start (see hifiasm_overlap_t::chain_offset). May be
 * NULL/0 when no chain was captured. */
void hifiasm_ovlp_sink_push(uint32_t q_id, uint32_t t_id,
                            uint32_t q_start, uint32_t q_end,
                            uint32_t t_start, uint32_t t_end,
                            uint32_t n_match, uint32_t block_len,
                            uint32_t shared_seed,
                            uint8_t is_same_strand,
                            const uint16_t *cigar, uint32_t cigar_len,
                            uint32_t cigar_t_start,
                            const uint64_t *chain, uint32_t chain_len);

/* Snapshot the read-name table into the sink. Called once by
 * ha_detect_candidates() just before it tears down R_INF, so the names outlive
 * the read store. name/name_index/total_name_length/total_reads come straight
 * from R_INF. No-op if the sink is not active. */
void hifiasm_ovlp_sink_capture_names(const char *name,
                                     const uint64_t *name_index,
                                     uint64_t total_reads,
                                     uint64_t total_name_length);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HIFIASM_OVERLAPS_INTERNAL_H */
