#ifndef MYLOASM_FFI_H
#define MYLOASM_FFI_H

/*
 * C ABI for myloasm's in-process SNPmer detector (external/myloasm submodule,
 * src/ffi.rs). Lets hifiasm obtain the gated SNPmer set directly instead of
 * parsing snpmers.tsv. Structs here MUST match the #[repr(C)] definitions in
 * ffi.rs field-for-field.
 *
 * Bit order: allele k-mers are 2-bit-packed in MYLOASM order, i.e. the base at
 * read position i is in bits [2i, 2i+1] (first base in the LOW bits), A=0 C=1
 * G=2 T=3. hifiasm's sketch loop rolls k-mers the other way (first base in the
 * HIGH bits), so the loader re-packs each k-mer on ingest before
 * canonicalising. See hetmer_load_from_myloasm().
 */

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint64_t allele0_kmer;
    uint64_t allele1_kmer;
    uint32_t allele0_count;
    uint32_t allele1_count;
} MyloSnpmer;

typedef struct {
    MyloSnpmer *ptr;
    size_t      len;
    size_t      capacity;
    int         k;
} MyloSnpmerSet;

/* Detect SNPmers from n_paths NUL-terminated read-file paths.
 * kmer_size: pass 0 for myloasm's default (21). threads: pass 0 for default.
 * Returns 0 on success (out populated), non-zero on error (out zeroed). */
int  myloasm_detect_snpmers(const char *const *paths, size_t n_paths,
                            int kmer_size, int threads, MyloSnpmerSet *out);

/* Free a set from myloasm_detect_snpmers. Safe on a zeroed set. */
void myloasm_snpmers_free(MyloSnpmerSet *set);

/* ---------------------------------------------------------------------------
 * Per-read marker indexing (open syncmers + SNPmers) for the fake-chain pass.
 * ---------------------------------------------------------------------------
 *
 * myloasm_index_reads() detects SNPmers and then indexes every read with its
 * open-syncmer positions and SNPmer positions (myloasm's get_twin_read_syncmer,
 * k=21, c=11 by default). Each marker carries the read's CANONICAL k-mer key at
 * that position (Kmer48::to_u64: min of forward / reverse-complement under the
 * middle-masked comparison). Because the key is canonical it is strand-
 * invariant: the same genomic locus produces the same key on any read. A
 * consumer can therefore match markers between two reads by key equality alone,
 * with no bit-order or reverse-complement handling.
 *
 * Key layout (k=21): base i at bits [2i, 2i+1]; middle base (index 10) at bits
 * [20,21]. Clipping to k=20 by masking the low 40 bits keeps the middle base,
 * so a SNPmer's two alleles remain distinct after the clip.
 *
 * Structs MUST match the #[repr(C)] definitions in ffi.rs field-for-field.
 */

typedef struct {
    uint32_t pos;   /* 0-based start position on the forward strand */
    uint32_t _pad;  /* explicit padding so key is 8-byte aligned    */
    uint64_t key;   /* canonical k-mer key at pos                    */
} MyloMarker;

typedef struct {
    const char *name;       /* read base id (first token), NOT NUL-terminated */
    size_t      name_len;
    const MyloMarker *syncmers;   /* into shared arena; do not free           */
    size_t      n_syncmers;
    const MyloMarker *snpmers;    /* into shared arena; do not free           */
    size_t      n_snpmers;
} MyloReadMarkers;

/* Opaque-owning result. Only `reads`, `n_reads` and `k` are meant to be read by
 * the host; the remaining fields are private allocation bookkeeping for
 * myloasm_read_index_free(). Layout must still match ffi.rs exactly. */
typedef struct {
    MyloReadMarkers *reads;
    size_t           n_reads;
    size_t           _reads_cap;

    MyloMarker      *_sync_arena;
    size_t           _sync_len;
    size_t           _sync_cap;

    MyloMarker      *_snp_arena;
    size_t           _snp_len;
    size_t           _snp_cap;

    uint8_t         *_name_arena;
    size_t           _name_len;
    size_t           _name_cap;

    int              k;   /* marker k-mer length used (default 21) */
} MyloReadIndex;

/* Detect SNPmers and index all reads. kmer_size/c/threads: pass 0 for defaults
 * (21 / 11 / auto). Returns 0 on success (out populated), non-zero on error
 * (out zeroed). Free with myloasm_read_index_free(). */
int  myloasm_index_reads(const char *const *paths, size_t n_paths,
                         int kmer_size, int c, int threads,
                         MyloReadIndex *out);

/* Free an index from myloasm_index_reads. Safe on a zeroed index. */
void myloasm_read_index_free(MyloReadIndex *idx);

/* ---------------------------------------------------------------------------
 * Anchor chaining (myloasm's DP, exposed for the second pass).
 * ---------------------------------------------------------------------------
 *
 * The host matches markers between a read pair (canonical key equality, m:n
 * within the candidate interval), merges syncmer + SNPmer anchors into one set,
 * and calls myloasm_chain() to run myloasm's own chaining DP over them. Only the
 * chainer is exposed; the caller owns matching and anchor construction.
 *
 * Strand: encode each anchor's query position with the relative strand in bit 31
 * (0 = forward/same-strand, 1 = reverse-complement); pos2 is the raw target
 * position. Within one candidate the strand is uniform (hifiasm already decided
 * it), so all anchors share the same bit.
 *
 * Structs MUST match the #[repr(C)] definitions in ffi.rs field-for-field.
 */

typedef struct {
    uint32_t pos1;   /* query position, relative strand in bit 31 */
    uint32_t pos2;   /* raw target position                       */
    uint32_t slot;   /* opaque caller index (recovers type tag/identity) */
    uint32_t _pad;
} MyloAnchor;

typedef struct {
    uint32_t slot;   /* the input anchor's slot, carried through the DP */
    uint32_t qpos;   /* strand-stripped query position                  */
    uint32_t tpos;   /* raw target position                             */
    uint32_t _pad;
} MyloChainAnchor;

/* Chain a set of anchors with myloasm's DP; writes the single best chain into a
 * caller-provided buffer (ordered along the chain, increasing query position).
 *
 *   anchors/n_anchors : input anchors (need not be sorted).
 *   gap_cost          : c; pass 0 for the default (11).
 *   match_score       : per-anchor score; pass 0 for the default (1).
 *   band              : predecessor-iteration bound = max_mult * 20; pass 0 to
 *                       use n_anchors (unbounded within the set).
 *   min_chain_length  : pass 0 for the default (3).
 *   out/out_cap       : caller buffer (out_cap >= n_anchors is always enough).
 *   out_n             : receives the number of anchors written.
 *   out_score         : receives the chain score.
 *   out_is_reverse    : receives 1 if the best chain is reverse-strand, else 0.
 *
 * Returns 0 on success (including "no chain": *out_n = 0), non-zero on error.
 * No Rust-side allocation is returned, so there is nothing to free. */
int  myloasm_chain(const MyloAnchor *anchors, size_t n_anchors,
                   int gap_cost, int match_score, size_t band,
                   size_t min_chain_length,
                   MyloChainAnchor *out, size_t out_cap,
                   size_t *out_n, int *out_score, int *out_is_reverse);

#ifdef __cplusplus
}
#endif

#endif
