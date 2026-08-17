#ifndef HETMER_H
#define HETMER_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Het-mer (SNP-discriminating k-mer) support for seeding.
 *
 * A het-mer is a k-mer whose middle base ((k-1)/2) marks a heterozygous site.
 * They are produced externally (e.g. Gene Myers' Phasemer/FastK) and loaded
 * here as a set of canonical 2-bit-packed k-mers. During sketching we scan the
 * RAW read (independently of hifiasm's HPC k=51 minimizer loop) with a k-bit
 * rolling k-mer and force-emit a minimizer wherever a het-mer is recognised,
 * so the chainer anchors on heterozygous positions. This mirrors myloasm's
 * get_twin_read_syncmer, which collects SNPmer positions alongside minimizers
 * in a single read scan.
 *
 * Canonicalisation matches detection: we compare the forward and
 * reverse-complement k-mer with the middle base masked out (the "split" form),
 * pick the smaller as canonical, and store the full canonical k-mer. Both
 * alleles of a site are stored, so recognition is a single set lookup.
 */
typedef struct hetmer_set_s hetmer_set_t;

/* Load het-mers from a Phasemer -Ls listing (lines "  <id>/<base>: <kmer>").
 * Only k-mers of length exactly `k` are used. Returns NULL on open failure.
 * `k` must be odd and <= 31. */
hetmer_set_t *hetmer_load(const char *path, int k);

/* Detect het-mers in-process via myloasm's C ABI (external/myloasm) from the
 * given read files, and build a het-mer set from the gated SNPmers. `kmer_size`
 * may be 0 to use myloasm's default (21). Returns NULL on failure. */
hetmer_set_t *hetmer_load_from_myloasm(const char *const *paths, size_t n_paths,
                                       int kmer_size, int threads);

/* Number of distinct canonical het-mers stored. */
uint64_t hetmer_size(const hetmer_set_t *hs);

/* k used when the set was built. */
int hetmer_k(const hetmer_set_t *hs);

/* Test membership of a canonical (split-min) 2-bit-packed k-mer. */
int hetmer_contains(const hetmer_set_t *hs, uint64_t canon_kmer);

/* Compute the split-canonical form of a FORWARD 2-bit-packed k-mer (hifiasm
 * bit order: first base in the high bits) for this set's k. This is the exact
 * key hetmer_contains() expects and the value stored in the set, so two reads
 * carrying the same allele produce the same canonical value. Use it to both
 * test membership and match SNPmer occurrences between reads by value. */
uint64_t hetmer_canonicalize(const hetmer_set_t *hs, uint64_t fwd_kmer);

void hetmer_destroy(hetmer_set_t *hs);

/* Global set consulted by the sketcher. NULL disables het-mer seeding.
 * Set once at startup after CLI parsing; read-only during sketching. */
extern hetmer_set_t *g_hetmer_set;

/* Diagnostic: total het-mer minimizers force-emitted across all sketch calls
 * (atomic-free; approximate under threads, adequate for a debug counter). */
extern uint64_t g_hetmer_emitted;

#ifdef __cplusplus
}
#endif

#endif
