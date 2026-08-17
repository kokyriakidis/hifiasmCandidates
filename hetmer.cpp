#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "hetmer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "khashl.h"
#include "htab.h"          /* seq_nt4_table */
#include "myloasm_ffi.h"   /* myloasm_detect_snpmers */

/* uint64 hash set (open-addressing, from khashl). */
KHASHL_SET_INIT(static klib_unused, hetset_t, hetset, uint64_t,
                kh_hash_uint64, kh_eq_generic)

struct hetmer_set_s {
    hetset_t *h;
    int       k;          /* k-mer length (odd, <= 31) */
    uint64_t  mask;       /* 2*k low bits set */
    uint64_t  split_mask; /* mask with the middle 2 bits cleared */
};

hetmer_set_t *g_hetmer_set = NULL;
uint64_t      g_hetmer_emitted = 0;

int hetmer_k(const hetmer_set_t *hs) { return hs ? hs->k : 0; }
uint64_t hetmer_size(const hetmer_set_t *hs) {
    return hs ? (uint64_t)kh_size(hs->h) : 0;
}

/* Canonicalise a forward 2-bit k-mer: compute the reverse complement, then
 * pick whichever of (fwd, rev) is smaller when the middle base is masked out
 * (the "split" comparison). Returns the full canonical k-mer (middle base
 * intact). This must match the recognition rule in the sketch loop exactly. */
static inline uint64_t canon_split(uint64_t fwd, int k, uint64_t mask,
                                   uint64_t split_mask, uint64_t *rev_out) {
    uint64_t rev = 0, x = fwd;
    for (int i = 0; i < k; ++i) {
        rev = (rev << 2) | (3ULL - (x & 3ULL)); /* complement + reverse */
        x >>= 2;
    }
    rev &= mask;
    if (rev_out) *rev_out = rev;
    uint64_t sf = fwd & split_mask;
    uint64_t sr = rev & split_mask;
    return (sf < sr) ? fwd : rev;
}

/* Insert a forward 2-bit k-mer (hifiasm bit order: first base in high bits)
 * into the set as its split-canonical form. Returns 1 if newly added. */
static inline int insert_fwd(hetmer_set_t *hs, uint64_t fwd) {
    uint64_t rev;
    uint64_t canon = canon_split(fwd & hs->mask, hs->k, hs->mask,
                                 hs->split_mask, &rev);
    int absent;
    hetset_put(hs->h, canon, &absent);
    return absent;
}

/* Re-pack a k-mer from myloasm bit order (base i at bits [2i,2i+1], first base
 * low) to hifiasm bit order (first base high). This reverses the base order in
 * 2-bit units, which is exactly the reversed-base-order convention seen in
 * myloasm's TSV. */
static inline uint64_t mylo_to_hifi(uint64_t km, int k) {
    uint64_t out = 0;
    for (int i = 0; i < k; ++i) {
        out = (out << 2) | (km & 3ULL);
        km >>= 2;
    }
    return out;
}

static hetmer_set_t *hetmer_new(int k) {
    hetmer_set_t *hs = (hetmer_set_t*)calloc(1, sizeof(hetmer_set_t));
    hs->h = hetset_init();
    hs->k = k;
    hs->mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    int mid = (k - 1) / 2;
    hs->split_mask = hs->mask & ~(3ULL << (2 * mid));
    return hs;
}

hetmer_set_t *hetmer_load_from_myloasm(const char *const *paths, size_t n_paths,
                                       int kmer_size, int threads) {
    MyloSnpmerSet set;
    int rc = myloasm_detect_snpmers(paths, n_paths, kmer_size, threads, &set);
    if (rc != 0) {
        fprintf(stderr, "[hetmer] myloasm_detect_snpmers failed (rc=%d)\n", rc);
        return NULL;
    }
    int k = set.k;
    if (k <= 0 || k > 31 || (k & 1) == 0) {
        fprintf(stderr, "[hetmer] myloasm returned invalid k=%d\n", k);
        myloasm_snpmers_free(&set);
        return NULL;
    }
    hetmer_set_t *hs = hetmer_new(k);
    uint64_t added = 0;
    for (size_t i = 0; i < set.len; ++i) {
        uint64_t a0 = mylo_to_hifi(set.ptr[i].allele0_kmer, k);
        uint64_t a1 = mylo_to_hifi(set.ptr[i].allele1_kmer, k);
        added += insert_fwd(hs, a0);
        added += insert_fwd(hs, a1);
    }
    fprintf(stderr,
            "[hetmer] myloasm detected %zu SNPmer sites -> %llu distinct"
            " canonical %d-mers\n",
            set.len, (unsigned long long)added, k);
    myloasm_snpmers_free(&set);
    return hs;
}

hetmer_set_t *hetmer_load(const char *path, int k) {
    if (!path || k <= 0 || k > 31 || (k & 1) == 0) {
        fprintf(stderr, "[hetmer] invalid k=%d (must be odd, 1..31)\n", k);
        return NULL;
    }
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[hetmer] cannot open %s\n", path);
        return NULL;
    }

    hetmer_set_t *hs = hetmer_new(k);

    char *line = NULL;
    size_t cap = 0;
    ssize_t nread;
    uint64_t n_lines = 0, n_added = 0, n_bad = 0;

    while ((nread = getline(&line, &cap, fp)) != -1) {
        /* Phasemer -Ls line: "   <id>/<base>: <kmer>"; also accept a bare kmer
         * or "<kmer> <count>". Extract the last whitespace-delimited token that
         * looks like a k-mer of length k. */
        char *colon = strrchr(line, ':');
        char *tok = colon ? colon + 1 : line;
        /* trim leading space */
        while (*tok == ' ' || *tok == '\t') ++tok;
        /* cut at first whitespace/newline */
        char *e = tok;
        while (*e && *e != ' ' && *e != '\t' && *e != '\n' && *e != '\r') ++e;
        int len = (int)(e - tok);
        if (len != k) { ++n_bad; continue; }

        uint64_t fwd = 0;
        int ok = 1;
        for (int i = 0; i < k; ++i) {
            int c = seq_nt4_table[(unsigned char)tok[i]];
            if (c >= 4) { ok = 0; break; }
            fwd = (fwd << 2) | (uint64_t)c;
        }
        if (!ok) { ++n_bad; continue; }

        if (insert_fwd(hs, fwd)) ++n_added;
        ++n_lines;
    }
    free(line);
    fclose(fp);

    fprintf(stderr,
            "[hetmer] loaded %s: %llu lines, %llu distinct canonical %d-mers"
            " (%llu skipped)\n",
            path, (unsigned long long)n_lines,
            (unsigned long long)n_added, k, (unsigned long long)n_bad);
    return hs;
}

uint64_t hetmer_canonicalize(const hetmer_set_t *hs, uint64_t fwd_kmer) {
    if (!hs) return fwd_kmer;
    return canon_split(fwd_kmer & hs->mask, hs->k, hs->mask, hs->split_mask,
                       NULL);
}

int hetmer_contains(const hetmer_set_t *hs, uint64_t canon_kmer) {
    if (!hs) return 0;
    khint_t itr = hetset_get(hs->h, canon_kmer);
    return itr != kh_end(hs->h);
}

void hetmer_destroy(hetmer_set_t *hs) {
    if (!hs) return;
    hetset_destroy(hs->h);
    free(hs);
}
