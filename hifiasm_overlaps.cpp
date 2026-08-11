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
