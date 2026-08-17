/*
 * Experiment 1 instrumentation -- NOT for upstream, NOT committed to the fork.
 *
 * Records the length of every exact-match query issued against the GLOBAL graph
 * FM index. The question it answers: HISAT2 searches the global index with
 * ftabLoHi followed by repeated mapLF/mapGLF, i.e. a maximal exact match whose
 * length is bounded only by the read (maxHitLen defaults to INDEX_MAX). If an
 * order-k index (GCSA2-style) is to replace the exact one, queries longer than k
 * can return false positives -- so the design hinges on this distribution.
 *
 * Header-only so no Makefile change is needed. Enabled by setting
 * HISAT2_QLEN_OUT to an output path; otherwise every call is a predictable
 * branch on a cached bool.
 */
#ifndef QLEN_PROBE_H_
#define QLEN_PROBE_H_

#include <cstdio>
#include <cstdlib>

#define QLEN_MAX 2048

inline unsigned long long* qlen_hist() {          // global graph FM index
    static unsigned long long h[QLEN_MAX + 2] = {0};
    return h;
}

inline unsigned long long* qlen_hist_local() {    // 57,344 bp local indexes
    static unsigned long long h[QLEN_MAX + 2] = {0};
    return h;
}

inline const char* qlen_out_path() {
    static const char* p = getenv("HISAT2_QLEN_OUT");
    return p;
}

struct QlenDumper {
    ~QlenDumper() {
        const char* p = qlen_out_path();
        if(!p) return;
        FILE* f = fopen(p, "w");
        if(!f) return;
        unsigned long long* h = qlen_hist();
        unsigned long long* l = qlen_hist_local();
        fprintf(f, "# scope\tlength\tcount\n");
        for(int i = 0; i <= QLEN_MAX + 1; i++) {
            if(h[i]) fprintf(f, "global\t%d\t%llu\n", i, h[i]);
        }
        for(int i = 0; i <= QLEN_MAX + 1; i++) {
            if(l[i]) fprintf(f, "local\t%d\t%llu\n", i, l[i]);
        }
        fclose(f);
    }
};

inline void qlen_record(unsigned long v) {
    if(!qlen_out_path()) return;
    static QlenDumper dumper;   // constructed on first use, dumps at exit
    (void)dumper;
    if(v > QLEN_MAX) v = QLEN_MAX + 1;
    qlen_hist()[v]++;
}

inline void qlen_record_local(unsigned long v) {
    if(!qlen_out_path()) return;
    static QlenDumper dumper;
    (void)dumper;
    if(v > QLEN_MAX) v = QLEN_MAX + 1;
    qlen_hist_local()[v]++;
}

#endif /*QLEN_PROBE_H_*/
