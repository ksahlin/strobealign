// based on https://github.com/lh3/ksw2/blob/master/cli.c

#include <string>
#include <cstdint>
#include <cstring>  // memset
#include <sstream>
#include <algorithm>

#include "ksw2wrap.hpp"
#include "ksw2/ksw2.h"
#include "aligner.hpp"
#include "cigar.hpp"

namespace {

unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
    int i, j;
    a = a < 0? -a : a;
    b = b > 0? -b : b;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j? a : b;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}

void global_aln(
    const std::string& qseq_, const std::string& tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez
) {
    int i, qlen, tlen;
    uint8_t *qseq, *tseq;
    ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
    ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
    ez->n_cigar = 0;
    qlen = qseq_.length();
    tlen = tseq_.length();
    qseq = (uint8_t*)calloc(qlen + 33, 1);
    tseq = (uint8_t*)calloc(tlen + 33, 1);
    for (i = 0; i < qlen; ++i)
        qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
    for (i = 0; i < tlen; ++i)
        tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];

    ksw_extz2_sse(nullptr, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, 0, flag, ez);

    free(qseq);
    free(tseq);
}

}  // namespace


aln_info ksw_extend(const std::string& query, const std::string& ref, int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, bool reverse_cigar) {
    int w = -1; // band width; -1 is inf
    int zdrop = -1; // -1 to disable
    int flag = KSW_EZ_EXTZ_ONLY;
    if (reverse_cigar) {
        flag |= KSW_EZ_REV_CIGAR;
    }

    //int c, i, pair = 1;
    //char *algo = "extd", *s;
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));

    int8_t mat[25];
    ksw_gen_simple_mat(5, mat, match, -mismatch);
    global_aln(query, ref, 5, mat, gap_open, gap_extend, w, zdrop, flag, &ez);

    aln_info info;
    auto cigar = Cigar(ez.cigar, ez.n_cigar).to_eqx(query, ref);
    info.cigar = cigar.to_string();
    info.edit_distance = cigar.edit_distance();
    info.ref_start = 0;
    info.ref_end = ez.max_t + 1;
    info.query_end = ez.max_q + 1;
    info.query_start = 0;
    info.sw_score = ez.max;

    kfree(km, ez.cigar);
    return info;
}
