#ifndef STROBEALIGN_CHAIN_HPP
#define STROBEALIGN_CHAIN_HPP
#include <algorithm>
#include <vector>
#include "index.hpp"
#include "mappingparameters.hpp"

static constexpr int N_PRECOMPUTED = 1024;

struct Anchor {
    uint query_start;
    uint ref_start;
    uint ref_id;

    bool operator<(const Anchor& other) const {
        return (ref_id < other.ref_id) || (ref_id == other.ref_id && ref_start < other.ref_start) || (ref_id == other.ref_id && ref_start == other.ref_start && query_start < other.query_start);
    }

    bool operator==(const Anchor& other) const {
        return (ref_id == other.ref_id) && ref_start == other.ref_start && query_start == other.query_start;
    }
};

struct Chain {
    uint id;
    uint ref_id;
    float score;
    std::vector<Anchor> anchors;
    uint query_start;
    uint query_end;
    uint ref_start;
    uint ref_end;
    bool is_revcomp;

    int ref_span() const {
        return ref_end - ref_start;
    }

    int query_span() const {
        return query_end - query_start;
    }

    int projected_ref_start() const {
        return std::max(0, int(ref_start) - int(query_start));
    }
};

struct Chainer {
    Chainer(ChainingParameters chaining_params, int k)
        : k(k)
        , chaining_params(chaining_params)
    {
        precomputed_scores[0] = 0;
        for (size_t i = 1; i < N_PRECOMPUTED; i++) {
            precomputed_scores[i] = compute_score_uncached(i, i);
        }
    }

    std::vector<Chain> get_chains(
        const std::array<std::vector<QueryRandstrobe>, 2>& query_randstrobes,
        const StrobemerIndex& index,
        AlignmentStatistics& statistics,
        Details& details,
        const MappingParameters& map_param
    ) const;

  private:
    const int k;
    const ChainingParameters chaining_params;
    float precomputed_scores[N_PRECOMPUTED];
    float collinear_chaining(
        const std::vector<Anchor>& anchors,
        std::vector<float>& dp,
        std::vector<int>& predecessors
    ) const;

    float compute_score(const int dq, const int dr) const;
    float compute_score_uncached(const int dq, const int dr) const;
};

/**
 * Fast log2 function for x >= 1.
 * Copied from https://github.com/lh3/minimap2/blob/master/mmpriv.h
 */
static inline float mg_log2(float x) {
    union { float f; uint32_t i; } z = { x };
    float log_2 = ((z.i >> 23) & 255) - 128;
    z.i &= ~(255 << 23);
    z.i += 127 << 23;
    log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
    return log_2;
}

std::ostream& operator<<(std::ostream& os, const Chain& chain);
std::ostream& operator<<(std::ostream& os, const Anchor& anchor);

#endif
