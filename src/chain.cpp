#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <utility>

#include "chain.hpp"
#include "pdqsort/pdqsort.h"
#include "timer.hpp"

void add_to_anchors_full(
    std::vector<Anchor>& anchors,
    int query_start,
    int query_end,
    const StrobemerIndex& index,
    size_t position
) {
    int min_diff = std::numeric_limits<int>::max();
    for (const auto hash = index.get_hash(position); index.get_hash(position) == hash; ++position) {
        int ref_start = index.get_strobe1_position(position);
        int ref_end = ref_start + index.strobe2_offset(position) + index.k();
        int diff = std::abs((query_end - query_start) - (ref_end - ref_start));
        if (diff <= min_diff) {
            int ref_idx = index.reference_index(position);
            anchors.push_back({uint(query_start), uint(ref_start), uint(ref_idx)});
            anchors.push_back({uint(query_end - index.k()), uint(ref_end - index.k()), uint(ref_idx)});
            min_diff = diff;
        }
    }
}

void add_to_anchors_partial(
    std::vector<Anchor>& anchors,
    int query_start,
    const StrobemerIndex& index,
    size_t position
) {
    for (const auto hash = index.get_main_hash(position); index.get_main_hash(position) == hash; ++position) {
        auto [ref_start, ref_end] = index.strobe_extent_partial(position);
        int ref_idx = index.reference_index(position);
        anchors.push_back({uint(query_start), uint(ref_start), uint(ref_idx)});
    }
}

void add_hits_to_anchors(
    const std::vector<Hit>& hits,
    const StrobemerIndex& index,
    std::vector<Anchor>& anchors
) {
    for (const Hit& hit : hits) {
        if (hit.is_partial) {
            add_to_anchors_partial(anchors, hit.query_start, index, hit.position);
        } else {
            add_to_anchors_full(anchors, hit.query_start, hit.query_end, index, hit.position);
        }
    }
}

std::tuple<int, int> find_anchors_rescue(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    McsStrategy mcs_strategy,
    std::vector<Anchor>& anchors
) {
    struct RescueHit {
        size_t position;
        unsigned int count;
        unsigned int query_start;
        unsigned int query_end;
        bool is_partial;

        bool operator<(const RescueHit& rhs) const {
            return std::tie(count, query_start, query_end) <
                   std::tie(rhs.count, rhs.query_start, rhs.query_end);
        }
    };

    int n_hits = 0;
    int partial_hits = 0;
    std::vector<RescueHit> rescue_hits;
    for (auto& qr : query_randstrobes) {
        size_t position = index.find_full(qr.hash);
        if (position != index.end()) {
            unsigned int count = index.get_count_full(position);
            size_t position_revcomp = index.find_full(qr.hash_revcomp);
            if (position_revcomp != index.end()) {
                count += index.get_count_full(position_revcomp);
            }
            rescue_hits.push_back({position, count, qr.start, qr.end, false});
        } else if (mcs_strategy == McsStrategy::Always) {
            size_t partial_pos = index.find_partial(qr.hash);
            if (partial_pos != index.end()) {
                unsigned int partial_count = index.get_count_partial(partial_pos);
                size_t position_revcomp = index.find_partial(qr.hash_revcomp);
                if (position_revcomp != index.end()) {
                    partial_count += index.get_count_partial(position_revcomp);
                }
                rescue_hits.push_back({partial_pos, partial_count, qr.start, qr.start + index.k(), true});
                partial_hits++;
            }
        }
    }

    std::sort(rescue_hits.begin(), rescue_hits.end());

    int cnt = 0;
    for (auto& rh : rescue_hits) {
        if ((rh.count > rescue_cutoff && cnt >= 5) || rh.count > 1000) {
            break;
        }
        if (rh.is_partial) {
            partial_hits++;
            add_to_anchors_partial(anchors, rh.query_start, index, rh.position);
        } else {
            add_to_anchors_full(anchors, rh.query_start, rh.query_end, index, rh.position);
        }
        cnt++;
        n_hits++;
    }

    return {n_hits, partial_hits};
}

/**
 * @brief Compute the chaining score between two anchors based on their distance.
 *
 * This function calculates the score contribution when chaining two anchors,
 * considering their relative positions on the query and reference sequences.
 * It penalizes large gaps and diagonal differences to encourage collinear chains.
 *
 * @param dq Difference in query start positions between anchor i and anchor j
 *           (i.e., ai.query_start - aj.query_start). Must be > 0.
 * @param dr Difference in reference start positions between anchor i and anchor j
 *           (i.e., ai.ref_start - aj.ref_start). Must be > 0.
 * @param k  Length of the k-mer used to form the anchor.
 * @param chaining_params Parameters controlling penalties:
 *        - diag_diff_penalty: Multiplier for the absolute difference |dr - dq|, penalizing non-diagonal moves.
 *        - gap_length_penalty: Multiplier for min(dq, dr), penalizing longer gaps.
 *
 * @return A float representing the score for chaining the two anchors.
 */
float Chainer::compute_score(const int dq, const int dr) const {
    // dq == dr is usually the most common case
    if (dq == dr && dq < N_PRECOMPUTED) {
        return precomputed_scores[dq];
    } else {
        return compute_score_uncached(dq, dr);
    }
}

float Chainer::compute_score_uncached(const int dq, const int dr) const {
    const int dd = std::abs(dr - dq);
    const int dg = std::min(dq, dr);
    float score = std::min(k, dg);

    const float lin_penalty = chaining_params.diag_diff_penalty * dd + chaining_params.gap_length_penalty * dg;
    const float log_penalty = dd >= 1 ? mg_log2(dd + 1) : 0.0f;
    score -= lin_penalty + 0.5 * log_penalty;

    return score;
}

float Chainer::collinear_chaining(
    const std::vector<Anchor>& anchors,
    std::vector<float>& dp,
    std::vector<int>& predecessors
) const {
    const size_t n = anchors.size();
    if (n == 0) {
        return 0;
    }

    dp.assign(n, k);
    predecessors.assign(n, -1);
    float best_score = 0;

    for (size_t i = 0; i < n; ++i) {
        const int lookup_end = std::max(0, static_cast<int>(i) - chaining_params.max_lookback);

        for (int j = i - 1; j >= lookup_end; --j) {
            const Anchor& ai = anchors[i];
            const Anchor& aj = anchors[j];

            if (ai.ref_id != aj.ref_id) {
                break;
            }

            const int dq = ai.query_start - aj.query_start;
            const int dr = ai.ref_start - aj.ref_start;

            if (dr >= chaining_params.max_ref_gap) {
                break;
            }
            if (dq <= 0 || dr <= 0) {
                continue;
            }

            const float score = compute_score(dq, dr);

            const float new_score = dp[j] + score;
            if (new_score > dp[i]) {
                dp[i] = new_score;
                predecessors[i] = j;

                // Runtime heuristic: If the predecessor is on the same diagonal,
                // assume that it is the best one and skip the remaining ones.
                if (dq == dr) {
                    break;
                }
            }
        }
        if (dp[i] > best_score) {
            best_score = dp[i];
        }
    }
    return best_score;
}

void extract_chains_from_dp(
    const std::vector<Anchor>& anchors,
    const std::vector<float>& dp,
    const std::vector<int>& predecessors,
    float best_score,
    const int k,
    bool is_revcomp,
    std::vector<Nam>& chains,
    const ChainingParameters& chaining_params
) {
    const size_t n = anchors.size();
    const float valid_score = best_score * chaining_params.valid_score_threshold;
    std::vector<bool> used(n, false);

    std::vector<std::pair<int, float>> candidates;
    for (size_t i = 0; i < n; ++i) {
        if (dp[i] >= valid_score) {
            candidates.push_back(std::make_pair(i, dp[i]));
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.second > b.second; });

    for (const auto& [i, score] : candidates) {
        if (used[i]) {
            continue;
        }

        int j = i;
        int c = 1;
        bool overlaps = false;

        while (predecessors[j] >= 0) {
            j = predecessors[j];
            if (used[j]) {
                overlaps = true;
                break;
            }
            used[j] = true;
            c++;
        }

        if (overlaps) {
            continue;
        }

        const Anchor& first = anchors[j];
        const Anchor& last = anchors[i];

        chains.push_back(
            Nam{int(chains.size()),
            int(first.query_start),
            int(last.query_start + k),
            -1,
            int(first.ref_start),
            int(last.ref_start + k),
            -1,
            c,
            int(last.ref_id),
            score + c * chaining_params.matches_weight,
            is_revcomp
            }
        );
    }
}

std::vector<Nam> Chainer::get_chains(
    const std::array<std::vector<QueryRandstrobe>, 2>& query_randstrobes,
    const StrobemerIndex& index,
    AlignmentStatistics& statistics,
    Details& details,
    const MappingParameters& map_param
) const {
    Timer hits_timer;

    std::array<std::vector<Hit>, 2> hits;
    for (int is_revcomp : {0, 1}) {
        bool sorting_needed1;
        HitsDetails hits_details1;
        std::tie(hits_details1, sorting_needed1, hits[is_revcomp]) =
            find_hits(query_randstrobes[is_revcomp], index, map_param.mcs_strategy);
        details.hits += hits_details1;
    }
    uint total_hits = details.hits.total_hits();
    int nonrepetitive_hits = hits[0].size() + hits[1].size();
    float nonrepetitive_fraction = total_hits > 0 ? ((float) nonrepetitive_hits) / ((float) total_hits) : 1.0;
    statistics.time_hit_finding += hits_timer.duration();

    std::vector<Nam> chains;

    for (int is_revcomp : {0, 1}) {
        float best_score = 0.0f;
        std::vector<Anchor> anchors;
        std::vector<float> dp;
        std::vector<int> predecessors;

        // Rescue if requested and needed
        if (map_param.rescue_level > 1 && (nonrepetitive_hits == 0 || nonrepetitive_fraction < 0.7)) {
            Timer rescue_timer;
            auto [n_hits, n_partial_hits] = find_anchors_rescue(
                query_randstrobes[is_revcomp], index, map_param.rescue_cutoff, map_param.mcs_strategy, anchors
            );
            statistics.n_rescue_hits += n_hits;
            statistics.n_rescue_partial_hits += n_partial_hits;
            details.rescue_nams += chains.size();
            details.nam_rescue = true;
            statistics.tot_time_rescue += rescue_timer.duration();
        } else {
            Timer hits_timer;
            add_hits_to_anchors(hits[is_revcomp], index, anchors);
            statistics.time_hit_finding += hits_timer.duration();
        }
        Timer chaining_timer;
        pdqsort(anchors.begin(), anchors.end());
        anchors.erase(
            std::unique(anchors.begin(), anchors.end()), anchors.end()
        );
        float score = collinear_chaining(anchors, dp, predecessors);
        best_score = score;

        extract_chains_from_dp(
            anchors, dp, predecessors, best_score,
            index.k(), is_revcomp, chains, map_param.chaining_params
        );
        statistics.time_chaining += chaining_timer.duration();
    }
    details.nams += chains.size();

    return chains;
}
