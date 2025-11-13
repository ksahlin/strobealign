#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <utility>

#include "chain.hpp"
#include "pdqsort/pdqsort.h"
#include "timer.hpp"
#include "logger.hpp"

static Logger& logger = Logger::get();

void add_to_anchors_full(
    std::vector<Anchor>& anchors,
    int query_start,
    int query_end,
    const StrobemerIndex& index,
    size_t position,
    uint prev_query_start,
    int bonus
) {
    int min_diff = std::numeric_limits<int>::max();
    for (const auto hash = index.get_hash(position); index.get_hash(position) == hash; ++position) {
        int ref_start = index.get_strobe1_position(position);
        int ref_end = ref_start + index.strobe2_offset(position) + index.k();
        int diff = std::abs((query_end - query_start) - (ref_end - ref_start));
        if (diff <= min_diff) {
            int ref_idx = index.reference_index(position);
            anchors.push_back({uint(query_start), uint(ref_start), uint(ref_idx), prev_query_start, bonus});
            anchors.push_back({uint(query_end - index.k()), uint(ref_end - index.k()), uint(ref_idx), uint(query_start), 0});
            min_diff = diff;
        }
    }
}

void add_to_anchors_partial(
    std::vector<Anchor>& anchors,
    int query_start,
    const StrobemerIndex& index,
    size_t position,
    uint prev_query_start,
    int bonus
) {
    for (const auto hash = index.get_main_hash(position); index.get_main_hash(position) == hash; ++position) {
        auto [ref_start, ref_end] = index.strobe_extent_partial(position);
        int ref_idx = index.reference_index(position);
        anchors.push_back({uint(query_start), uint(ref_start), uint(ref_idx), prev_query_start, bonus});
    }
}

void add_hits_to_anchors(
    const std::vector<Hit>& hits,
    const StrobemerIndex& index,
    std::vector<Anchor>& anchors
) {
    if (hits.empty()) {
        return;
    }

    // Split hits up into their constituent strobes
    struct Strobe {
        size_t query_start;
        size_t hit_index; // the original hit this was split from

        bool operator<(const Strobe& other) const {
            return query_start < other.query_start;
        }
    };

    std::vector<Strobe> strobes;
    for (size_t i = 0; i < hits.size(); ++i) {
        const Hit& hit = hits[i];
        if (hit.is_partial) {
            strobes.push_back(Strobe{hit.query_start, i});
        } else {
            // Full hits result in two constituent hits
            strobes.push_back(Strobe{hit.query_start, i});
            strobes.push_back(Strobe{hit.query_end - index.k(), i});
        }
    }
    std::sort(strobes.begin(), strobes.end());

    for (auto strobe : strobes) {
        std::cerr << " - Strobe ";
        if (hits[strobe.hit_index].is_partial) {
            std::cerr << "p";
        } else {
            if (hits[strobe.hit_index].query_start == strobe.query_start) {
                std::cerr << "1";
            } else {
                std::cerr << "2";

            }
        }
        std::cerr << " at query_start=" << strobe.query_start << " from hit " << strobe.hit_index;
        std::cerr << '\n';
    }


    // Map hit indices back to strobe indices
    std::vector<size_t> hits_to_strobes1;
    std::vector<size_t> hits_to_strobes2;
    hits_to_strobes1.assign(hits.size(), 0);
    hits_to_strobes2.assign(hits.size(), 0);

    for (size_t i = 0; i < strobes.size(); ++i) {
        size_t hit_index = strobes[i].hit_index;
        if (hits[hit_index].is_partial) {
            hits_to_strobes1[hit_index] = i;
            // TODO remove this if if below assertion holds
            assert(hits[hit_index].query_start == strobes[i].query_start);
        } else {
            // full hit
            if (hits[hit_index].query_start == strobes[i].query_start) {
                hits_to_strobes1[hit_index] = i;
            } else {
                hits_to_strobes2[hit_index] = i;
            }
        }
    }

    std::vector<uint> prev_query_starts;
    std::vector<uint> bonuses;

    prev_query_starts.assign(strobes.size(), 0);
    bonuses.assign(strobes.size(), 0);
    uint filtered_length = 0;
    size_t last_unfiltered = 0;
    for (size_t i = 1; i < strobes.size(); i++) {
        Strobe& strobe = strobes[i];
        if (hits[strobe.hit_index].is_filtered) {
            filtered_length += strobe.query_start + index.k() - std::max(strobe.query_start, strobes[i-1].query_start);
        } else {
            prev_query_starts[i] = strobes[last_unfiltered].query_start;
            bonuses[i] = filtered_length;
            filtered_length = 0;
            last_unfiltered = i;
        }
    }

    for (size_t hit_index = 0; hit_index < hits.size(); hit_index++) {
        const Hit& hit = hits[hit_index];
        if (hit.is_filtered) {
            continue;
        }
        if (hit.is_partial) {
            // add_to_anchors_partial(anchors, hit.query_start, index, hit.position, prev_unfiltered->query_start, filtered_length);
            int query_start = hit.query_start;
            size_t position = hit.position;
            uint prev_query_start = prev_query_starts[hits_to_strobes1[hit_index]];
            int bonus = bonuses[hits_to_strobes1[hit_index]];

            for (const auto hash = index.get_main_hash(position); index.get_main_hash(position) == hash; ++position) {
                auto [ref_start, ref_end] = index.strobe_extent_partial(position);
                int ref_idx = index.reference_index(position);
                logger.trace() << "Adding p at query_start=" << query_start << " bonus=" << bonus << " prev_query_start=" << prev_query_start << '\n';

                anchors.push_back({uint(query_start), uint(ref_start), uint(ref_idx), prev_query_start, bonus});
            }
        } else {
            // add_to_anchors_full(anchors, hit.query_start, hit.query_end, index, hit.position, prev_unfiltered->query_start, filtered_length);
            int query_start = hit.query_start;
            int query_end = hit.query_end;
            size_t position = hit.position;

            uint prev_query_start1 = prev_query_starts[hits_to_strobes1[hit_index]];
            uint prev_query_start2 = prev_query_starts[hits_to_strobes2[hit_index]];
            int bonus1 = bonuses[hits_to_strobes1[hit_index]];
            int bonus2 = bonuses[hits_to_strobes2[hit_index]];
            logger.trace() << "Adding 1 at query_start=" << query_start << " bonus=" << bonus1 << " prev_query_start=" << prev_query_start1 << '\n';
            logger.trace() << "Adding 2 at query_start=" << query_end - index.k() << " bonus=" << bonus2 << " prev_query_start=" << prev_query_start2 << '\n';
            assert(prev_query_start1 < query_start);
            assert(prev_query_start2 < query_end - index.k());

            int min_diff = std::numeric_limits<int>::max();
            for (const auto hash = index.get_hash(position); index.get_hash(position) == hash; ++position) {
                int ref_start = index.get_strobe1_position(position);
                int ref_end = ref_start + index.strobe2_offset(position) + index.k();
                int diff = std::abs((query_end - query_start) - (ref_end - ref_start));
                if (diff <= min_diff) {
                    int ref_idx = index.reference_index(position);
                    anchors.push_back({uint(query_start), uint(ref_start), uint(ref_idx), prev_query_start1, bonus1});
                    anchors.push_back({uint(query_end - index.k()), uint(ref_end - index.k()), uint(ref_idx), prev_query_start2, bonus2});
                    min_diff = diff;
                }
            }
        }
    }
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

            float new_score = dp[j] + score;

            // If we are chaining to a direct predecessor on the query and are
            // are on the same diagonal, add a bonus to the score that tries to
            // compensate for skipped repetitive hits between the two anchors
            if (dq == dr && ai.prev_query_start == aj.query_start) {
                //logger.trace() << "Adding bonus " << ai.skipped_filtered_bonus << "\n";
                new_score += ai.skipped_filtered_bonus;
            }

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

        chains.push_back(Nam{
            int(chains.size()),        // nam_id
            int(first.query_start),    // query_start
            int(last.query_start + k), // query_end
            -1,                        // query_prev_match_startpos
            int(first.ref_start),      // ref_start
            int(last.ref_start + k),   // ref_end
            -1,                        // ref_prev_match_startpos
            c,                         // n_matches
            int(last.ref_id),          // ref_id
            score + c * chaining_params.matches_weight, // score
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
    std::array<HitsDetails, 2> hits_details;
    for (int is_revcomp : {0, 1}) {
        bool sorting_needed1;
        std::tie(hits_details[is_revcomp], sorting_needed1, hits[is_revcomp]) =
            find_hits(query_randstrobes[is_revcomp], index, map_param.mcs_strategy, map_param.rescue_threshold);
        details.hits += hits_details[is_revcomp];
    }
    statistics.time_hit_finding += hits_timer.duration();

    std::vector<Nam> chains;

    // Runtime heuristic: If one orientation appears to have many more hits
    // than the other, we assume it is the correct one and do not check the
    // other.
    std::vector<int> orientations;
    if (hits_details[0].is_better_than(hits_details[1])) {
        orientations.push_back(0);
    } else if (hits_details[1].is_better_than(hits_details[0])) {
        orientations.push_back(1);
    } else {
        orientations.push_back(0);
        orientations.push_back(1);
    }

    for (int is_revcomp : orientations) {
        float best_score = 0.0f;
        std::vector<Anchor> anchors;
        std::vector<float> dp;
        std::vector<int> predecessors;

        Timer hits_timer;
        add_hits_to_anchors(hits[is_revcomp], index, anchors);
        statistics.time_hit_finding += hits_timer.duration();
        statistics.n_anchors += anchors.size();
        Timer chaining_timer;
        logger.trace() << "Chaining " << anchors.size() << " anchors\n";
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

std::ostream& operator<<(std::ostream& os, const Anchor& anchor) {
    os << "Anchor(ref_id=" << anchor.ref_id << ", ref_start=" << anchor.ref_start << ", query_start=" << anchor.query_start << ")";
    return os;
}
