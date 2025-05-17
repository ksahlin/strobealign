#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <set>
#include <sstream>
#include <utility>

#include "aln.hpp"
#include "chain.hpp"
#include "randstrobes.hpp"
#include "robin_hood.h"

inline void add_anchor(std::vector<Anchor>& anchors, int query_start, int ref_start, int ref_index) {
    anchors.push_back({query_start, ref_start, ref_index});
}

inline void add_to_anchors_vector_full(
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
            add_anchor(anchors, query_start, ref_start, ref_idx);
            add_anchor(anchors, query_end - index.k(), ref_end - index.k(), ref_idx);
            min_diff = diff;
        }
    }
}

inline void add_to_anchors_vector_partial(
    std::vector<Anchor>& anchors,
    int query_start,
    const StrobemerIndex& index,
    size_t position
) {
    for (const auto hash = index.get_main_hash(position); index.get_main_hash(position) == hash; ++position) {
        auto [ref_start, ref_end] = index.strobe_extent_partial(position);
        int ref_idx = index.reference_index(position);
        add_anchor(anchors, query_start, ref_start, ref_idx);
    }
}

void hits_to_anchors_vector(const std::vector<Hit>& hits, const StrobemerIndex& index, std::vector<Anchor>& anchors) {
    for (const Hit& hit : hits) {
        if (hit.is_partial) {
            add_to_anchors_vector_partial(anchors, hit.query_start, index, hit.position);
        } else {
            add_to_anchors_vector_full(anchors, hit.query_start, hit.query_end, index, hit.position);
        }
    }
    std::sort(anchors.begin(), anchors.end());
}

std::tuple<int, int> find_anchors_rescue(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    bool use_mcs,
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
        } else if (use_mcs) {
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
            add_to_anchors_vector_partial(anchors, rh.query_start, index, rh.position);
        } else {
            add_to_anchors_vector_full(anchors, rh.query_start, rh.query_end, index, rh.position);
        }
        cnt++;
        n_hits++;
    }

    std::sort(anchors.begin(), anchors.end());
    return {n_hits, partial_hits};
}

static inline float
compute_score(const Anchor& ai, const Anchor& aj, const int k, const ChainingPrameters& ch_params) {
    if (ai.ref_id != aj.ref_id)
        return FLT_MIN; 

    int dq = ai.query_start - aj.query_start;
    int dr = ai.ref_start - aj.ref_start;

    if (dq <= 0 || dr <= 0)
        return FLT_MIN;

    int dd = std::abs(dr - dq);
    int dg = std::min(dq, dr);
    float score = std::min(k, dg);

    float penalty = ch_params.gd * dd + ch_params.gl * dg;
    score -= penalty;

    return score;
}

void collinear_chaining(
    std::vector<Anchor>& anchors,
    const int k,
    bool is_revcomp,
    std::vector<Nam>& chains,
    const ChainingPrameters& ch_params
) {
    const int n = anchors.size();
    if (n == 0) {
        return;
    }

    std::vector<float> dp(n, k);
    std::vector<int> backtrack(n, -1);

    int best_score = 0;

    for (int i = 0; i < n; ++i) {
        const int lookup_end = std::max(0, i - ch_params.h);

        for (int j = i - 1; j >= lookup_end; --j) {
            float score = compute_score(anchors[i], anchors[j], k, ch_params);
            if (score == FLT_MIN)
                continue;

            int new_score = dp[j] + score;
            if (new_score > dp[i]) {
                dp[i] = new_score;
                backtrack[i] = j;
            }
        }
        if (dp[i] > best_score) {
            best_score = dp[i];
        }
    }

    const float valid_score = best_score * ch_params.vp;
    std::vector<bool> used(n, false);

    for (int i = n - 1; i >= 0; --i) {
        if (dp[i] < valid_score || used[i]) {
            continue;
        }

        int j = i;
        int c = 0;

        while (backtrack[j] >= 0) {
            j = backtrack[j];
            used[j] = true;
            c++;
        }

        const Anchor& first = anchors[j];
        const Anchor& last = anchors[i];
        const float score = dp[i];

        chains.push_back(Nam{
            int(chains.size()), first.query_start, last.query_start + k, -1, first.ref_start,
            last.ref_start + k, -1, c, last.ref_id, score, is_revcomp

        });
    }
}

template <typename T>
bool by_score(const T& a, const T& b) {
    return a.score > b.score;
}

void shuffle_top_chains(std::vector<Nam>& chains, std::minstd_rand& random_engine) {
    if (chains.empty()) {
        return;
    }
    auto best_score = chains[0].score;
    auto it = std::find_if(chains.begin(), chains.end(), [&](const Nam& chain) {
        return chain.score != best_score;
    });
    if (it > chains.begin() + 1) {
        std::shuffle(chains.begin(), it, random_engine);
    }
}

std::vector<Nam> get_chains(
    const klibpp::KSeq& record,
    const StrobemerIndex& index,
    // AlignmentStatistics& statistics,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    std::minstd_rand& random_engine
) {
    // Compute randstrobes
    // Timer strobe_timer;
    auto query_randstrobes = randstrobes_query(record.seq, index_parameters);
    // statistics.n_randstrobes += query_randstrobes[0].size() + query_randstrobes[1].size();
    // statistics.tot_construct_strobemers += strobe_timer.duration();

    // Find NAMs
    // Timer nam_timer;

    int total_hits = 0;
    // int partial_hits = 0;
    std::array<std::vector<Hit>, 2> hits;
    for (int is_revcomp : {0, 1}) {
        int total_hits1, partial_hits1;
        bool sorting_needed1;
        std::tie(total_hits1, partial_hits1, sorting_needed1, hits[is_revcomp]) =
            find_hits(query_randstrobes[is_revcomp], index, map_param.use_mcs);
        total_hits += total_hits1;
        // partial_hits += partial_hits1;
    }
    int nonrepetitive_hits = hits[0].size() + hits[1].size();
    float nonrepetitive_fraction = total_hits > 0 ? ((float) nonrepetitive_hits) / ((float) total_hits) : 1.0;
    // statistics.n_hits += nonrepetitive_hits;
    // statistics.n_partial_hits += partial_hits;

    std::vector<Nam> chains;
    std::vector<Anchor> anchors_vector;
    anchors_vector.reserve(100);
    // Rescue if requested and needed
    if (map_param.rescue_level > 1 && (nonrepetitive_hits == 0 || nonrepetitive_fraction < 0.7)) {
        // Timer rescue_timer;
        chains.clear();
        // int n_rescue_hits{0};
        // int n_partial_hits{0};
        for (int is_revcomp : {0, 1}) {
            anchors_vector.clear();
            auto [n_rescue_hits_oriented, n_partial_hits_oriented] = find_anchors_rescue(
                query_randstrobes[is_revcomp], index, map_param.rescue_cutoff, map_param.use_mcs, anchors_vector
            );
            std::sort(anchors_vector.begin(), anchors_vector.end());
            anchors_vector.erase(std::unique(anchors_vector.begin(), anchors_vector.end()), anchors_vector.end());
            collinear_chaining(anchors_vector, index.k(), is_revcomp, chains, map_param.ch_params);
        

            // n_rescue_hits += n_rescue_hits_oriented;
            // n_partial_hits += n_partial_hits_oriented;
        }
        // statistics.n_rescue_hits += n_rescue_hits;
        // statistics.n_partial_hits += partial_hits;
        // details.rescue_nams = nams.size();
        // details.nam_rescue = true;
        // statistics.tot_time_rescue += rescue_timer.duration();
    } else {
        for (int is_revcomp : {0, 1}) {
            anchors_vector.clear();
            hits_to_anchors_vector(hits[is_revcomp], index, anchors_vector);
            std::sort(anchors_vector.begin(), anchors_vector.end());
            anchors_vector.erase(std::unique(anchors_vector.begin(), anchors_vector.end()), anchors_vector.end());
            collinear_chaining(anchors_vector, index.k(), is_revcomp, chains, map_param.ch_params);
            // details.nams = nams.size();
            // statistics.tot_find_nams += nam_timer.duration();
        }
    }

    // Sort by score
    // Timer nam_sort_timer;
    std::sort(chains.begin(), chains.end(), by_score<Nam>);
    shuffle_top_chains(chains, random_engine);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

#ifdef TRACE
    std::cerr << "Query: " << record.name << '\n';
    std::cerr << "Found " << chains.size() << " NAMs\n";
    for (const auto& nam : chains) {
        std::cerr << "- " << nam << '\n';
    }
#endif

    return chains;
}
