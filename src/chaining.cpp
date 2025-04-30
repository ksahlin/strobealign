#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#include "fastq.hpp"
#include "index.hpp"
#include "indexparameters.hpp"
#include "logger.hpp"
#include "nam.hpp"
#include "randstrobes.hpp"
#include "readlen.hpp"
#include "robin_hood.h"

struct Anchor {
    int query_start;
    int ref_start;
};

inline void add_to_anchors_map_full(
    robin_hood::unordered_map<unsigned int, std::vector<Anchor>>& anchors_map,
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
            anchors_map[index.reference_index(position)].push_back(Anchor{query_start, ref_start});
            anchors_map[index.reference_index(position)].push_back(
                Anchor{query_end - index.k(), ref_end - index.k()}
            );
            min_diff = diff;
        }
    }
}

inline void add_to_anchors_map_partial(
    robin_hood::unordered_map<unsigned int, std::vector<Anchor>>& anchors_map,
    int query_start,
    const StrobemerIndex& index,
    size_t position
) {
    for (const auto hash = index.get_main_hash(position); index.get_main_hash(position) == hash; ++position) {
        auto [ref_start, ref_end] = index.strobe_extent_partial(position);
        anchors_map[index.reference_index(position)].push_back(Anchor{query_start, ref_start});
    }
}

robin_hood::unordered_map<unsigned int, std::vector<Anchor>>
hits_to_anchors(const std::vector<Hit>& hits, const StrobemerIndex& index) {
    robin_hood::unordered_map<unsigned int, std::vector<Anchor>> anchor_map;
    anchor_map.reserve(100);  // idk why ?

    for (const Hit& hit : hits) {
        if (hit.is_partial) {
            add_to_anchors_map_partial(anchor_map, hit.query_start, index, hit.position);
        } else {
            add_to_anchors_map_full(anchor_map, hit.query_start, hit.query_end, index, hit.position);
        }
    }

    return anchor_map;
}

std::tuple<int, int, robin_hood::unordered_map<unsigned int, std::vector<Anchor>>> find_anchors_rescue(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    bool use_mcs
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
    robin_hood::unordered_map<unsigned int, std::vector<Anchor>> anchors_map;
    anchors_map.reserve(100);
    int n_hits = 0;
    int partial_hits = 0;
    std::vector<RescueHit> rescue_hits;
    rescue_hits.reserve(5000);
    for (auto& qr : query_randstrobes) {
        size_t position = index.find_full(qr.hash);
        if (position != index.end()) {
            unsigned int count = index.get_count_full(position);

            size_t position_revcomp = index.find_full(qr.hash_revcomp);
            if (position_revcomp != index.end()) {
                count += index.get_count_full(position_revcomp);
            }
            RescueHit rh{position, count, qr.start, qr.end, false};
            rescue_hits.push_back(rh);
        } else if (use_mcs) {
            size_t partial_pos = index.find_partial(qr.hash);
            if (partial_pos != index.end()) {
                unsigned int partial_count = index.get_count_partial(partial_pos);
                size_t position_revcomp = index.find_partial(qr.hash_revcomp);
                if (position_revcomp != index.end()) {
                    partial_count += index.get_count_partial(position_revcomp);
                }
                RescueHit rh{partial_pos, partial_count, qr.start, qr.start + index.k(), true};
                rescue_hits.push_back(rh);
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
            add_to_anchors_map_partial(anchors_map, rh.query_start, index, rh.position);
        } else {
            add_to_anchors_map_full(anchors_map, rh.query_start, rh.query_end, index, rh.position);
        }
        cnt++;
        n_hits++;
    }

    return {n_hits, partial_hits, anchors_map};
}

struct Chain {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;

    // for debugging :p
    std::vector<Anchor> anchors;

    int ref_id;
    float score;
    bool is_revcomp = false;

    int ref_span() const { return ref_end - ref_start; }

    int query_span() const { return query_end - query_start; }

    int projected_ref_start() const { return std::max(0, ref_start - query_start); }
};

static inline int compute_score(const Anchor& ai, const Anchor& aj, const int k) {
    int dq = ai.query_start - aj.query_start;
    int dr = ai.ref_start - aj.ref_start;
    if (dq < 0 || dr < 0)
        return INT32_MIN;

    int dd = std::abs(dr - dq);
    int dg = std::min(dq, dr);
    int score = std::min(k, dg);

    const float gap_open = 1.0f;
    const float gap_extend = 0.05f;
    float penalty = gap_open * dd + gap_extend * dg;
    score -= static_cast<int>(penalty);

    return score;
}

void collinear_chaining(
    int ref_id,
    std::vector<Anchor>& anchors,
    const int k,
    bool sort,
    bool is_revcomp,
    std::vector<Chain>& chains  // inout
) {
    int n = anchors.size();
    if (n == 0) {
        return;
    }

    if (sort) {
        std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
            return (a.ref_start < b.ref_start) ||
                   (a.ref_start == b.ref_start && a.query_start < b.query_start);
        });
    }

    std::vector<int> dp(n, 0);
    std::vector<int> backtrack(n, -1);

    int best_score = 0;
    int best_index = -1;

    for (int i = 0; i < n; ++i) {
        dp[i] = k;

        for (int j = 0; j < i; ++j) {
            int score = compute_score(anchors[i], anchors[j], k);
            if (score == INT32_MIN)
                continue;

            int new_score = dp[j] + score;
            if (new_score > dp[i]) {
                dp[i] = new_score;
                backtrack[i] = j;
            }
        }
        if (dp[i] > best_score) {
            best_score = dp[i];
            best_index = i;
        }
    }

    if (best_index < 0)
        return;

    std::vector<Anchor> chain_anchors;
    for (int i = best_index; i >= 0; i = backtrack[i]) {
        chain_anchors.push_back(anchors[i]);
    }
    std::reverse(chain_anchors.begin(), chain_anchors.end());

    const Anchor& first = chain_anchors.front();
    const Anchor& last = chain_anchors.back();

    chains.push_back(Chain{
        first.query_start, last.query_start + k, first.ref_start, last.ref_start + k,
        std::move(chain_anchors), ref_id, float(best_score), is_revcomp

    });
}

template <typename T>
bool by_score(const T& a, const T& b) {
    return a.score > b.score;
}

void shuffle_top_chains(std::vector<Chain>& chains, std::minstd_rand& random_engine) {
    if (chains.empty()) {
        return;
    }
    auto best_score = chains[0].score;
    auto it = std::find_if(chains.begin(), chains.end(), [&](const Chain& chain) {
        return chain.score != best_score;
    });
    if (it > chains.begin() + 1) {
        std::shuffle(chains.begin(), it, random_engine);
    }
}

std::vector<Chain> get_chains(
    const klibpp::KSeq& record,
    const StrobemerIndex& index,
    // AlignmentStatistics& statistics,
    // Details& details,
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
    bool sorting_needed = false;
    std::array<std::vector<Hit>, 2> hits;
    for (int is_revcomp : {0, 1}) {
        int total_hits1, partial_hits1;
        bool sorting_needed1;
        std::tie(total_hits1, partial_hits1, sorting_needed1, hits[is_revcomp]) =
            find_hits(query_randstrobes[is_revcomp], index, map_param.use_mcs);
        sorting_needed = sorting_needed || sorting_needed1;
        total_hits += total_hits1;
        // partial_hits += partial_hits1;
    }
    int nonrepetitive_hits = hits[0].size() + hits[1].size();
    float nonrepetitive_fraction = total_hits > 0 ? ((float) nonrepetitive_hits) / ((float) total_hits) : 1.0;
    // statistics.n_hits += nonrepetitive_hits;
    // statistics.n_partial_hits += partial_hits;

    std::vector<Chain> chains;

    // Rescue if requested and needed
    if (map_param.rescue_level > 1 && (nonrepetitive_hits == 0 || nonrepetitive_fraction < 0.7)) {
        // Timer rescue_timer;
        chains.clear();
        // int n_rescue_hits{0};
        // int n_partial_hits{0};
        for (int is_revcomp : {0, 1}) {
            auto [n_rescue_hits_oriented, n_partial_hits_oriented, anchors_map] = find_anchors_rescue(
                query_randstrobes[is_revcomp], index, map_param.rescue_cutoff, map_param.use_mcs
            );

            for (auto& [ref_id, anchors] : anchors_map) {
                collinear_chaining(ref_id, anchors, index.k(), sorting_needed, is_revcomp, chains);
            }

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
            auto anchors_map = hits_to_anchors(hits[is_revcomp], index);

            for (auto& [ref_id, anchors] : anchors_map) {
                collinear_chaining(ref_id, anchors, index.k(), sorting_needed, is_revcomp, chains);
            }
            // details.nams = nams.size();
            // statistics.tot_find_nams += nam_timer.duration();
        }
    }

    // Sort by score
    // Timer nam_sort_timer;
    std::sort(chains.begin(), chains.end(), by_score<Chain>);
    shuffle_top_chains(chains, random_engine);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

#ifdef TRACE
    std::cerr << "Query: " << record.name << '\n';
    std::cerr << "Found " << nams.size() << " NAMs\n";
    for (const auto& nam : nams) {
        std::cerr << "- " << nam << '\n';
    }
#endif

    return chains;
}

int main() {
    Logger& logger = Logger::get();
    logger.set_level(LOG_INFO);

    const std::string fasta_file = "/home/nico/data/strobealign/drosophila.fa";
    const std::string fastq_file = "/home/nico/data/strobealign/sim3-drosophila-500.1.fastq";

    References refs = References::from_fasta(fasta_file);
    logger.info() << "Loaded " << refs.size() << " reference(s).\n";

    InputBuffer input_buffer(fastq_file, "", 100000, false);
    uint64_t read_len = estimate_read_length(input_buffer);
    logger.info() << "Estimated read length: " << read_len << " bp\n";

    IndexParameters index_parameters = IndexParameters::from_read_length(read_len);
    StrobemerIndex index(refs, index_parameters);
    index.populate(0.0002f, 1);
    logger.info() << "Index built: " << index.stats.tot_strobemer_count << " strobemers.\n";

    std::unique_ptr<RewindableFile> fq = open_fastq(const_cast<std::string&>(fastq_file));
    std::vector<klibpp::KSeq> records = fq->stream().read(100000);

    std::minstd_rand random_engine;
    MappingParameters map_param;
    std::ofstream out("alignment-col.txt");

    for (const klibpp::KSeq& record : records) {
        std::vector<Chain> chains = get_chains(record, index, map_param, index_parameters, random_engine);

        logger.info() << "Read: " << record.name << "\tnb chains: " << chains.size();
        out << "Read: " << record.name << "\tnb chains: " << chains.size();

        if (chains.size() != 0) {
            const Chain best_chain = chains[0];

            logger.info() << "  Best Chain at " << refs.names[best_chain.ref_id]
                          << "\tori: " << (best_chain.is_revcomp ? "rev" : "fwd")
                          << "\tscore: " << best_chain.score << "\tanchors: " << best_chain.anchors.size()
                          << "\tq: " << best_chain.query_start << "-" << best_chain.query_end
                          << "\tr: " << best_chain.ref_start << "-" << best_chain.ref_end;

            out << "  Best Chain at " << refs.names[best_chain.ref_id]
                << "\tori: " << (best_chain.is_revcomp ? "rev" : "fwd") << "\tscore: " << best_chain.score
                << "\tanchors: " << best_chain.anchors.size() << "\tq: " << best_chain.query_start << "-"
                << best_chain.query_end << "\tr: " << best_chain.ref_start << "-" << best_chain.ref_end;

            // for (const Anchor& anchor : best_chain.anchors) {
            //     logger.info() << "    Anchor: q=" << anchor.query_start << " r=" << anchor.ref_start << "\n";
            // }
        }

        logger.info() << "\n";
        out << "\n";
    }

    out.close();

    return 0;
}
