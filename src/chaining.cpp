#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "fastq.hpp"
#include "index.hpp"
#include "indexparameters.hpp"
#include "logger.hpp"
#include "nam.hpp"
#include "randstrobes.hpp"

static inline int
compute_score_match(const Match& ai, const Match& aj, float gap_diff_penalty, float gap_span_penalty) {
    int q_gain = ai.query_end - aj.query_end;
    int r_gain = ai.ref_end - aj.ref_end;
    if (q_gain <= 0 || r_gain <= 0)
        return INT32_MIN;

    int gap_diff = std::abs(r_gain - q_gain);
    int gap_span = std::min(q_gain, r_gain);

    int score = gap_span;

    if (gap_diff > 0 || gap_span > (aj.query_end - aj.query_start)) {
        float lin_pen = gap_diff_penalty * gap_diff + gap_span_penalty * gap_span;
        float log_pen = gap_diff >= 1 ? std::log2(gap_diff + 1.0f) : 0.0f;
        score -= static_cast<int>(lin_pen + 0.5f * log_pen);
    }
    return score;
}

std::vector<Match> chain_matches(std::vector<Match> matches) {
    if (matches.empty())
        return {};

    std::sort(matches.begin(), matches.end(), [](Match& a, Match& b) {
        if (a.ref_start == b.ref_start) {
            return a.query_start < b.query_start;
        }
        return a.ref_start < b.ref_start;
    });

    const float gap_diff_pen = 0.5f;
    const float gap_span_pen = 0.1f;
    int n = static_cast<int>(matches.size());

    std::vector<int> dp(n), prev(n, -1);
    int best_i = 0, best_score = INT32_MIN;

    for (int i = 0; i < n; ++i) {
        dp[i] = matches[i].query_end - matches[i].query_start;

        for (int j = 0; j < i; ++j) {
            int score = compute_score_match(matches[i], matches[j], gap_diff_pen, gap_span_pen);
            if (score == INT32_MIN)
                continue;

            score += dp[j];

            if (score > dp[i]) {
                dp[i] = score;
                prev[i] = j;
            }
        }
        if (dp[i] > best_score) {
            best_score = dp[i];
            best_i = i;
        }
    }

    std::vector<Match> chain;
    for (int cur = best_i; cur >= 0; cur = prev[cur]) {
        chain.push_back(matches[cur]);
    }
    std::reverse(chain.begin(), chain.end());
    return chain;
}

int main() {
    Logger& logger = Logger::get();
    logger.set_level(LOG_INFO);

    const std::string fasta_file = "/home/nico/data/strobealign/drosophila.fa";
    References refs = References::from_fasta(fasta_file);
    logger.info() << "Loaded " << refs.size() << " reference(s).";

    int read_len = 150;
    IndexParameters idx_params = IndexParameters::from_read_length(read_len);
    StrobemerIndex index(refs, idx_params);
    index.populate(0.0002f, 1);
    logger.info() << "Index built: " << index.stats.tot_strobemer_count << " strobemers.";

    const std::string fastq_file = "/home/nico/data/strobealign/sim3-drosophila-500.1.fastq";
    std::unique_ptr<RewindableFile> fq = open_fastq(const_cast<std::string&>(fastq_file));
    std::vector<klibpp::KSeq> records = fq->stream().read(1000);

    size_t best_chain_len = 0;
    std::string best_read;
    int best_ori = -1;
    std::vector<Match> best_chain;

    for (klibpp::KSeq& rec : records) {
        // logger.info() << "Read: " << rec.name << "\n";

        std::array<std::vector<QueryRandstrobe>, 2> qr = randstrobes_query(rec.seq, idx_params);
        for (int ori : {0, 1}) {
            std::tuple<uint64_t, uint64_t, bool, std::vector<Hit>> result = find_hits(qr[ori], index, false);
            std::vector<Hit>& hits = std::get<3>(result);
            // logger.info() << " ori=" << (ori ? "rev" : "fwd") << " hits=" << hits.size() << "\n";

            robin_hood::unordered_map<uint, std::vector<Match>> matches_map = hits_to_matches(hits, index);

            std::vector<Match> all_matches;
            for (robin_hood::pair<uint, std::vector<Match>>& kv : matches_map) {
                std::vector<Match>& mv = kv.second;
                all_matches.insert(all_matches.end(), mv.begin(), mv.end());
            }

            std::vector<Match> chain = chain_matches(std::move(all_matches));
            // logger.info() << "  chain_size=" << chain.size() << "\n";

            if (chain.size() > best_chain_len) {
                best_chain_len = chain.size();
                best_chain = chain;
                best_read = rec.name;
                best_ori = ori;
            }
        }
    }

    logger.info() << "\n=== Best match chain ===\n"
                  << " Read: " << best_read << "\n"
                  << " Ori : " << (best_ori ? "rev" : "fwd") << "\n"
                  << " Len : " << best_chain_len << "\n";
    for (const Match& m : best_chain) {
        logger.info() << "    q=[" << m.query_start << "," << m.query_end << "]  r=[" << m.ref_start << ","
                      << m.ref_end << "]" << "\n";
    }
    return 0;
}
