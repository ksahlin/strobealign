#ifndef STROBEALIGN_STATISTICS_HPP
#define STROBEALIGN_STATISTICS_HPP

#include <chrono>
#include <cstdint>

#include "hits.hpp"

/* Details about aligning a single or paired-end read */
struct Details {
    bool nam_rescue{false}; // find_nams_rescue() was needed
    uint32_t nams{0};  // No. of NAMs found
    uint32_t rescue_nams{0}; // No. of NAMs found during rescue
    uint32_t inconsistent_nams{0};
    uint32_t mate_rescue{0}; // No. of times rescue by local alignment was attempted
    uint32_t tried_alignment{0}; // No. of computed alignments (get_alignment or rescue_mate)
    uint32_t gapped{0};  // No. of gapped alignments computed (in get_alignment)
    uint32_t best_alignments{0}; // No. of best alignments with same score

    HitsDetails hits;

    Details& operator+=(const Details& other) {
        nam_rescue = nam_rescue || other.nam_rescue;
        nams += other.nams;
        rescue_nams += other.rescue_nams;
        inconsistent_nams += other.inconsistent_nams;
        mate_rescue += other.mate_rescue;
        tried_alignment += other.tried_alignment;
        gapped += other.gapped;
        best_alignments += other.best_alignments;
        hits += other.hits;
        return *this;
    }
};

struct AlignmentStatistics {
    std::chrono::duration<double> tot_read_file{0};
    std::chrono::duration<double> tot_construct_strobemers{0};
    std::chrono::duration<double> time_hit_finding{0};
    std::chrono::duration<double> time_chaining{0};
    std::chrono::duration<double> tot_time_rescue{0};
    std::chrono::duration<double> tot_sort_nams{0};
    std::chrono::duration<double> tot_extend{0};

    uint64_t n_reads{0};
    uint64_t n_randstrobes{0};
    uint64_t n_rescue_hits{0};
    uint64_t n_rescue_partial_hits{0};
    uint64_t n_anchors{0};
    uint64_t n_nams{0};
    uint64_t n_rescue_nams{0};
    uint64_t tot_aligner_calls{0};
    uint64_t tot_rescued{0};
    uint64_t tried_alignment{0};
    uint64_t inconsistent_nams{0};
    uint64_t nam_rescue{0};

    HitsDetails hits;

    AlignmentStatistics operator+=(const AlignmentStatistics& other) {
        this->tot_read_file += other.tot_read_file;
        this->tot_construct_strobemers += other.tot_construct_strobemers;
        this->time_hit_finding += other.time_hit_finding;
        this->time_chaining += other.time_chaining;
        this->tot_time_rescue += other.tot_time_rescue;
        this->tot_sort_nams += other.tot_sort_nams;
        this->tot_extend += other.tot_extend;
        this->n_reads += other.n_reads;
        this->n_randstrobes += other.n_randstrobes;
        this->n_rescue_hits += other.n_rescue_hits;
        this->n_rescue_partial_hits += other.n_rescue_partial_hits;
        this->n_anchors += other.n_anchors;
        this->n_nams += other.n_nams;
        this->n_rescue_nams += other.n_rescue_nams;
        this->tot_aligner_calls += other.tot_aligner_calls;
        this->tot_rescued += other.tot_rescued;
        this->tried_alignment += other.tried_alignment;
        this->inconsistent_nams += other.inconsistent_nams;
        this->nam_rescue += other.nam_rescue;
        this->hits += other.hits;
        return *this;
    }

    AlignmentStatistics operator+=(const Details& details) {
        this->n_nams += details.nams;
        this->n_rescue_nams += details.rescue_nams;
        this->nam_rescue += details.nam_rescue;
        this->tot_rescued += details.mate_rescue;
        this->tried_alignment += details.tried_alignment;
        this->inconsistent_nams += details.inconsistent_nams;
        this->hits += details.hits;

        return *this;
    }
};

#endif
