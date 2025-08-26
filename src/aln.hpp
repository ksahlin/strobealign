#ifndef STROBEALIGN_ALN_HPP
#define STROBEALIGN_ALN_HPP

#include <string>
#include <vector>
#include <random>
#include "kseq++/kseq++.hpp"
#include "index.hpp"
#include "refs.hpp"
#include "sam.hpp"
#include "aligner.hpp"
#include "insertsizedistribution.hpp"
#include "statistics.hpp"
#include "nam.hpp"
#include "mcsstrategy.hpp"
#include "mappingparameters.hpp"

void align_or_map_paired(
    const klibpp::KSeq& record1,
    const klibpp::KSeq& record2,
    Sam& sam,
    std::string& outstring,
    AlignmentStatistics& statistics,
    InsertSizeDistribution& isize_est,
    const Aligner& aligner,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    std::minstd_rand& random_engine,
    std::vector<double> &abundances
);

void align_or_map_single(
    const klibpp::KSeq& record,
    Sam& sam,
    std::string& outstring,
    AlignmentStatistics& statistics,
    const Aligner& aligner,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    std::minstd_rand& random_engine,
    std::vector<double> &abundances
);

// Private declarations, only needed for tests

template <typename T>
bool by_score(const T& a, const T& b);

bool has_shared_substring(const std::string& read_seq, const std::string& ref_seq, int k);

void shuffle_top_nams(std::vector<Nam>& nams, std::minstd_rand& random_engine);

#endif
