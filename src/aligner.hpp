#ifndef STROBEALIGN_ALIGNER_HPP
#define STROBEALIGN_ALIGNER_HPP

#include <string>
#include <tuple>
#include <optional>
#include "alignment.hpp"
#include "ssw/ssw_cpp.h"
#include "chain.hpp"
#include "piecewisealigner.hpp"

struct Aligner {
public:
    Aligner(AlignmentParameters parameters, int k)
        : parameters(parameters)
        , ssw_aligner(StripedSmithWaterman::Aligner(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_extend))
        , piecewise_aligner(Piecewise::Aligner(parameters, k))
    { }

    std::optional<AlignmentInfo> align(const std::string &query, const std::string &ref) const;
    AlignmentInfo align_piecewise(
        const std::string &query,
        const std::string &ref,
        std::vector<Anchor>& anchors,
        int padding
    ) const;

    AlignmentParameters parameters;

    unsigned calls_count() {
        return m_align_calls;
    }

private:
    const StripedSmithWaterman::Aligner ssw_aligner;
    const StripedSmithWaterman::Filter filter;
    const Piecewise::Aligner piecewise_aligner;
    mutable unsigned m_align_calls{0};  // no. of calls to the align() method
};

#endif
