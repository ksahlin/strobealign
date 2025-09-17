#ifndef PIECEWISE_ALIGNER_H
#define PIECEWISE_ALIGNER_H

#include <string_view>
#include <cstddef>
#include <vector>
#include "block-aligner/c/block_aligner.h"
#include "alignment.hpp"
#include "chain.hpp"

namespace Piecewise {

struct AlignmentResult {
    int score;
    size_t query_start;
    size_t query_end;
    size_t ref_start;
    size_t ref_end;
    Cigar cigar;
};

class Aligner {
private:
    AlignmentParameters params;
    int k;
    SizeRange range;
    Gaps gaps;
    AAMatrix* matrix;
    int32_t x_drop_threshold;

    Cigar build_cigar(const Cigar* cigar, size_t cigar_len) const;
    Cigar build_cigar_swap_ID(const Cigar* cigar, size_t cigar_len) const;
    Cigar build_cigar_reverse_swap_ID(const Cigar* cigar, size_t cigar_len) const;
    
    AlignmentResult global_alignment(const std::string_view& query, const std::string_view& ref) const;
    AlignmentResult xdrop_alignment(const std::string_view& query, const std::string_view& ref, bool reverse = false) const;
    
    void align_before_first_anchor(
        const std::string& reference,
        const std::string& query,
        const Anchor& first_anchor, 
        const int padding,
        AlignmentInfo* result
    ) const;
    
    void align_after_last_anchor(
        const std::string& reference,
        const std::string& query,
        const Anchor& last_anchor,
        const int padding,
        AlignmentInfo* result
    ) const;

public:
    explicit Aligner(const AlignmentParameters& params, int k);
    
    ~Aligner();
    
    AlignmentInfo piecewise_extension_alignment(
        const std::string& reference,
        const std::string& query,
        const std::vector<Anchor>& anchors,
        const int padding
    ) const;

    const AlignmentParameters& get_parameters() const { return params; }
};

} // namespace Piecewise

#endif
