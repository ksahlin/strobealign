#ifndef PIECEWISE_ALIGNER_H
#define PIECEWISE_ALIGNER_H

#include <string_view>
#include <cstddef>
#include "block-aligner/c/block_aligner.h"
#include "aligner.hpp"
#include "cigar.hpp"

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
    SizeRange range;
    Gaps gaps;
    AAMatrix* matrix;
    int32_t x_drop_threshold;

    Cigar build_cigar(const Cigar* cigar, size_t cigar_len) const;
    Cigar build_cigar_swap_ID(const Cigar* cigar, size_t cigar_len) const;
    Cigar build_cigar_reverse_swap_ID(const Cigar* cigar, size_t cigar_len) const;

public:
    explicit Aligner(const AlignmentParameters& params);
    
    ~Aligner();
    
    Aligner(const Aligner&) = delete;
    Aligner& operator=(const Aligner&) = delete;
    
    Aligner(Aligner&& other) noexcept;
    Aligner& operator=(Aligner&& other) noexcept;
    
    AlignmentResult global_alignment(const std::string_view& query, const std::string_view& ref) const;
    AlignmentResult xdrop_alignment(const std::string_view& query, const std::string_view& ref, bool reverse = false) const;

    const AlignmentParameters& get_parameters() const { return params; }
};

} // namespace Piecewise

#endif
