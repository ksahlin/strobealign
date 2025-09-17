#include <string>
#include <vector>
#include "aligner.hpp"
#include "piecewisealigner.hpp"
#include "chain.hpp"
#include "sam.hpp"

using namespace klibpp;

void align_single_piecewse(
    const Aligner &aligner,
    const AlignmentParameters& scoring_params,
    Sam& sam,
    const std::vector<Chain>& chains,
    const KSeq& record,
    const int k,
    const References& references,
    Details& details,
    float dropoff_threshold,
    int max_tries,
    unsigned max_secondary,
    std::minstd_rand& random_engine
); 
