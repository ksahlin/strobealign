#ifndef STROBEALIGN_MCSSTRATEGY_HPP
#define STROBEALIGN_MCSSTRATEGY_HPP

#include <sstream>

// Different ways in which multi-context seeds (MCS) are used
enum class McsStrategy {
    // For each strobemer, do a full lookup. If that did not generate a hit,
    // try a partial lookup.
    Always,

    // For each strobemer, do a full lookup. If after processing the entire
    // query, no hits were generated, do partial lookups of each strobemer.
    Rescue,

    // Do full lookups only
    Off,
};

std::ostream& operator<<(std::ostream& os, const McsStrategy& mcs_strategy);

#endif
