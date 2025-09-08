#include "mcsstrategy.hpp"

std::ostream& operator<<(std::ostream& os, const McsStrategy& mcs_strategy) {
    if (mcs_strategy == McsStrategy::Always) { os <<  "always"; }
    else if (mcs_strategy == McsStrategy::Rescue) { os << "rescue"; }

    return os;
}
