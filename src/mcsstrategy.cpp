#include "mcsstrategy.hpp"

std::ostream& operator<<(std::ostream& os, const McsStrategy& mcs_strategy) {
    if (mcs_strategy == McsStrategy::Always) { os <<  "always"; }
    else if (mcs_strategy == McsStrategy::Rescue) { os << "rescue"; }
    else if (mcs_strategy == McsStrategy::FirstStrobe) { os << "first strobe"; }
    else if (mcs_strategy == McsStrategy::Off) { os << "off"; }

    return os;
}

std::ostream& operator<<(std::ostream& os, const SearchStrategy& search_strategy) {
    if (search_strategy == SearchStrategy::First) { os <<  "first only"; }
    else if (search_strategy == SearchStrategy::Second) { os << "second only"; }
    else if (search_strategy == SearchStrategy::Rescue) { os << "second as rescue"; }

    return os;
}
