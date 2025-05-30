#include <vector>

#include "index.hpp"
#include "indexparameters.hpp"
#include "nam.hpp"
#include "readlen.hpp"

struct Anchor {
    int query_start;
    int ref_start;
    int ref_id;

    bool operator<(const Anchor& other) const {
        return (ref_id < other.ref_id) || (ref_id == other.ref_id && ref_start < other.ref_start) || (ref_id == other.ref_id && ref_start == other.ref_start && query_start < other.query_start);
    }

    bool operator==(const Anchor& other) const {
        return (ref_id == other.ref_id) && ref_start == other.ref_start && query_start == other.query_start;
    }
};

std::vector<Nam> get_chains(
    const klibpp::KSeq& record,
    const StrobemerIndex& index,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    std::minstd_rand& random_engine
);
