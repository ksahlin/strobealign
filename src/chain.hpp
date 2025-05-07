#include <vector>

#include "index.hpp"
#include "indexparameters.hpp"
#include "nam.hpp"
#include "readlen.hpp"

struct Anchor {
    int query_start;
    int ref_start;

    bool operator<(const Anchor& other) const {
        return (ref_start < other.ref_start) ||
               (ref_start == other.ref_start && query_start < other.query_start);
    }

    bool operator==(const Anchor& other) const {
        return ref_start == other.ref_start && query_start == other.query_start;
    }
};

// struct Chain {
//     int query_start;
//     int query_end;
//     int ref_start;
//     int ref_end;
//     std::vector<Anchor> anchors;
//     int ref_id;
//     float score;
//     bool is_revcomp;
//
//     int ref_span() const;
//     int query_span() const;
//     int projected_ref_start() const;
// };

std::vector<Nam> get_chains(
    const klibpp::KSeq& record,
    const StrobemerIndex& index,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    std::minstd_rand& random_engine
);
