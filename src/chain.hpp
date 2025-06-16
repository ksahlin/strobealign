#include <vector>

#include "index.hpp"
#include "indexparameters.hpp"
#include "nam.hpp"
#include "readlen.hpp"

struct Anchor {
    uint query_start;
    uint ref_start;
    uint ref_id;

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

static inline float mg_log2(float x)
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}
