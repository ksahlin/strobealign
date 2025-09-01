#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/array.h>
#include <nanobind/make_iterator.h>
#include <iostream>
#include <sstream>
#include "randstrobes.hpp"
#include "revcomp.hpp"
#include "refs.hpp"
#include "index.hpp"
#include "nam.hpp"
#include "logger.hpp"


namespace nb = nanobind;
using namespace nb::literals;

class SyncmerIteratorWrapper {
public:
    explicit SyncmerIteratorWrapper(const std::string& seq, SyncmerParameters parameters)
        : m_seq(seq), it(m_seq, parameters) { }

    Syncmer next() {
        return it.next();
    }

private:
    std::string m_seq;
    SyncmerIterator it;
};

class Record {
public:
    Record(const std::string& name, const std::string& sequence)
        : m_name(name)
        , m_sequence(sequence)
    { }

    const std::string& name() const {
        return m_name;
    }

    const std::string& sequence() const {
        return m_sequence;
    }
private:
    const std::string& m_name;
    const std::string& m_sequence;
};

NB_MODULE(strobealign_extension, m_) {
    (void) m_;
    // Add definitions to the *parent* module
    nb::module_ m = nb::module_::import_("strobealign");

    m.doc() = "strobealign aligns short reads using dynamic seed size with strobemers";
    nb::enum_<LOG_LEVELS>(m, "LOG_LEVELS", nb::is_arithmetic())
        .value("LOG_DEBUG", LOG_LEVELS::LOG_DEBUG)
        .value("LOG_INFO", LOG_LEVELS::LOG_INFO)
        .value("LOG_WARNING", LOG_LEVELS::LOG_WARNING)
        .value("LOG_ERROR", LOG_LEVELS::LOG_ERROR)
        .export_values();

    nb::class_<Logger>(m, "Logger")
        .def_static("get", &Logger::get, nb::rv_policy::reference)
        .def("set_level", &Logger::set_level)
    ;
    nb::class_<SyncmerParameters>(m, "SyncmerParameters")
        .def_ro("k", &SyncmerParameters::k)
        .def_ro("s", &SyncmerParameters::s)
        .def_ro("t", &SyncmerParameters::t_syncmer)
    ;
    nb::class_<RandstrobeParameters>(m, "RandstrobeParameters")
        .def_ro("max_dist", &RandstrobeParameters::max_dist)
        .def_ro("w_min", &RandstrobeParameters::w_min)
        .def_ro("w_max", &RandstrobeParameters::w_max)
        .def_ro("q", &RandstrobeParameters::q)
    ;
    nb::class_<SyncmerIteratorWrapper>(m, "SyncmerIterator")
        .def(nb::init<const std::string&, SyncmerParameters>())
        .def("__next__", [](SyncmerIteratorWrapper& siw) -> std::pair<size_t, randstrobe_hash_t> {
            auto syncmer = siw.next();
            if (syncmer.is_end()) {
                throw nb::stop_iteration();
            } else {
                return std::make_pair(syncmer.position, syncmer.hash);
            }
        })
        .def("__iter__", [](SyncmerIteratorWrapper& siw) -> SyncmerIteratorWrapper { return siw; })
    ;
    nb::class_<References>(m, "References")
        .def(nb::init())
        .def("add", &References::add)
        .def_static("from_fasta", &References::from_fasta)
        .def("__getitem__", [](const References& refs, size_t i) {
            return Record(refs.names[i], refs.sequences[i]);
        })
        .def("__len__", [](const References& refs) { return refs.sequences.size(); })
    ;
    nb::class_<Record>(m, "Record", "FASTA record")
        .def_prop_ro("name", &Record::name)
        .def_prop_ro("sequence", &Record::sequence)
    ;
    nb::class_<IndexParameters>(m, "IndexParameters")
        .def_static("from_read_length", [](int r, int k = IndexParameters::DEFAULT, int s = IndexParameters::DEFAULT, int l = IndexParameters::DEFAULT, int u = IndexParameters::DEFAULT) {
            return IndexParameters::from_read_length(r, k, s, l, u);
        }, "r"_a, "k"_a = IndexParameters::DEFAULT, "s"_a = IndexParameters::DEFAULT, "l"_a = IndexParameters::DEFAULT, "u"_a = IndexParameters::DEFAULT)
        .def_ro("syncmer", &IndexParameters::syncmer)
        .def_ro("randstrobe", &IndexParameters::randstrobe)
    ;
    nb::class_<RefRandstrobe>(m, "RefRandstrobe", "Randstrobe on a reference")
        .def_prop_ro("position", &RefRandstrobe::position)
        .def_prop_ro("hash", &RefRandstrobe::hash)
        .def_prop_ro("reference_index", &RefRandstrobe::reference_index)
        .def_prop_ro("strobe2_offset", &RefRandstrobe::strobe2_offset)
    ;
    nb::bind_vector<std::vector<RefRandstrobe>>(m, "RefRandstrobeVector");
    nb::class_<StrobemerIndex>(m, "StrobemerIndex")
        .def(nb::init<References&, IndexParameters&>())
        .def("find", [](const StrobemerIndex& index, uint64_t key) -> std::vector<RefRandstrobe> {
            std::vector<RefRandstrobe> v;
            auto position = index.find_full(key);
            while (position != index.end() && index.get_hash(position) == key) {
                v.push_back(index.get_randstrobe(position));
                position++;
            }
            return v;
        })
        .def("populate", &StrobemerIndex::populate, "f"_a = 0.0002, "threads"_a = 1)
        .def("reference_index", &StrobemerIndex::reference_index)
        .def("get_strobe1_position", &StrobemerIndex::get_strobe1_position)
        .def("strobe2_offset", &StrobemerIndex::strobe2_offset)
        .def_ro("filter_cutoff", &StrobemerIndex::filter_cutoff)
        .def_prop_ro("k", &StrobemerIndex::k)
    ;

    nb::class_<QueryRandstrobe>(m, "QueryRandstrobe")
        .def_ro("hash", &QueryRandstrobe::hash)
        .def_ro("start", &QueryRandstrobe::start)
        .def_ro("end", &QueryRandstrobe::end)
        .def("__repr__", [](const QueryRandstrobe& qr) {
            std::stringstream s; s << qr; return s.str();
        })
    ;
    nb::bind_vector<std::vector<QueryRandstrobe>>(m, "QueryRandstrobeVector");

    nb::class_<Nam>(m, "Nam")
        .def_ro("query_start", &Nam::query_start)
        .def_ro("query_end", &Nam::query_end)
        .def_ro("ref_start", &Nam::ref_start)
        .def_ro("ref_end", &Nam::ref_end)
        .def_ro("score", &Nam::score)
        .def_ro("n_hits", &Nam::n_matches)
        .def_ro("reference_index", &Nam::ref_id)
        .def_rw("is_revcomp", &Nam::is_revcomp)
        .def_prop_ro("ref_span", &Nam::ref_span)
        .def_prop_ro("query_span", &Nam::query_span)
        .def("__repr__", [](const Nam& nam) {
            std::stringstream s; s << nam; return s.str();
        })
    ;
    nb::bind_vector<std::vector<Nam>>(m, "NamVector");
    nb::bind_vector<std::vector<Hit>>(m, "HitVector");
    nb::bind_vector<std::vector<Match>>(m, "MatchVector");

    nb::class_<Hit>(m, "Hit")
        .def_ro("position", &Hit::position)
        .def_ro("query_start", &Hit::query_start)
        .def_ro("query_end", &Hit::query_end)
        .def_ro("is_partial", &Hit::is_partial)
        .def("__repr__", [](const Hit& hit) {
            std::stringstream s; s << hit; return s.str();
        })
    ;

    nb::class_<Match>(m, "Match")
        .def_ro("query_start", &Match::query_start)
        .def_ro("query_end", &Match::query_end)
        .def_ro("ref_start", &Match::ref_start)
        .def_ro("ref_end", &Match::ref_end)
        .def("__repr__", [](const Match& match) {
            std::stringstream s; s << match; return s.str();
        })
    ;

    m.def("reverse_complement", &reverse_complement);
    m.def("randstrobes_query", &randstrobes_query);

    nb::enum_<McsStrategy>(m, "McsStrategy")
        .value("Always", McsStrategy::Always)
        .value("Rescue", McsStrategy::Rescue)
        .value("Off", McsStrategy::Off)
        .value("FirstStrobe", McsStrategy::FirstStrobe)
        .export_values();

    m.def("find_hits", [](const std::vector<QueryRandstrobe>& query_randstrobes, const StrobemerIndex& index, McsStrategy mcs_strategy) -> std::vector<Hit> {
        auto [total_hits, partial_hits, sorting_needed, hits] = find_hits(query_randstrobes, index, mcs_strategy);
        return hits;
    }, nb::arg("query_randstrobes"), nb::arg("index"), nb::arg("mcs_strategy"));

    m.def("hits_to_matches", [](const std::vector<Hit>& hits, const StrobemerIndex& index) -> std::unordered_map<unsigned int, std::vector<Match>> {
        auto rhmap = hits_to_matches(hits, index);
        std::unordered_map<unsigned int, std::vector<Match>> map;
        for (const auto& [key, value] : rhmap)
            map.insert({key, value});

        return map;
    }, nb::arg("hits"), nb::arg("index"));

    m.def("merge_matches_into_nams", [](
            const std::unordered_map<unsigned int, std::vector<Match>>& matches_map,
            int k,
            bool sort,
            bool is_revcomp
        ) -> std::vector<Nam> {
        std::vector<Nam> nams;

        robin_hood::unordered_map<unsigned int, std::vector<Match>> rhmap;
        for (const auto& [key, value] : matches_map)
            rhmap.insert({key, value});
        merge_matches_into_nams(rhmap, k, sort, is_revcomp, nams);

        return nams;
    }, nb::arg("matches_map"), nb::arg("k"), nb::arg("sort"), nb::arg("is_revcomp"));
}
