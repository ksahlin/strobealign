#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/make_iterator.h>
#include <iostream>
#include <sstream>
#include "randstrobes.hpp"
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
        .def_ro("position", &RefRandstrobe::position)
        .def_ro("hash", &RefRandstrobe::hash)
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
        .def_ro("filter_cutoff", &StrobemerIndex::filter_cutoff)
        .def_prop_ro("k", &StrobemerIndex::k)
    ;

    nb::class_<QueryRandstrobe>(m, "QueryRandstrobe")
        .def_ro("hash", &QueryRandstrobe::hash)
        .def_ro("start", &QueryRandstrobe::start)
        .def_ro("end", &QueryRandstrobe::end)
        .def_ro("is_reverse", &QueryRandstrobe::is_reverse)
        .def("__repr__", [](const QueryRandstrobe& qr) {
            std::stringstream s; s << qr; return s.str();
        })
    ;
    nb::bind_vector<QueryRandstrobeVector>(m, "QueryRandstrobeVector");

    nb::class_<Nam>(m, "Nam")
        .def_ro("query_start", &Nam::query_start)
        .def_ro("query_end", &Nam::query_end)
        .def_ro("ref_start", &Nam::ref_start)
        .def_ro("ref_end", &Nam::ref_end)
        .def_ro("score", &Nam::score)
        .def_ro("n_hits", &Nam::n_matches)
        .def_ro("reference_index", &Nam::ref_id)
        .def_rw("is_rc", &Nam::is_rc)
        .def_prop_ro("ref_span", &Nam::ref_span)
        .def_prop_ro("query_span", &Nam::query_span)
        .def("__repr__", [](const Nam& nam) {
            std::stringstream s; s << nam; return s.str();
        })
    ;
    nb::bind_vector<std::vector<Nam>>(m, "NamVector");

    m.def("randstrobes_query", &randstrobes_query);
    m.def("find_nams", [](const QueryRandstrobeVector &query_randstrobes, const StrobemerIndex& index, bool use_mcs) -> std::vector<Nam> {
        auto [nonrepetitive_fraction, n_hits, nams] = find_nams(query_randstrobes, index, use_mcs);
        return nams;
    }, nb::arg("query_randstrobes"), nb::arg("index"), nb::arg("use_mcs"));
}
