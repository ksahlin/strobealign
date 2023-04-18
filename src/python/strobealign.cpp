#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
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
    explicit SyncmerIteratorWrapper(const std::string& seq, size_t k, size_t s, size_t t)
        : m_seq(seq), it(m_seq, k, s, t) { }

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

    nb::enum_<LOG_LEVELS>(m, "LOG_LEVELS")
        .value("LOG_DEBUG", LOG_LEVELS::LOG_DEBUG)
        .value("LOG_INFO", LOG_LEVELS::LOG_INFO)
        .value("LOG_WARNING", LOG_LEVELS::LOG_WARNING)
        .value("LOG_ERROR", LOG_LEVELS::LOG_ERROR)
        .export_values();

    nb::class_<Logger>(m, "Logger")
        .def_static("get", &Logger::get, nb::rv_policy::reference)
        .def("set_level", &Logger::set_level)
    ;
    nb::class_<SyncmerIteratorWrapper>(m, "SyncmerIterator")
        .def(nb::init<const std::string&, size_t, size_t, size_t>())
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
    nb::class_<Record>(m, "Record")
        .def_prop_ro("name", &Record::name)
        .def_prop_ro("sequence", &Record::sequence)
    ;
    nb::class_<IndexParameters>(m, "IndexParameters")
        .def_static("from_read_length", [](int r) { return IndexParameters::from_read_length(r); })
        .def_ro("k", &IndexParameters::k)
        .def_ro("s", &IndexParameters::s)
        .def_ro("t", &IndexParameters::t_syncmer)
        .def_ro("max_dist", &IndexParameters::max_dist)
        .def_ro("w_min", &IndexParameters::w_min)
        .def_ro("w_max", &IndexParameters::w_max)
        .def_ro("q", &IndexParameters::q)
    ;
    nb::class_<RefRandstrobe>(m, "RefRandstrobe")
        .def_ro("position", &RefRandstrobe::position)
        .def_prop_ro("reference_index", &RefRandstrobe::reference_index)
        .def_prop_ro("strobe2_offset", &RefRandstrobe::strobe2_offset)
    ;
    nb::bind_vector<RefRandstrobeVector>(m, "RefRandstrobeVector");
    nb::class_<StrobemerIndex>(m, "StrobemerIndex")
        .def(nb::init<References&, IndexParameters&>())
        .def("find", [](const StrobemerIndex& index, uint64_t key) -> RefRandstrobeVector {
            RefRandstrobeVector v;
            auto it = index.find(key);
            if (it != index.end()) {
                it->second.for_each(index.flat_vector, [&v](const RefRandstrobe& rr) {
                    v.push_back(rr);
                });
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
        .def_ro("query_start", &Nam::query_s)
        .def_ro("query_end", &Nam::query_e)
        .def_ro("ref_start", &Nam::ref_s)
        .def_ro("ref_end", &Nam::ref_e)
        .def_ro("score", &Nam::score)
        .def_ro("n_hits", &Nam::n_hits)
        .def_ro("reference_index", &Nam::ref_id)
        .def_ro("is_rc", &Nam::is_rc)
        .def_prop_ro("ref_span", &Nam::ref_span)
        .def_prop_ro("query_span", &Nam::query_span)
        .def("__repr__", [](const Nam& nam) {
            std::stringstream s; s << nam; return s.str();
        })
    ;
    nb::bind_vector<std::vector<Nam>>(m, "NamVector");

    m.def("randstrobes_query", &randstrobes_query);
    m.def("find_nams", [](const QueryRandstrobeVector &query_randstrobes, const StrobemerIndex& index) -> std::vector<Nam> {
        return find_nams(query_randstrobes, index).second;
    });
}
