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


NB_MODULE(pystrobealign, m) {
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
        .def("__iter__", [](SyncmerIteratorWrapper& siw) -> SyncmerIteratorWrapper { return siw; });

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

    nb::class_<StrobemerIndex>(m, "StrobemerIndex")
        .def(nb::init<References&, IndexParameters&>())
        .def("populate", &StrobemerIndex::populate, "f"_a = 0.0002, "threads"_a = 1)
        ;

    m.def("randstrobes_query", &randstrobes_query);

    nb::bind_vector<QueryRandstrobeVector>(m, "QueryRandstrobeVector");

    nb::bind_vector<std::vector<Nam>>(m, "NamVector");

    nb::class_<Nam>(m, "Nam")
        .def_ro("query_start", &Nam::query_s)
        .def_ro("query_end", &Nam::query_e)
        .def_ro("ref_start", &Nam::ref_s)
        .def_ro("ref_end", &Nam::ref_e)
        .def_ro("score", &Nam::score)
        .def("__repr__", [](const Nam& nam) { std::stringstream s; s << nam; return s.str(); })
    ;
    m.def("find_nams", [](const QueryRandstrobeVector &query_randstrobes, const StrobemerIndex& index) -> std::vector<Nam> {
        return find_nams(query_randstrobes, index).second;
    });
}
