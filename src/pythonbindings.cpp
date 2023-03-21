#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/make_iterator.h>
#include <iostream>
#include "randstrobes.hpp"
#include "refs.hpp"

namespace nb = nanobind;
//using namespace nb::literals;

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
        ;

    nb::class_<Record>(m, "Record")
        .def_prop_ro("name", &Record::name)
        .def_prop_ro("sequence", &Record::sequence)
        ;
}
