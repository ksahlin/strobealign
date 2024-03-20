use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::time::Instant;
use pyo3::prelude::*;
use pyo3::types::PyString;
use ::rstrobes::fasta;
use ::rstrobes::fasta::RefSequence;

#[pyclass(name = "RefSequence")]
struct PyRefSequence {
    refseq: RefSequence,
}

#[pyclass(name = "References")]
struct PyReferences {
    references: Vec<RefSequence>,
}

#[pymethods]
impl PyReferences {
    #[staticmethod]
    fn from_fasta(path: &str) -> PyResult<PyReferences> {
        let f = File::open(path)?;
        let mut reader = BufReader::new(f);
        let references = fasta::read_fasta(&mut reader).unwrap();

        Ok(Self { references })
    }
    
    fn __len__(&self) -> usize {
        self.references.len()
    }
}

/*impl PyRefSequence {
    fn new(refseq: RefSequence) -> PyRefSequence {
        Self { refseq }
    }
}*/
/*
#[pymethods]
impl PyRefSequence {
    #[new]
    fn new(refseq: RefSequence) -> PyRefSequence {
        Self { refseq }
    }
    
    /*.def(nb::init())
    .def("add", &References::add)
    .def_static("from_fasta", &References::from_fasta)
    .def("__getitem__", [](const References& refs, size_t i) {
    return Record(refs.names[i], refs.sequences[i]);
    })
    .def("__len__", [](const References& refs) { return refs.sequences.size(); })
    */
    
}

#[pyfunction]
fn read_fasta<'py>(path: &Bound<'py, PyString>) -> PyResult<Bound<'py, PyReferences>> {
    let f = File::open(path.to_string())?;
    let mut reader = BufReader::new(f);
    let references = fasta::read_fasta(&mut reader).unwrap();

    Ok(PyReferences { references })
}

fn hello() -> PyResult<String> {
    Ok("Hello".to_string())
}
*/
#[pymodule]
fn rstrobes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyReferences>()?;
    //m.add_function(wrap_pyfunction!(hello, m)?)?;
    Ok(())
}
