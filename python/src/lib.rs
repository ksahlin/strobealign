use pyo3;

#[pyo3::pymodule]
mod strobealign {
    use pyo3::prelude::*;

    use ::strobealign::revcomp;

    #[pyfunction]
    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        revcomp::reverse_complement(seq)
    }
}
