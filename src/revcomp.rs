// a, A -> T
// c, C -> G
// g, G -> C
// t, T, u, U -> A
const REVCOMP_TABLE: [u8; 256] = [
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'T', b'N', b'G',  b'N', b'N', b'N', b'C',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'A', b'A', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'T', b'N', b'G',  b'N', b'N', b'N', b'C',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'A', b'A', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',
    b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N',  b'N', b'N', b'N', b'N'
];

pub fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|ch| REVCOMP_TABLE[*ch as usize]).collect()
}
