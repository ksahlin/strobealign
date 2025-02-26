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


#[cfg(test)]
mod test {
    use crate::revcomp::reverse_complement;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b""), b"");
        assert_eq!(reverse_complement(b""), b"");
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"C"), b"G");
        assert_eq!(reverse_complement(b"G"), b"C");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"TG"), b"CA");
        assert_eq!(reverse_complement(b"AC"), b"GT");
        assert_eq!(reverse_complement(b"ACG"), b"CGT");
        assert_eq!(reverse_complement(b"AACGT"), b"ACGTT");
    }
}
