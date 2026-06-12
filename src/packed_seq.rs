/// 2-bit packed nucleotide sequence.
///
/// Encoding: A=0, C=1, G=2, T/U=3.  Ambiguous/unknown bases (N, R, Y, ...)
/// are replaced by a deterministic pseudo-random 2-bit value derived from
/// their position (following BWA's convention), so that syncmers spanning
/// an ambiguous region have effectively random hashes and will not
/// systematically collide with real read k-mers.
///
/// 32 bases are packed per u64, little-endian within each word
/// (base i occupies bits 2*(i%32)+1 : 2*(i%32) of word i/32).
#[derive(Debug, Clone, Default)]
pub struct PackedSeq {
    data: Vec<u64>,
    len: usize,
}

const fn make_encode_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}

/// IS_WILDCARD[c] = 1 for any ASCII byte c that is not a standard ACGTU nucleotide.
const fn make_is_wildcard_table() -> [u8; 256] {
    let mut t = [1u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 0;
    t[b'c' as usize] = 0;
    t[b'G' as usize] = 0;
    t[b'g' as usize] = 0;
    t[b'T' as usize] = 0;
    t[b't' as usize] = 0;
    t[b'U' as usize] = 0;
    t[b'u' as usize] = 0;
    t
}

static ENCODE: [u8; 256] = make_encode_table();
static IS_WILDCARD: [u8; 256] = make_is_wildcard_table();
static DECODE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Deterministic pseudo-random 2-bit value for position `pos`.
///
/// This uses a Fibonacci/Knuth multiplicative hash so every ambiguous base at a
/// distinct position maps to a random-looking 2-bit value, making syncmers
/// that span N runs have effectively random hashes (following BWA's approach).
#[inline]
fn encode_wildcard(pos: usize) -> u64 {
    let h = (pos as u64).wrapping_mul(0x9e3779b97f4a7c15);
    (h ^ (h >> 30)) & 3
}

impl PackedSeq {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a `PackedSeq` by packing every byte of a slice.
    pub fn from_slice(s: &[u8]) -> Self {
        let mut seq = PackedSeq::new();
        for &c in s {
            seq.push(c);
        }
        seq
    }

    /// Append one byte (upper- or lowercase ASCII nucleotide).
    /// ACGTU → packed 2-bit value; any other byte (N, R, Y, ...) →
    /// deterministic pseudo-random 2-bit value derived from its position,
    /// so that syncmers spanning ambiguous regions get random hashes
    pub fn push(&mut self, c: u8) {
        let bits = if IS_WILDCARD[c as usize] == 0 {
            ENCODE[c as usize] as u64
        } else {
            encode_wildcard(self.len)
        };
        let bit = (self.len & 31) * 2;
        if bit == 0 {
            self.data.push(bits);
        } else {
            *self.data.last_mut().unwrap() |= bits << bit;
        }
        self.len += 1;
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Return the raw 2-bit value of base `i` (0=A, 1=C, 2=G, 3=T).
    /// Ambiguous positions return a pseudo-random value in 0..=3 (never 4).
    /// No table lookups - use this in inner loops.
    pub fn nucleotide_bits(&self, i: usize) -> u8 {
        let word = i >> 5;
        let bit = (i & 31) * 2;
        ((self.data[word] >> bit) & 3) as u8
    }

    /// Decode base `i` as an ASCII byte (A/C/G/T).
    /// Ambiguous positions decode to a pseudo-random nucleotide, not 'N'.
    pub fn decode_at(&self, i: usize) -> u8 {
        DECODE[self.nucleotide_bits(i) as usize]
    }

    /// Decode bases `start..end` as ASCII bytes.
    pub fn decode(&self, start: usize, end: usize) -> Vec<u8> {
        (start..end).map(|i| self.decode_at(i)).collect()
    }

    /// Decode the full sequence as ASCII bytes.
    pub fn decode_all(&self) -> Vec<u8> {
        self.decode(0, self.len)
    }
}
