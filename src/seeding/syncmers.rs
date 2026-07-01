use std::cmp::min;
use std::collections::VecDeque;

use super::InvalidSeedingParameter;
use super::hash::xxh64;
use crate::packed_seq::{PackedSeq, PackedSeqSlice};

/// Trait for types that can serve as a sequence source for `SyncmerIterator`.
///
/// `nucleotide_bits` returns the nucleotide as a 2-bit value (0=A, 1=C, 2=G, 3=T,
/// 4=N/ambiguous). This is what the syncmer hash machinery needs directly,
/// so implementations should avoid going through ASCII as an intermediate.
#[allow(clippy::len_without_is_empty)]
pub trait SeqAccess {
    fn nucleotide_bits(&self, i: usize) -> u8;
    fn len(&self) -> usize;
}

impl SeqAccess for &[u8] {
    /// ASCII byte → 2-bit (0-3) or 4 for N.
    fn nucleotide_bits(&self, i: usize) -> u8 {
        NUCLEOTIDES[self[i] as usize]
    }

    fn len(&self) -> usize {
        (**self).len()
    }
}

impl<'a> SeqAccess for &'a PackedSeqSlice<'a> {
    /// Direct 2-bit extract - no table lookups, never returns 4.
    /// `**self` reaches `PackedSeq` so Rust resolves the inherent method,
    /// not this trait method (no recursion).
    fn nucleotide_bits(&self, i: usize) -> u8 {
        (**self).nucleotide_bits(i)
    }

    fn len(&self) -> usize {
        (**self).len()
    }
}

impl SeqAccess for &PackedSeq {
    fn nucleotide_bits(&self, i: usize) -> u8 {
        (**self).nucleotide_bits(i)
    }

    fn len(&self) -> usize {
        (**self).len()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Syncmer {
    hash: u64,
    is_forward: bool,
    pub position: usize,
}

impl Syncmer {
    pub fn hash(&self) -> u64 {
        self.hash
    }

    pub fn is_forward(&self) -> bool {
        self.is_forward
    }

    pub fn toggle_orientation(&mut self) {
        self.is_forward = !self.is_forward;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SyncmerParameters {
    pub k: usize,
    pub s: usize,
    pub t: usize,
}

impl SyncmerParameters {
    pub fn try_new(k: usize, s: usize) -> Result<Self, InvalidSeedingParameter> {
        if !(8..=32).contains(&k) {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "k must be at least 8 and at most 32",
            ));
        }
        if s > k {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "s must not be larger than k",
            ));
        }
        if !(k - s).is_multiple_of(2) {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "(k - s) must be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...",
            ));
        }
        Ok(SyncmerParameters {
            k,
            s,
            t: (k - s) / 2 + 1,
        })
    }
}

pub struct SyncmerIterator<S: SeqAccess> {
    seq: S,
    k: usize,
    s: usize,
    t: usize,
    kmask: u64,
    smask: u64,
    kshift: usize,
    sshift: usize,
    qs: VecDeque<u64>, // s-mer hashes
    qs_min_val: u64,
    l: usize,
    xk: [u64; 2],
    xs: [u64; 2],
    i: usize,
    exhausted: bool,
}

impl<S: SeqAccess> SyncmerIterator<S> {
    pub fn new(seq: S, k: usize, s: usize, t: usize) -> SyncmerIterator<S> {
        SyncmerIterator {
            seq,
            k,
            s,
            t,
            kmask: (1 << (2 * k)) - 1,
            smask: (1 << (2 * s)) - 1,
            kshift: (k - 1) * 2,
            sshift: (s - 1) * 2,
            qs: VecDeque::new(),
            qs_min_val: u64::MAX,
            l: 0,
            xk: [0, 0],
            xs: [0, 0],
            i: 0,
            exhausted: false,
        }
    }
}

// a, A -> 0
// c, C -> 1
// g, G -> 2
// t, T, u, U -> 3
#[rustfmt::skip]
static NUCLEOTIDES: [u8; 256] = [
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
];

fn syncmer_kmer_hash(value: u64) -> u64 {
    xxh64(value)
}

fn syncmer_smer_hash(value: u64) -> u64 {
    xxh64(value)
}

impl<S: SeqAccess> Iterator for SyncmerIterator<S> {
    type Item = Syncmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        for i in self.i..self.seq.len() {
            let c = self.seq.nucleotide_bits(i);
            if c < 4 {
                // not an "N" base
                self.xk[0] = ((self.xk[0] << 2) | (c as u64)) & self.kmask; // forward strand
                self.xk[1] = (self.xk[1] >> 2) | (((3 - c) as u64) << self.kshift); // reverse strand
                self.xs[0] = ((self.xs[0] << 2) | (c as u64)) & self.smask; // forward strand
                self.xs[1] = (self.xs[1] >> 2) | (((3 - c) as u64) << self.sshift); // reverse strand
                self.l += 1;
                if self.l < self.s {
                    continue;
                }
                // we find an s-mer
                let ys = min(self.xs[0], self.xs[1]);
                let hash_s = syncmer_smer_hash(ys);
                self.qs.push_back(hash_s);
                // not enough hashes in the queue, yet
                if self.qs.len() < self.k - self.s + 1 {
                    continue;
                }
                if self.qs.len() == (self.k - self.s + 1) {
                    // We are at the last s-mer within the first k-mer, need to decide if we add it
                    for j in 0..self.qs.len() {
                        self.qs_min_val = min(self.qs[j], self.qs_min_val);
                    }
                } else {
                    // update queue and current minimum and position
                    let front = self.qs.pop_front().unwrap();
                    if front == self.qs_min_val {
                        // we popped a minimum, find new brute force
                        self.qs_min_val = u64::MAX;
                        for j in 0..self.qs.len() {
                            self.qs_min_val = min(self.qs[j], self.qs_min_val);
                        }
                    } else if hash_s < self.qs_min_val {
                        // the new value added to queue is the new minimum
                        self.qs_min_val = hash_s;
                    }
                }
                if self.qs[self.t - 1] == self.qs_min_val {
                    // occurs at t:th position in k-mer
                    let yk = min(self.xk[0], self.xk[1]);
                    let syncmer = Syncmer {
                        hash: syncmer_kmer_hash(yk),
                        is_forward: self.xk[0] <= self.xk[1],
                        position: i + 1 - self.k,
                    };
                    self.i = i + 1;
                    return Some(syncmer);
                }
            } else {
                // if there is an "N", restart
                self.qs_min_val = u64::MAX;
                self.l = 0;
                self.xs = [0, 0];
                self.xk = [0, 0];
                self.qs.clear();
            }
        }
        self.exhausted = true;

        None
    }
}

#[cfg(test)]
mod test {
    use crate::{io::fasta::read_ref, revcomp::reverse_complement};

    use super::*;

    #[test]
    fn invalid_parameters() {
        // k-s not even
        assert!(SyncmerParameters::try_new(20, 17).is_err());
        // s larger than k
        assert!(SyncmerParameters::try_new(20, 25).is_err());
        // k out of range
        assert!(SyncmerParameters::try_new(34, 16).is_err());
        assert!(SyncmerParameters::try_new(6, 25).is_err());
    }

    #[test]
    fn syncmer_iterator() {
        let seq = "AAAAAAAAAAAAAAAAAAAA";
        assert_eq!(seq.len(), 20);

        let parameters = SyncmerParameters::try_new(8, 4).unwrap();
        assert_eq!(parameters.t, 3);

        let mut iterator =
            SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let syncmer = iterator.next().unwrap();
        assert_eq!(syncmer.position, 0);
    }

    fn syncmers_of(seq: &[u8], parameters: &SyncmerParameters) -> Vec<Syncmer> {
        SyncmerIterator::new(seq, parameters.k, parameters.s, parameters.t).collect()
    }

    #[test]
    fn canonical_syncmers() {
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        let seq: Vec<u8> = read_ref("tests/phix.fasta").unwrap().contig(0).decode_all();
        let sequences: [Vec<u8>; 2] = [b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(), seq];
        for s in sequences {
            let seq_revcomp = reverse_complement(&s);
            let syncmers_forward = syncmers_of(&s, &parameters);
            let mut syncmers_reverse = syncmers_of(&seq_revcomp, &parameters);
            syncmers_reverse.reverse();
            for syncmer_rev in &mut syncmers_reverse {
                syncmer_rev.position = s.len() - parameters.k - syncmer_rev.position;
            }
            assert_eq!(syncmers_forward.len(), syncmers_reverse.len());
            for (sf, sr) in syncmers_forward.iter().zip(syncmers_reverse.iter()) {
                assert_eq!(sf.hash(), sr.hash());
                assert_eq!(sf.position, sr.position);
                assert_ne!(sf.is_forward(), sr.is_forward());
            }
        }
    }
}
