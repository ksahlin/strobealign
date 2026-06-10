use std::cmp::min;
use std::collections::VecDeque;
use std::marker::PhantomData;

use super::InvalidSeedingParameter;
use super::hash::xxh64;

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

/// Encoding strategy for the [`SyncmerIterator`]. An implementation defines how
/// nucleotides are accumulated into k-mer and s-mer values, how the canonical
/// (forward vs. reverse) orientation is chosen, and how a finished syncmer is
/// produced. [`KmerEncoding`] is the standard nucleotide encoding; further
/// encodings (e.g. for aDNA) plug in by implementing this trait.
pub trait SyncmerEncoding {
    type Syncmer: Copy;
    type KmerValue: Copy + Default;
    type SmerValue: Copy + Default;

    fn encode_nucleotide(ch: u8) -> Option<u8>;
    fn kmer_mask(k: usize) -> Self::KmerValue;
    fn kmer_shift(k: usize) -> usize;
    fn smer_mask(s: usize) -> Self::SmerValue;
    fn smer_shift(s: usize) -> usize;
    fn update_kmer_forward(
        current: Self::KmerValue,
        encoded: u8,
        mask: Self::KmerValue,
    ) -> Self::KmerValue;
    fn update_kmer_reverse(current: Self::KmerValue, encoded: u8, shift: usize) -> Self::KmerValue;
    fn update_smer_forward(
        current: Self::SmerValue,
        encoded: u8,
        mask: Self::SmerValue,
    ) -> Self::SmerValue;
    fn update_smer_reverse(current: Self::SmerValue, encoded: u8, shift: usize) -> Self::SmerValue;
    fn canonical_smer(forward: Self::SmerValue, reverse: Self::SmerValue) -> Self::SmerValue;
    /// Returns `(canonical_value, is_forward)` where `is_forward = (forward <= reverse)`.
    fn canonical_kmer(
        forward: Self::KmerValue,
        reverse: Self::KmerValue,
    ) -> (Self::KmerValue, bool);
    fn hash_smer(value: Self::SmerValue) -> u64;
    fn make_syncmer(
        kmer: Self::KmerValue,
        position: usize,
        k: usize,
        ry_len: usize,
        is_forward: bool,
    ) -> Self::Syncmer;
}

pub struct KmerEncoding;

impl SyncmerEncoding for KmerEncoding {
    type Syncmer = Syncmer;
    type KmerValue = u64;
    type SmerValue = u64;

    #[inline]
    fn encode_nucleotide(ch: u8) -> Option<u8> {
        let c = NUCLEOTIDES[ch as usize];
        if c < 4 { Some(c) } else { None }
    }

    #[inline]
    fn kmer_mask(k: usize) -> u64 {
        (1u64 << (2 * k)) - 1
    }

    #[inline]
    fn kmer_shift(k: usize) -> usize {
        (k - 1) * 2
    }

    #[inline]
    fn smer_mask(s: usize) -> u64 {
        (1u64 << (2 * s)) - 1
    }

    #[inline]
    fn smer_shift(s: usize) -> usize {
        (s - 1) * 2
    }

    #[inline]
    fn update_kmer_forward(current: u64, encoded: u8, mask: u64) -> u64 {
        ((current << 2) | (encoded as u64)) & mask
    }

    #[inline]
    fn update_kmer_reverse(current: u64, encoded: u8, shift: usize) -> u64 {
        (current >> 2) | (((3 - encoded) as u64) << shift)
    }

    #[inline]
    fn update_smer_forward(current: u64, encoded: u8, mask: u64) -> u64 {
        ((current << 2) | (encoded as u64)) & mask
    }

    #[inline]
    fn update_smer_reverse(current: u64, encoded: u8, shift: usize) -> u64 {
        (current >> 2) | (((3 - encoded) as u64) << shift)
    }

    #[inline]
    fn canonical_smer(forward: u64, reverse: u64) -> u64 {
        min(forward, reverse)
    }

    #[inline]
    fn canonical_kmer(forward: u64, reverse: u64) -> (u64, bool) {
        (min(forward, reverse), forward <= reverse)
    }

    #[inline]
    fn hash_smer(value: u64) -> u64 {
        xxh64(value)
    }

    #[inline]
    fn make_syncmer(
        kmer: u64,
        position: usize,
        _k: usize,
        _ry_len: usize,
        is_forward: bool,
    ) -> Syncmer {
        Syncmer {
            hash: xxh64(kmer),
            is_forward,
            position,
        }
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
        if k < 8 || k > 32 {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "k must be at least 8 and at most 32",
            ));
        }
        if s > k {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "s must not be larger than k",
            ));
        }
        if (k - s) % 2 != 0 {
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

pub struct SyncmerIterator<'a, E: SyncmerEncoding> {
    seq: &'a [u8],
    k: usize,
    s: usize,
    t: usize,
    ry_len: usize,
    kmask: E::KmerValue,
    smask: E::SmerValue,
    kshift: usize,
    sshift: usize,
    qs: VecDeque<u64>, // s-mer hashes
    qs_min_val: u64,
    l: usize,
    xk: [E::KmerValue; 2],
    xs: [E::SmerValue; 2],
    i: usize,
    exhausted: bool,
    _phantom: PhantomData<E>,
}

impl<'a, E: SyncmerEncoding> SyncmerIterator<'a, E> {
    pub fn new(seq: &'a [u8], k: usize, s: usize, t: usize) -> Self {
        Self::with_ry_len(seq, k, s, t, k)
    }

    pub fn with_ry_len(seq: &'a [u8], k: usize, s: usize, t: usize, ry_len: usize) -> Self {
        SyncmerIterator {
            seq,
            k,
            s,
            t,
            ry_len,
            kmask: E::kmer_mask(k),
            smask: E::smer_mask(s),
            kshift: E::kmer_shift(k),
            sshift: E::smer_shift(s),
            qs: VecDeque::new(),
            qs_min_val: u64::MAX,
            l: 0,
            xk: [E::KmerValue::default(), E::KmerValue::default()],
            xs: [E::SmerValue::default(), E::SmerValue::default()],
            i: 0,
            exhausted: false,
            _phantom: PhantomData,
        }
    }
}

pub type KmerSyncmerIterator<'a> = SyncmerIterator<'a, KmerEncoding>;

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

impl<E: SyncmerEncoding> Iterator for SyncmerIterator<'_, E> {
    type Item = E::Syncmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        for i in self.i..self.seq.len() {
            let ch = self.seq[i];
            if let Some(c) = E::encode_nucleotide(ch) {
                // not an "N" base
                self.xk[0] = E::update_kmer_forward(self.xk[0], c, self.kmask); // forward strand
                self.xk[1] = E::update_kmer_reverse(self.xk[1], c, self.kshift); // reverse strand
                self.xs[0] = E::update_smer_forward(self.xs[0], c, self.smask); // forward strand
                self.xs[1] = E::update_smer_reverse(self.xs[1], c, self.sshift); // reverse strand
                self.l += 1;
                if self.l < self.s {
                    continue;
                }
                // we find an s-mer
                let ys = E::canonical_smer(self.xs[0], self.xs[1]);
                let hash_s = E::hash_smer(ys);
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
                    let (yk, is_forward) = E::canonical_kmer(self.xk[0], self.xk[1]);
                    let syncmer =
                        E::make_syncmer(yk, i + 1 - self.k, self.k, self.ry_len, is_forward);
                    self.i = i + 1;
                    return Some(syncmer);
                }
            } else {
                // if there is an "N", restart
                self.qs_min_val = u64::MAX;
                self.l = 0;
                self.xs = [E::SmerValue::default(), E::SmerValue::default()];
                self.xk = [E::KmerValue::default(), E::KmerValue::default()];
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
    fn test_invalid_parameters() {
        // k-s not even
        assert!(SyncmerParameters::try_new(20, 17).is_err());
        // s larger than k
        assert!(SyncmerParameters::try_new(20, 25).is_err());
        // k out of range
        assert!(SyncmerParameters::try_new(34, 16).is_err());
        assert!(SyncmerParameters::try_new(6, 25).is_err());
    }

    #[test]
    fn test_syncmer_iterator() {
        let seq = "AAAAAAAAAAAAAAAAAAAA";
        assert_eq!(seq.len(), 20);

        let parameters = SyncmerParameters::try_new(8, 4).unwrap();
        assert_eq!(parameters.t, 3);

        let mut iterator =
            KmerSyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let syncmer = iterator.next().unwrap();
        assert_eq!(syncmer.position, 0);
    }

    fn syncmers_of(seq: &[u8], parameters: &SyncmerParameters) -> Vec<Syncmer> {
        KmerSyncmerIterator::new(seq, parameters.k, parameters.s, parameters.t).collect()
    }

    #[test]
    fn test_canonical_syncmers() {
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        let seq = read_ref("tests/phix.fasta")
            .unwrap()
            .pop()
            .unwrap()
            .sequence;
        let sequences = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".into(), seq];
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
