use std::cmp::min;
use std::collections::VecDeque;
use std::marker::PhantomData;

use crate::hash::{xxh32, xxh64};
use crate::index::InvalidIndexParameter;

pub trait SyncmerLike: Copy {
    fn hash(&self) -> u64;
    fn position(&self) -> usize;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Syncmer {
    pub hash: u64,
    pub position: usize,
}

impl SyncmerLike for Syncmer {
    #[inline]
    fn hash(&self) -> u64 {
        self.hash
    }

    #[inline]
    fn position(&self) -> usize {
        self.position
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RymerSyncmer {
    pub hash1: u32,
    pub hash2: u32,
    pub position: usize,
}

impl SyncmerLike for RymerSyncmer {
    #[inline]
    fn hash(&self) -> u64 {
        ((self.hash1 as u64) << 32) | (self.hash2 as u64)
    }

    #[inline]
    fn position(&self) -> usize {
        self.position
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SyncmerParameters {
    pub k: usize,
    pub s: usize,
    pub t: usize,
}

impl SyncmerParameters {
    pub fn try_new(k: usize, s: usize) -> Result<Self, InvalidIndexParameter> {
        if k < 8 || k > 32 {
            return Err(InvalidIndexParameter::InvalidParameter(
                "k must be at least 8 and at most 32",
            ));
        }
        if s > k {
            return Err(InvalidIndexParameter::InvalidParameter(
                "s must not be larger than k",
            ));
        }
        if (k - s) % 2 != 0 {
            return Err(InvalidIndexParameter::InvalidParameter(
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

pub trait SyncmerEncoding {
    type Syncmer: SyncmerLike;
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
    fn canonical_kmer(forward: Self::KmerValue, reverse: Self::KmerValue) -> Self::KmerValue;
    fn hash_smer(value: Self::SmerValue) -> u64;
    fn make_syncmer(kmer: Self::KmerValue, position: usize) -> Self::Syncmer;
}

pub struct KmerEncoding;

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
    fn canonical_kmer(forward: u64, reverse: u64) -> u64 {
        min(forward, reverse)
    }

    #[inline]
    fn hash_smer(value: u64) -> u64 {
        xxh64(value)
    }

    #[inline]
    fn make_syncmer(kmer: u64, position: usize) -> Syncmer {
        Syncmer {
            hash: xxh64(kmer),
            position,
        }
    }
}

pub struct RymerEncoding;

// S1: a, A -> 0; c, C -> 1; g, G -> 0; t, T, u, U -> 1
#[rustfmt::skip]
static RYMER_S1: [u8; 256] = [
        0, 1, 0, 1,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 0,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  1, 1, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 0,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  1, 1, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
];

// S2: a, A -> 0; c, C -> 0; g, G -> 1; t, T, u, U -> 1
#[rustfmt::skip]
static RYMER_S2: [u8; 256] = [
        0, 0, 1, 1,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 0,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  1, 1, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 0,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  1, 1, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
];

impl SyncmerEncoding for RymerEncoding {
    type Syncmer = RymerSyncmer;
    type KmerValue = (u32, u32);
    type SmerValue = (u32, u32);

    #[inline]
    fn encode_nucleotide(ch: u8) -> Option<u8> {
        let s1 = RYMER_S1[ch as usize];
        if s1 >= 4 {
            return None;
        }
        let s2 = RYMER_S2[ch as usize];
        Some(s1 | (s2 << 1))
    }

    #[inline]
    fn kmer_mask(k: usize) -> (u32, u32) {
        let mask = (1u32 << k) - 1;
        (mask, mask)
    }

    #[inline]
    fn kmer_shift(k: usize) -> usize {
        k - 1
    }

    #[inline]
    fn smer_mask(s: usize) -> (u32, u32) {
        let mask = (1u32 << s) - 1;
        (mask, mask)
    }

    #[inline]
    fn smer_shift(s: usize) -> usize {
        s - 1
    }

    #[inline]
    fn update_kmer_forward(current: (u32, u32), encoded: u8, mask: (u32, u32)) -> (u32, u32) {
        let s1 = encoded & 1;
        let s2 = (encoded >> 1) & 1;
        (
            ((current.0 << 1) | (s1 as u32)) & mask.0,
            ((current.1 << 1) | (s2 as u32)) & mask.1,
        )
    }

    #[inline]
    fn update_kmer_reverse(current: (u32, u32), encoded: u8, shift: usize) -> (u32, u32) {
        let s1 = encoded & 1;
        let s2 = (encoded >> 1) & 1;
        let s1_comp = 1 - s1;
        let s2_comp = 1 - s2;
        (
            (current.0 >> 1) | ((s1_comp as u32) << shift),
            (current.1 >> 1) | ((s2_comp as u32) << shift),
        )
    }

    #[inline]
    fn update_smer_forward(current: (u32, u32), encoded: u8, mask: (u32, u32)) -> (u32, u32) {
        let s1 = encoded & 1;
        let s2 = (encoded >> 1) & 1;
        (
            ((current.0 << 1) | (s1 as u32)) & mask.0,
            ((current.1 << 1) | (s2 as u32)) & mask.1,
        )
    }

    #[inline]
    fn update_smer_reverse(current: (u32, u32), encoded: u8, shift: usize) -> (u32, u32) {
        let s1 = encoded & 1;
        let s2 = (encoded >> 1) & 1;
        let s1_comp = 1 - s1;
        let s2_comp = 1 - s2;
        (
            (current.0 >> 1) | ((s1_comp as u32) << shift),
            (current.1 >> 1) | ((s2_comp as u32) << shift),
        )
    }

    #[inline]
    fn canonical_smer(forward: (u32, u32), reverse: (u32, u32)) -> (u32, u32) {
        let fwd_combined = ((forward.0 as u64) << 32) | (forward.1 as u64);
        let rev_combined = ((reverse.0 as u64) << 32) | (reverse.1 as u64);
        if fwd_combined <= rev_combined {
            forward
        } else {
            reverse
        }
    }

    #[inline]
    fn canonical_kmer(forward: (u32, u32), reverse: (u32, u32)) -> (u32, u32) {
        let fwd_combined = ((forward.0 as u64) << 32) | (forward.1 as u64);
        let rev_combined = ((reverse.0 as u64) << 32) | (reverse.1 as u64);
        if fwd_combined <= rev_combined {
            forward
        } else {
            reverse
        }
    }

    #[inline]
    fn hash_smer(value: (u32, u32)) -> u64 {
        xxh32(value.0) as u64
    }

    #[inline]
    fn make_syncmer(kmer: (u32, u32), position: usize) -> RymerSyncmer {
        RymerSyncmer {
            hash1: xxh32(kmer.0),
            hash2: xxh32(kmer.1),
            position,
        }
    }
}

pub struct SyncmerIterator<'a, E: SyncmerEncoding> {
    seq: &'a [u8],
    k: usize,
    s: usize,
    t: usize,
    kmask: E::KmerValue,
    smask: E::SmerValue,
    kshift: usize,
    sshift: usize,
    qs: VecDeque<u64>,
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
        SyncmerIterator {
            seq,
            k,
            s,
            t,
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
                self.xk[0] = E::update_kmer_forward(self.xk[0], c, self.kmask);
                self.xk[1] = E::update_kmer_reverse(self.xk[1], c, self.kshift);
                self.xs[0] = E::update_smer_forward(self.xs[0], c, self.smask);
                self.xs[1] = E::update_smer_reverse(self.xs[1], c, self.sshift);
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
                    let yk = E::canonical_kmer(self.xk[0], self.xk[1]);
                    let syncmer = E::make_syncmer(yk, i + 1 - self.k);
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

pub type KmerSyncmerIterator<'a> = SyncmerIterator<'a, KmerEncoding>;
pub type RymerSyncmerIterator<'a> = SyncmerIterator<'a, RymerEncoding>;

#[cfg(test)]
mod test {
    use super::{
        KmerSyncmerIterator, NUCLEOTIDES, RYMER_S1, RYMER_S2, RymerSyncmerIterator,
        SyncmerIterator, SyncmerLike, SyncmerParameters,
    };

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

        let mut iterator: KmerSyncmerIterator =
            SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let syncmer = iterator.next().unwrap();
        assert_eq!(syncmer.position(), 0);
    }

    #[test]
    fn test_rymer_encoding_tables() {
        // A -> (0, 0)
        assert_eq!(RYMER_S1[b'A' as usize], 0);
        assert_eq!(RYMER_S2[b'A' as usize], 0);
        assert_eq!(RYMER_S1[b'a' as usize], 0);
        assert_eq!(RYMER_S2[b'a' as usize], 0);

        // C -> (1, 0)
        assert_eq!(RYMER_S1[b'C' as usize], 1);
        assert_eq!(RYMER_S2[b'C' as usize], 0);
        assert_eq!(RYMER_S1[b'c' as usize], 1);
        assert_eq!(RYMER_S2[b'c' as usize], 0);

        // G -> (0, 1)
        assert_eq!(RYMER_S1[b'G' as usize], 0);
        assert_eq!(RYMER_S2[b'G' as usize], 1);
        assert_eq!(RYMER_S1[b'g' as usize], 0);
        assert_eq!(RYMER_S2[b'g' as usize], 1);

        // T -> (1, 1)
        assert_eq!(RYMER_S1[b'T' as usize], 1);
        assert_eq!(RYMER_S2[b'T' as usize], 1);
        assert_eq!(RYMER_S1[b't' as usize], 1);
        assert_eq!(RYMER_S2[b't' as usize], 1);

        // N -> invalid (4)
        assert_eq!(RYMER_S1[b'N' as usize], 4);
        assert_eq!(RYMER_S2[b'N' as usize], 4);
    }

    #[test]
    fn test_kmer_encoding_table() {
        // A -> 0
        assert_eq!(NUCLEOTIDES[b'A' as usize], 0);
        assert_eq!(NUCLEOTIDES[b'a' as usize], 0);
        // C -> 1
        assert_eq!(NUCLEOTIDES[b'C' as usize], 1);
        assert_eq!(NUCLEOTIDES[b'c' as usize], 1);
        // G -> 2
        assert_eq!(NUCLEOTIDES[b'G' as usize], 2);
        assert_eq!(NUCLEOTIDES[b'g' as usize], 2);
        // T -> 3
        assert_eq!(NUCLEOTIDES[b'T' as usize], 3);
        assert_eq!(NUCLEOTIDES[b't' as usize], 3);
        // N -> invalid
        assert_eq!(NUCLEOTIDES[b'N' as usize], 4);
    }

    #[test]
    fn test_rymer_syncmer_iterator() {
        let seq = "ACGTACGTACGTACGTACGT";
        assert_eq!(seq.len(), 20);

        let parameters = SyncmerParameters::try_new(8, 4).unwrap();

        let mut iterator: RymerSyncmerIterator =
            SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let syncmer = iterator.next();
        assert!(syncmer.is_some());
    }

    #[test]
    fn test_syncmer_like_trait() {
        let seq = "ACGTACGTACGTACGTACGT";
        let parameters = SyncmerParameters::try_new(8, 4).unwrap();

        let mut kmer_iter: KmerSyncmerIterator =
            SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let kmer_syncmer = kmer_iter.next().unwrap();
        assert!(kmer_syncmer.hash() > 0);
        assert!(kmer_syncmer.position() < seq.len());

        let mut rymer_iter: RymerSyncmerIterator =
            SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let rymer_syncmer = rymer_iter.next().unwrap();
        assert!(rymer_syncmer.hash() > 0);
        assert!(rymer_syncmer.position() < seq.len());
    }
}
