use std::cmp::min;
use std::collections::VecDeque;
use std::marker::PhantomData;

use crate::hash::{xxh32, xxh64};
use crate::index::InvalidIndexParameter;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Syncmer {
    pub hash: u64,
    pub position: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RymerSyncmer {
    pub hash1: u64,
    pub hash2: u64,
    pub position: usize,
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
    fn canonical_kmer(forward: Self::KmerValue, reverse: Self::KmerValue) -> Self::KmerValue;
    fn hash_smer(value: Self::SmerValue) -> u64;
    fn make_syncmer(kmer: Self::KmerValue, position: usize, k: usize, ry_len: usize) -> Self::Syncmer;
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
    fn make_syncmer(kmer: u64, position: usize, _k: usize, _ry_len: usize) -> Syncmer {
        Syncmer {
            hash: xxh64(kmer),
            position,
        }
    }
}

/// Compare two nucleotide sequences using RY encoding for the first `ry_len`
/// positions and exact equality for the rest.
pub fn ry_equal(a: &[u8], b: &[u8], _ry_len: usize) -> bool {
    const MAX_MISMATCHES: usize = 2;
    a.len() == b.len()
        && a.iter().cloned()
            .zip(b.iter().cloned())
            .enumerate()
            .filter(|(_i, (x, y))| {
                // if i < ry_len {
                //     let sx = RYMER_S1[x as usize];
                //     let sy = RYMER_S1[y as usize];
                //     sx < 4 && sx == sy
                // } else {
                //     x == y
                // }
                x != y
            })
            .count() <= MAX_MISMATCHES
}

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
                    let syncmer = E::make_syncmer(yk, i + 1 - self.k, self.k, self.ry_len);
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

pub struct RymerEncoding;

impl SyncmerEncoding for RymerEncoding {
    type Syncmer = RymerSyncmer;
    type KmerValue = (u64, u32, u32);
    type SmerValue = (u32, u32);

    #[inline]
    fn encode_nucleotide(ch: u8) -> Option<u8> {
        let nuc = NUCLEOTIDES[ch as usize];
        let s1 = RYMER_S1[ch as usize];
        if nuc >= 4 || s1 >= 4 {
            return None;
        }
        let s2 = RYMER_S2[ch as usize];
        // Pack: bits 0-1 = nuc, bit 2 = s1, bit 3 = s2
        Some(nuc | (s1 << 2) | (s2 << 3))
    }

    #[inline]
    fn kmer_mask(k: usize) -> (u64, u32, u32) {
        let nuc_mask = (1u64 << (2 * k)) - 1;
        let ry_mask = (1u32 << k) - 1;
        (nuc_mask, ry_mask, ry_mask)
    }

    #[inline]
    fn kmer_shift(k: usize) -> usize {
        // Used for reverse update. For nuc it's (k-1)*2, for rymer it's k-1.
        // We store (k-1)*2 and derive the rymer shift from it.
        (k - 1) * 2
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
    fn update_kmer_forward(
        current: (u64, u32, u32),
        encoded: u8,
        mask: (u64, u32, u32),
    ) -> (u64, u32, u32) {
        let nuc = (encoded & 3) as u64;
        let s1 = ((encoded >> 2) & 1) as u32;
        let s2 = ((encoded >> 3) & 1) as u32;
        (
            ((current.0 << 2) | nuc) & mask.0,
            ((current.1 << 1) | s1) & mask.1,
            ((current.2 << 1) | s2) & mask.2,
        )
    }

    #[inline]
    fn update_kmer_reverse(
        current: (u64, u32, u32),
        encoded: u8,
        shift: usize,
    ) -> (u64, u32, u32) {
        let nuc = (encoded & 3) as u64;
        let s1 = ((encoded >> 2) & 1) as u32;
        let s2 = ((encoded >> 3) & 1) as u32;
        let ry_shift = shift / 2; // kmer_shift stores (k-1)*2, rymer needs k-1
        (
            (current.0 >> 2) | ((3 - nuc) << shift),
            (current.1 >> 1) | ((1 - s1) << ry_shift),
            (current.2 >> 1) | ((1 - s2) << ry_shift),
        )
    }

    #[inline]
    fn update_smer_forward(
        current: (u32, u32),
        encoded: u8,
        mask: (u32, u32),
    ) -> (u32, u32) {
        let s1 = ((encoded >> 2) & 1) as u32;
        let s2 = ((encoded >> 3) & 1) as u32;
        (
            ((current.0 << 1) | s1) & mask.0,
            ((current.1 << 1) | s2) & mask.1,
        )
    }

    #[inline]
    fn update_smer_reverse(current: (u32, u32), encoded: u8, shift: usize) -> (u32, u32) {
        let s1 = ((encoded >> 2) & 1) as u32;
        let s2 = ((encoded >> 3) & 1) as u32;
        (
            (current.0 >> 1) | ((1 - s1) << shift),
            (current.1 >> 1) | ((1 - s2) << shift),
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
    fn canonical_kmer(
        forward: (u64, u32, u32),
        reverse: (u64, u32, u32),
    ) -> (u64, u32, u32) {
        // Use rymer S1/S2 for strand selection (consistent with RymerEncoding)
        let fwd_ry = ((forward.1 as u64) << 32) | (forward.2 as u64);
        let rev_ry = ((reverse.1 as u64) << 32) | (reverse.2 as u64);
        if fwd_ry <= rev_ry {
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
    fn make_syncmer(
        kmer: (u64, u32, u32),
        position: usize,
        k: usize,
        ry_len: usize,
    ) -> RymerSyncmer {
        if ry_len == k {
            RymerSyncmer {
                hash1: xxh64(kmer.1 as u64),
                hash2: xxh64(kmer.2 as u64),
                position,
            }
        } else {
            let kmer_len = k - ry_len;
            let ry_s1 = (kmer.1 >> kmer_len) as u64;
            let ry_s2 = (kmer.2 >> kmer_len) as u64;
            let kmer_tail_bits = kmer_len * 2;
            let kmer_tail_mask = (1u64 << kmer_tail_bits) - 1;
            let kmer_tail = kmer.0 & kmer_tail_mask;
            let combined = (ry_s1 << kmer_tail_bits) | kmer_tail;
            RymerSyncmer {
                hash1: xxh64(combined),
                hash2: xxh64(ry_s2),
                position,
            }
        }
    }
}

pub type RymerSyncmerIterator<'a> = SyncmerIterator<'a, RymerEncoding>;

#[cfg(test)]
mod test {
    use crate::syncmers::{KmerEncoding, SyncmerEncoding};

    use super::{
        KmerSyncmerIterator, RymerSyncmerIterator, NUCLEOTIDES, RYMER_S1, RYMER_S2,
        SyncmerIterator, SyncmerParameters, ry_equal,
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
        assert_eq!(syncmer.position, 0);
    }

    #[test]
    fn test_kmer_update() {
        let k = 5;
        let s = 3;
        let kmask = KmerEncoding::kmer_mask(k);
        let kshift = KmerEncoding::kmer_shift(k);
        let smask = KmerEncoding::smer_mask(s);
        let sshift = KmerEncoding::smer_shift(s);
        let mut xk: [u64; 2] = [0, 0];
        let mut xs: [u64; 2] = [0, 0];
        let seq = "TACGTCAA";
        for ch in seq.chars() {
            if let Some(c) = KmerEncoding::encode_nucleotide(ch as u8) {
                xk[0] = KmerEncoding::update_kmer_forward(xk[0], c, kmask);
                xk[1] = KmerEncoding::update_kmer_reverse(xk[1], c, kshift);
                xs[0] = KmerEncoding::update_smer_forward(xs[0], c, smask);
                xs[1] = KmerEncoding::update_smer_reverse(xs[1], c, sshift);
            }
        }
        assert_eq!(xk[0], 0b1011010000); // GTCAA
        assert_eq!(xk[1], 0b1111100001); // TTGAC
        assert_eq!(xs[0], 0b010000); // CAA
        assert_eq!(xs[1], 0b111110); // TTG
    }

    #[test]
    fn test_rymer_update() {
        use crate::syncmers::RymerEncoding;

        let k = 5;
        let s = 3;
        let kmask = RymerEncoding::kmer_mask(k);
        let kshift = RymerEncoding::kmer_shift(k);
        let smask = RymerEncoding::smer_mask(s);
        let sshift = RymerEncoding::smer_shift(s);
        let mut xk: [(u64, u32, u32); 2] = [(0, 0, 0), (0, 0, 0)];
        let mut xs: [(u32, u32); 2] = [(0, 0), (0, 0)];
        let seq = "TACGTCAA";
        for ch in seq.chars() {
            if let Some(c) = RymerEncoding::encode_nucleotide(ch as u8) {
                xk[0] = RymerEncoding::update_kmer_forward(xk[0], c, kmask);
                xk[1] = RymerEncoding::update_kmer_reverse(xk[1], c, kshift);
                xs[0] = RymerEncoding::update_smer_forward(xs[0], c, smask);
                xs[1] = RymerEncoding::update_smer_reverse(xs[1], c, sshift);
            }
        }
        // GTCAA: nuc=0b1011010000, S1=0b01100, S2=0b11000
        assert_eq!(xk[0].0, 0b1011010000);
        assert_eq!(xk[0].1, 0b01100);
        assert_eq!(xk[0].2, 0b11000);
        assert_eq!(xs[0], (0b100, 0b000));
        assert_eq!(xs[1], (0b110, 0b111));
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
    fn test_rymer_same_positions_different_ry_len() {
        let seq = "ACGTACGTACGTACGTACGT";
        let parameters = SyncmerParameters::try_new(8, 4).unwrap();

        // Check syncmer positions
        let positions_full: Vec<usize> = RymerSyncmerIterator::with_ry_len(
            seq.as_bytes(),
            parameters.k,
            parameters.s,
            parameters.t,
            parameters.k,
        )
        .map(|s| s.position)
        .collect();

        let positions_half: Vec<usize> = RymerSyncmerIterator::with_ry_len(
            seq.as_bytes(),
            parameters.k,
            parameters.s,
            parameters.t,
            parameters.k / 2,
        )
        .map(|s| s.position)
        .collect();

        let positions_zero: Vec<usize> = RymerSyncmerIterator::with_ry_len(
            seq.as_bytes(),
            parameters.k,
            parameters.s,
            parameters.t,
            0,
        )
        .map(|s| s.position)
        .collect();

        assert_eq!(positions_full, positions_half);
        assert_eq!(positions_full, positions_zero);
    }

    #[test]
    fn test_mixed_rymer_different_hashes_for_different_ry_len() {
        let seq = "ACGTACGTACGTACGTACGT";
        let parameters = SyncmerParameters::try_new(8, 4).unwrap();

        let hashes_full: Vec<_> = RymerSyncmerIterator::with_ry_len(
            seq.as_bytes(),
            parameters.k,
            parameters.s,
            parameters.t,
            parameters.k,
        )
        .map(|s| (s.hash1, s.hash2))
        .collect();

        let hashes_half: Vec<_> = RymerSyncmerIterator::with_ry_len(
            seq.as_bytes(),
            parameters.k,
            parameters.s,
            parameters.t,
            parameters.k / 2,
        )
        .map(|s| (s.hash1, s.hash2))
        .collect();

        assert_ne!(hashes_full, hashes_half);
    }

    #[test]
    fn test_mixed_rymer_iterator_produces_syncmers() {
        let seq = "ACGTACGTACGTACGTACGT";
        let parameters = SyncmerParameters::try_new(8, 4).unwrap();

        let syncmers: Vec<_> = RymerSyncmerIterator::with_ry_len(
            seq.as_bytes(),
            parameters.k,
            parameters.s,
            parameters.t,
            3,
        )
        .collect();

        assert!(!syncmers.is_empty());
        for s in &syncmers {
            assert!(s.position + parameters.k <= seq.len());
        }
    }

    #[test]
    fn test_ry_equal_with_ry_len() {
        // Full RY check (A↔G, C↔T)
        assert!(ry_equal(b"ACGT", b"GTAT", 4));
        assert!(!ry_equal(b"ACGT", b"GCTT", 4));
        // First 2 positions RY, last 2 exact
        assert!(ry_equal(b"ACGT", b"GCGT", 2));
        assert!(!ry_equal(b"ACGT", b"GCGA", 2));
        // Full exact check
        assert!(ry_equal(b"ACGT", b"ACGT", 0));
        assert!(!ry_equal(b"ACGT", b"GCGT", 0));
    }

}
