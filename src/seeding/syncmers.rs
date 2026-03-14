use std::cmp::min;
use std::collections::VecDeque;

use super::InvalidSeedingParameter;
use super::hash::xxh64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Syncmer {
    /// Syncmer hash with canonicity stored in the least significant bit
    hash_and_canonical: u64,
    pub position: usize,
}

impl Syncmer {
    pub fn hash(&self) -> u64 {
        self.hash_and_canonical & !1
    }

    pub fn is_canonical(&self) -> bool {
        self.hash_and_canonical & 1 != 0
    }

    pub fn toggle_canonical(&mut self) {
        self.hash_and_canonical ^= 1;
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

pub struct SyncmerIterator<'a> {
    seq: &'a [u8],
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

impl<'a> SyncmerIterator<'a> {
    pub fn new(seq: &'a [u8], k: usize, s: usize, t: usize) -> SyncmerIterator<'a> {
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

impl Iterator for SyncmerIterator<'_> {
    type Item = Syncmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        for i in self.i..self.seq.len() {
            let ch = self.seq[i];
            let c = NUCLEOTIDES[ch as usize];
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
                        hash_and_canonical: (syncmer_kmer_hash(yk) & !1)
                            | ((self.xk[0] <= self.xk[1]) as u64),
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
    use std::{fs::File, io::BufReader};

    use crate::{io::fasta::read_fasta, revcomp::reverse_complement};

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
            SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let syncmer = iterator.next().unwrap();
        assert_eq!(syncmer.position, 0);
    }

    #[test]
    fn test_canonical_syncmers() {
        let mut f = BufReader::new(File::open("tests/phix.fasta").unwrap());
        let seq = read_fasta(&mut f).unwrap().pop().unwrap().sequence;
        let sequences = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".into(), seq];
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        for s in sequences {
            let revcomped = reverse_complement(&s);

            let syncmers_forward: Vec<_> =
                SyncmerIterator::new(&s, parameters.k, parameters.s, parameters.t).collect();

            let mut syncmers_reverse: Vec<_> =
                SyncmerIterator::new(&revcomped, parameters.k, parameters.s, parameters.t)
                    .collect();
            syncmers_reverse.reverse();
            for syncmer in &mut syncmers_reverse {
                syncmer.position = s.len() - parameters.k - syncmer.position;
            }

            assert_eq!(syncmers_forward, syncmers_reverse);
        }
    }
}
