use std::cmp::min;
use std::collections::VecDeque;
use crate::hash::xxh64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Syncmer {
    pub hash: u64,
    pub position: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SyncmerParameters {
    pub k: usize,
    pub s: usize,
    pub t: usize,
}

impl SyncmerParameters {
    pub fn new(k: usize, s: usize) -> Self {
        SyncmerParameters { k, s, t: (k - s) / 2 + 1 }
    }

// TODO
// void verify() const {
//     if (k <= 7 || k > 32) {
//         throw BadParameter("k not in [8,32]");
//     }
//     if (s > k) {
//         throw BadParameter("s is larger than k");
//     }
//     if ((k - s) % 2 != 0) {
//         throw BadParameter("(k - s) must be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...");
//     }
// }
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
            kmask: (1 << (2*k)) - 1,
            smask: (1 << (2*s)) - 1,
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

impl<'a> Iterator for SyncmerIterator<'a> {
    type Item = Syncmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        for i in self.i..self.seq.len() {
            let ch = self.seq[i];
            let c = NUCLEOTIDES[ch as usize];
            if c < 4 { // not an "N" base
                self.xk[0] = ((self.xk[0] << 2) | (c as u64)) & self.kmask;        // forward strand
                self.xk[1] = (self.xk[1] >> 2) | (((3 - c) as u64) << self.kshift);  // reverse strand
                self.xs[0] = ((self.xs[0] << 2) | (c as u64)) & self.smask;        // forward strand
                self.xs[1] = (self.xs[1] >> 2) | (((3 - c) as u64) << self.sshift);  // reverse strand
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
                if self.qs.len() == (self.k - self.s + 1) { // We are at the last s-mer within the first k-mer, need to decide if we add it
                    // TODO use min
                    for j in 0..self.qs.len() {
                        if self.qs[j] <= self.qs_min_val {
                            self.qs_min_val = self.qs[j];
                        }
                    }
                } else {
                    // update queue and current minimum and position
                    let front = self.qs.pop_front().unwrap();
                    if front == self.qs_min_val {
                        // we popped a minimum, find new brute force
                        self.qs_min_val = u64::MAX;
                        for j in 0..self.qs.len() {
                            if self.qs[j] <= self.qs_min_val {
                                self.qs_min_val = self.qs[j];
                            }
                        }
                    } else if hash_s < self.qs_min_val { // the new value added to queue is the new minimum
                        self.qs_min_val = hash_s;
                    }
                }
                if self.qs[self.t - 1] == self.qs_min_val { // occurs at t:th position in k-mer
                    let yk = min(self.xk[0], self.xk[1]);
                    let syncmer = Syncmer { hash: syncmer_kmer_hash(yk), position: i + 1 - self.k };
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
    use super::*;

    #[test]
    fn test_syncmer_iterator() {
        let seq = "AAAAAAAAAAAAAAAAAAAA";
        assert_eq!(seq.len(), 20);

        let parameters = SyncmerParameters::new(8, 4);
        assert_eq!(parameters.t, 3);

        let mut iterator = SyncmerIterator::new(seq.as_bytes(), parameters.k, parameters.s, parameters.t);
        let syncmer = iterator.next().unwrap();
        assert_eq!(syncmer.position, 0);
    }
}
