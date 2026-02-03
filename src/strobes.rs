use std::collections::VecDeque;

use crate::index::{InvalidIndexParameter, REF_RANDSTROBE_HASH_MASK};
use crate::syncmers::{RymerSyncmer, SyncmerLike};

pub const DEFAULT_AUX_LEN: u8 = 17;

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RandstrobeHashMethod {
    /// Combines two individual syncmer hashes into a randstrobe hash
    ///
    /// The first syncmer is designated as the "main", the other is
    /// the "auxiliary".
    /// The combined hash is obtained by setting the top bits to the bits of
    /// the main hash and the bottom bits to the bits of the auxiliary
    /// hash. Since entries in the index are sorted by randstrobe hash, this allows
    /// us to search for the main syncmer by masking out the lower bits.
    #[default]
    McsHash = 0,
}

impl RandstrobeHashMethod {
    #[inline]
    pub fn hash(&self, hash1: u64, hash2: u64, main_hash_mask: u64) -> u64 {
        match self {
            Self::McsHash => {
                ((hash1 & main_hash_mask) | (hash2 & !main_hash_mask)) & REF_RANDSTROBE_HASH_MASK
            }
        }
    }

    #[inline]
    pub fn hash_revcomp(&self, hash1: u64, hash2: u64, main_hash_mask: u64) -> u64 {
        self.hash(hash2, hash1, main_hash_mask)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RandstrobeParameters {
    pub w_min: usize,
    pub w_max: usize,
    pub q: u64,
    pub max_dist: u8,
    pub hash_method: RandstrobeHashMethod,

    /// Mask for bits of the hash that represent the main hash
    pub main_hash_mask: u64,
}

impl RandstrobeParameters {
    pub fn try_new(
        w_min: usize,
        w_max: usize,
        q: u64,
        max_dist: u8,
        hash_method: RandstrobeHashMethod,
        main_hash_mask: u64,
    ) -> Result<Self, InvalidIndexParameter> {
        if w_min > w_max {
            return Err(InvalidIndexParameter::InvalidParameter(
                "w_min is greater than w_max (choose different -l/-u parameters)",
            ));
        }
        Ok(RandstrobeParameters {
            w_min,
            w_max,
            q,
            max_dist,
            hash_method,
            main_hash_mask,
        })
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Randstrobe {
    pub hash: u64,
    pub hash_revcomp: u64,
    pub strobe1_pos: usize,
    pub strobe2_pos: usize,
}

impl Randstrobe {
    pub fn from_strobes<S: SyncmerLike>(
        hash_method: RandstrobeHashMethod,
        strobe1: S,
        strobe2: S,
        main_hash_mask: u64,
    ) -> Self {
        Randstrobe {
            hash: hash_method.hash(strobe1.hash(), strobe2.hash(), main_hash_mask),
            hash_revcomp: hash_method.hash_revcomp(strobe1.hash(), strobe2.hash(), main_hash_mask),
            strobe1_pos: strobe1.position(),
            strobe2_pos: strobe2.position(),
        }
    }
}

pub struct RandstrobeIterator<S: SyncmerLike, I: Iterator<Item = S>> {
    parameters: RandstrobeParameters,
    syncmers: VecDeque<S>,
    syncmer_iterator: I,
}

impl<S: SyncmerLike, I: Iterator<Item = S>> RandstrobeIterator<S, I> {
    pub fn new(syncmer_iterator: I, parameters: RandstrobeParameters) -> RandstrobeIterator<S, I> {
        RandstrobeIterator {
            parameters,
            syncmers: VecDeque::<S>::new(),
            syncmer_iterator,
        }
    }
}

impl<S: SyncmerLike, SI: Iterator<Item = S>> Iterator for RandstrobeIterator<S, SI> {
    type Item = Randstrobe;
    fn next(&mut self) -> Option<Self::Item> {
        while self.syncmers.len() <= self.parameters.w_max {
            if let Some(syncmer) = self.syncmer_iterator.next() {
                self.syncmers.push_back(syncmer);
            } else {
                break;
            }
        }
        if self.syncmers.is_empty() {
            return None;
        }
        let strobe1 = self.syncmers[0];
        let max_position = strobe1.position() + self.parameters.max_dist as usize;
        let mut min_val = u64::MAX;
        let mut strobe2 = self.syncmers[0]; // Defaults if no nearby syncmer

        for i in self.parameters.w_min..self.syncmers.len() {
            debug_assert!(i <= self.parameters.w_max);
            if self.syncmers[i].position() > max_position {
                break;
            }
            let b = (strobe1.hash() ^ self.syncmers[i].hash()) & self.parameters.q;
            let ones = b.count_ones() as u64;
            if ones < min_val {
                min_val = ones;
                strobe2 = self.syncmers[i];
            }
        }
        self.syncmers.pop_front();

        Some(Randstrobe::from_strobes(
            self.parameters.hash_method,
            strobe1,
            strobe2,
            self.parameters.main_hash_mask,
        ))
    }
}

pub struct RymerIterator<I: Iterator<Item = RymerSyncmer>> {
    syncmer_iterator: I,
    main_hash_mask: u64,
    hash_method: RandstrobeHashMethod,
}

impl<I: Iterator<Item = RymerSyncmer>> RymerIterator<I> {
    pub fn new(syncmer_iterator: I, main_hash_mask: u64) -> Self {
        RymerIterator {
            syncmer_iterator,
            main_hash_mask,
            hash_method: RandstrobeHashMethod::McsHash,
        }
    }
}

impl<I: Iterator<Item = RymerSyncmer>> Iterator for RymerIterator<I> {
    type Item = Randstrobe;

    fn next(&mut self) -> Option<Self::Item> {
        let syncmer = self.syncmer_iterator.next()?;
        let hash1 = syncmer.hash1 as u64;
        let hash2 = syncmer.hash2 as u64;

        Some(Randstrobe {
            hash: self.hash_method.hash(hash1, hash2, self.main_hash_mask),
            hash_revcomp: self
                .hash_method
                .hash_revcomp(hash1, hash2, self.main_hash_mask),
            strobe1_pos: syncmer.position,
            strobe2_pos: syncmer.position,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fasta::{RefSequence, read_fasta};
    use crate::index::IndexParameters;
    use crate::syncmers::{KmerSyncmerIterator, RymerSyncmerIterator, SyncmerParameters};
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn test_randstrobe_hash() {
        let main_hash_mask = 0xfffffffffc000000;
        let hash1 = 0x00000000ffffffff;
        let hash2 = 0xffffffff00000000;
        let hash_method = RandstrobeHashMethod::McsHash;

        assert!(hash_method.hash(hash1, hash2, main_hash_mask) == 0x00000000fc000000);
        assert!(hash_method.hash_revcomp(hash1, hash2, main_hash_mask) == 0xffffffff03ffff00);
    }

    fn read_phix() -> RefSequence {
        let f = File::open("tests/phix.fasta").unwrap();
        let mut reader = BufReader::new(f);

        read_fasta(&mut reader).unwrap().first().unwrap().clone()
    }

    #[test]
    fn test_randstrobe_iterator() {
        let refseq = read_phix().sequence;
        let parameters = IndexParameters::default_from_read_length(300);
        let syncmer_iter = KmerSyncmerIterator::new(
            &refseq,
            parameters.syncmer.k,
            parameters.syncmer.s,
            parameters.syncmer.t,
        );
        let randstrobe_iter = RandstrobeIterator::new(syncmer_iter, parameters.randstrobe.clone());

        for randstrobe in randstrobe_iter {
            assert!(randstrobe.hash > 0);
            assert!(randstrobe.strobe2_pos >= randstrobe.strobe2_pos);
            assert!(randstrobe.strobe1_pos < refseq.len());
            assert!(randstrobe.strobe2_pos < refseq.len());
        }
    }

    // Ensure SyncmerIterator and RandstrobeIterator return the same number of
    // items. We need this to hold for `count_randstrobes()`.
    #[test]
    fn test_syncmer_and_randstrobe_iterator_same_count() {
        let refseq = read_phix().sequence;
        let parameters = IndexParameters::default_from_read_length(100);
        let syncmer_iter = KmerSyncmerIterator::new(
            &refseq,
            parameters.syncmer.k,
            parameters.syncmer.s,
            parameters.syncmer.t,
        );
        let syncmer_count = syncmer_iter.count();

        let syncmer_iter = KmerSyncmerIterator::new(
            &refseq,
            parameters.syncmer.k,
            parameters.syncmer.s,
            parameters.syncmer.t,
        );
        let randstrobe_iter = RandstrobeIterator::new(syncmer_iter, parameters.randstrobe.clone());
        let randstrobe_count = randstrobe_iter.count();

        assert_eq!(randstrobe_count, syncmer_count);
    }

    #[test]
    fn test_rymer_iterator() {
        let refseq = read_phix().sequence;
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        let main_hash_mask = 0xfffffffffc000000u64;

        let syncmer_iter =
            RymerSyncmerIterator::new(&refseq, parameters.k, parameters.s, parameters.t);
        let rymer_iter = RymerIterator::new(syncmer_iter, main_hash_mask);

        for randstrobe in rymer_iter {
            assert!(randstrobe.strobe1_pos == randstrobe.strobe2_pos);
            assert!(randstrobe.strobe1_pos < refseq.len());
        }
    }

    #[test]
    fn test_rymer_iterator_count_matches_syncmer_count() {
        let refseq = read_phix().sequence;
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        let main_hash_mask = 0xfffffffffc000000u64;

        let syncmer_iter =
            RymerSyncmerIterator::new(&refseq, parameters.k, parameters.s, parameters.t);
        let syncmer_count = syncmer_iter.count();

        let syncmer_iter =
            RymerSyncmerIterator::new(&refseq, parameters.k, parameters.s, parameters.t);
        let rymer_iter = RymerIterator::new(syncmer_iter, main_hash_mask);
        let rymer_count = rymer_iter.count();

        assert_eq!(rymer_count, syncmer_count);
    }
}
