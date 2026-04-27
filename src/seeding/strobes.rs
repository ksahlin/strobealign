use std::collections::VecDeque;

use super::{InvalidSeedingParameter, Syncmer};
use crate::index::{REF_RANDSTROBE_HASH_MASK, STROBE2_OFFSET_BITS};

pub const DEFAULT_AUX_LEN: u8 = 17;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RandstrobeParameters {
    pub w_min: usize,
    pub w_max: usize,
    pub q: u64,
    pub max_dist: u8,

    /// Mask for bits of the hash that represent the main hash
    pub main_hash_mask: u64,

    /// Mask for the main hash including its orientation bit
    pub forward_main_hash_mask: u64,

    /// Bit position of the partial orientation bit
    pub partial_orientation_pos: u32,
}

impl RandstrobeParameters {
    pub fn try_new(
        w_min: usize,
        w_max: usize,
        q: u64,
        max_dist: u8,
        main_hash_mask: u64,
    ) -> Result<Self, InvalidSeedingParameter> {
        let partial_orientation_pos = main_hash_mask.trailing_zeros() - 1;
        RandstrobeParameters {
            w_min,
            w_max,
            q,
            max_dist,
            main_hash_mask,
            forward_main_hash_mask: main_hash_mask | (1u64 << partial_orientation_pos),
            partial_orientation_pos,
        }
        .with_window(Some(w_min), Some(w_max))
    }

    pub fn with_window(
        mut self,
        w_min: Option<usize>,
        w_max: Option<usize>,
    ) -> Result<Self, InvalidSeedingParameter> {
        let new_w_min = w_min.unwrap_or(self.w_min);
        let new_w_max = w_max.unwrap_or(self.w_max);

        if w_min > w_max {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "w_min is greater than w_max (choose different -l/-u parameters)",
            ));
        }

        self.w_min = new_w_min;
        self.w_max = new_w_max;

        Ok(self)
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
    pub fn from_strobes(
        strobe1: Syncmer,
        strobe2: Syncmer,
        parameters: &RandstrobeParameters,
    ) -> Self {
        let is_forward1 = strobe1.is_forward() as u64;
        let is_forward2 = strobe2.is_forward() as u64;
        Randstrobe {
            hash: Randstrobe::hash(
                strobe1.hash(),
                strobe2.hash(),
                is_forward1,
                is_forward2,
                parameters,
            ),
            hash_revcomp: Randstrobe::hash(
                strobe2.hash(),
                strobe1.hash(),
                is_forward1 ^ 1,
                is_forward2 ^ 1,
                parameters,
            ),
            strobe1_pos: strobe1.position,
            strobe2_pos: strobe2.position,
        }
    }

    /// Combines two individual syncmer hashes into a randstrobe hash
    ///
    /// The first syncmer is designated as the "main", the other is
    /// the "auxiliary".
    /// The combined hash is obtained by setting the top bits to the bits of
    /// the main hash and the bottom bits to the bits of the auxiliary
    /// hash. Since entries in the index are sorted by randstrobe hash, this allows
    /// us to search for the main syncmer by masking out the lower bits.
    /// The orientation bit position (between main and aux) is left as zero.
    pub fn hash(
        hash1: u64,
        hash2: u64,
        is_forward1: u64,
        is_forward2: u64,
        parameters: &RandstrobeParameters,
    ) -> u64 {
        ((hash1 & parameters.main_hash_mask)
            | (hash2 & !parameters.forward_main_hash_mask)
            | (is_forward1 << parameters.partial_orientation_pos)
            | (is_forward2 << STROBE2_OFFSET_BITS))
            & (REF_RANDSTROBE_HASH_MASK << 1)
    }
}

pub struct RandstrobeIterator<I: Iterator<Item = Syncmer>> {
    parameters: RandstrobeParameters,
    syncmers: VecDeque<Syncmer>,
    syncmer_iterator: I,
}

impl<I: Iterator<Item = Syncmer>> RandstrobeIterator<I> {
    pub fn new(syncmer_iterator: I, parameters: RandstrobeParameters) -> RandstrobeIterator<I> {
        RandstrobeIterator {
            parameters,
            syncmers: VecDeque::<Syncmer>::new(),
            syncmer_iterator,
        }
    }
}

impl<SI: Iterator<Item = Syncmer>> Iterator for RandstrobeIterator<SI> {
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
        let max_position = strobe1.position + self.parameters.max_dist as usize;
        let mut min_val = u64::MAX;
        let mut strobe2 = self.syncmers[0]; // Defaults if no nearby syncmer

        for i in self.parameters.w_min..self.syncmers.len() {
            debug_assert!(i <= self.parameters.w_max);
            if self.syncmers[i].position > max_position {
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

        Some(Randstrobe::from_strobes(strobe1, strobe2, &self.parameters))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::io::fasta::{RefSequence, read_ref};
    use crate::seeding::SeedingParameters;
    use crate::seeding::SyncmerIterator;

    fn read_phix() -> RefSequence {
        read_ref("tests/phix.fasta")
            .unwrap()
            .first()
            .unwrap()
            .clone()
    }

    #[test]
    fn test_randstrobe_iterator() {
        let refseq = read_phix().sequence;
        let parameters = SeedingParameters::new(300);
        let syncmer_iter = SyncmerIterator::new(
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
        let parameters = SeedingParameters::new(100);
        let syncmer_iter = SyncmerIterator::new(
            &refseq,
            parameters.syncmer.k,
            parameters.syncmer.s,
            parameters.syncmer.t,
        );
        let syncmer_count = syncmer_iter.count();

        let syncmer_iter = SyncmerIterator::new(
            &refseq,
            parameters.syncmer.k,
            parameters.syncmer.s,
            parameters.syncmer.t,
        );
        let randstrobe_iter = RandstrobeIterator::new(syncmer_iter, parameters.randstrobe.clone());
        let randstrobe_count = randstrobe_iter.count();

        assert_eq!(randstrobe_count, syncmer_count);
    }
}
