mod hash;
pub mod parameters;
pub mod strobes;
pub mod syncmers;

pub use parameters::{InvalidSeedingParameter, SeedingParameters};
pub use strobes::{DEFAULT_AUX_LEN, RandstrobeIterator, RandstrobeParameters};
pub use syncmers::{Syncmer, SyncmerIterator, SyncmerParameters};

use crate::revcomp::reverse_complement;

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub hash_revcomp: u64,
    pub start: usize,
    pub end: usize,
}

/// Generate randstrobes for a query sequence and its reverse complement.
pub fn randstrobes_query(seq: &[u8], parameters: &SeedingParameters) -> [Vec<QueryRandstrobe>; 2] {
    let mut randstrobes = {
        let expected = seq.len() / (parameters.syncmer.k - parameters.syncmer.s + 1);
        [Vec::with_capacity(expected), Vec::with_capacity(expected)]
    };
    if seq.len() < parameters.randstrobe.w_max {
        return randstrobes;
    }

    let seq_rc = reverse_complement(seq);
    for (s, is_revcomp) in [(seq, false), (&seq_rc, true)] {
        // Generate randstrobes for the forward sequence
        let syncmer_iter = SyncmerIterator::new(
            s,
            parameters.syncmer.k,
            parameters.syncmer.s,
            parameters.syncmer.t,
        );
        let randstrobe_iter = RandstrobeIterator::new(syncmer_iter, parameters.randstrobe.clone());

        for randstrobe in randstrobe_iter {
            randstrobes[is_revcomp as usize].push(QueryRandstrobe {
                hash: randstrobe.hash,
                hash_revcomp: randstrobe.hash_revcomp,
                start: randstrobe.strobe1_pos,
                end: randstrobe.strobe2_pos + parameters.syncmer.k,
            });
        }
    }

    randstrobes
}
