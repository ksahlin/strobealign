mod hash;
pub mod parameters;
pub mod strobes;
pub mod syncmers;

pub use parameters::{InvalidSeedingParameter, SeedingParameters};
pub use strobes::{DEFAULT_AUX_LEN, RandstrobeIterator, RandstrobeParameters};
pub use syncmers::{Syncmer, SyncmerIterator, SyncmerParameters};

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

    // Generate syncmers for the forward sequence
    let syncmer_iter = SyncmerIterator::new(
        seq,
        parameters.syncmer.k,
        parameters.syncmer.s,
        parameters.syncmer.t,
    );
    let mut syncmers: Vec<_> = syncmer_iter.collect();

    // Generate randstrobes for the forward sequence
    let randstrobe_iter =
        RandstrobeIterator::new(syncmers.iter().cloned(), parameters.randstrobe.clone());

    for randstrobe in randstrobe_iter {
        randstrobes[0].push(QueryRandstrobe {
            hash: randstrobe.hash,
            hash_revcomp: randstrobe.hash_revcomp,
            start: randstrobe.strobe1_pos,
            end: randstrobe.strobe2_pos + parameters.syncmer.k,
        });
    }

    // For the reverse complement, we can re-use the syncmers of the forward
    // sequence because canonical syncmers are invariant under reverse
    // complementing. Only the coordinates need to be adjusted.
    syncmers.reverse();
    for i in 0..syncmers.len() {
        syncmers[i].position = seq.len() - syncmers[i].position - parameters.syncmer.k;
    }

    // Randstrobes cannot be re-used for the reverse complement:
    // If in the forward direction, syncmer[i] and syncmer[j] were paired up, it
    // is not necessarily the case that syncmer[j] is going to be paired with
    // syncmer[i] in the reverse direction because i is fixed in the forward
    // direction and j is fixed in the reverse direction.
    let rc_randstrobe_iter =
        RandstrobeIterator::new(syncmers.into_iter(), parameters.randstrobe.clone());
    for randstrobe in rc_randstrobe_iter {
        randstrobes[1].push(QueryRandstrobe {
            hash: randstrobe.hash,
            hash_revcomp: randstrobe.hash_revcomp,
            start: randstrobe.strobe1_pos,
            end: randstrobe.strobe2_pos + parameters.syncmer.k,
        });
    }
    randstrobes
}
