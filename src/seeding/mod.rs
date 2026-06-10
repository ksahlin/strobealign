mod hash;
pub mod parameters;
pub mod strobes;
pub mod syncmers;

pub use parameters::{InvalidSeedingParameter, Profile, SeedingParameters};
pub use strobes::{DEFAULT_AUX_LEN, RandstrobeIterator, RandstrobeParameters, RymerIterator};

use strobes::Randstrobe;
pub use syncmers::{
    KmerSyncmerIterator, RymerSyncmer, RymerSyncmerIterator, Syncmer, SyncmerIterator,
    SyncmerParameters, ry_equal,
};

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

    // Forward and reverse-complement randstrobes. The forward syncmers are
    // re-used for the reverse complement, because canonical/rymer syncmers are
    // invariant under reverse complementing; only the coordinates and
    // orientation need to be adjusted.
    //
    // Note that the randstrobes themselves cannot be re-used: if syncmer[i] and
    // syncmer[j] were paired up in the forward direction, syncmer[j] is not
    // necessarily paired with syncmer[i] in the reverse direction because i is
    // fixed in the forward direction and j is fixed in the reverse direction.
    let sp = &parameters.syncmer;
    let rp = &parameters.randstrobe;
    if parameters.adna_mode {
        let mut syncmers: Vec<_> =
            RymerSyncmerIterator::with_ry_len(seq, sp.k, sp.s, sp.t, parameters.ry_len).collect();
        push_query_randstrobes(
            RymerIterator::new(syncmers.iter().copied(), rp.clone()),
            sp.k,
            &mut randstrobes[0],
        );
        syncmers.reverse();
        for syncmer in &mut syncmers {
            syncmer.position = seq.len() - syncmer.position - sp.k;
            syncmer.toggle_orientation();
        }
        push_query_randstrobes(
            RymerIterator::new(syncmers.into_iter(), rp.clone()),
            sp.k,
            &mut randstrobes[1],
        );
    } else {
        let mut syncmers: Vec<_> = KmerSyncmerIterator::new(seq, sp.k, sp.s, sp.t).collect();
        push_query_randstrobes(
            RandstrobeIterator::new(syncmers.iter().copied(), rp.clone()),
            sp.k,
            &mut randstrobes[0],
        );
        syncmers.reverse();
        for syncmer in &mut syncmers {
            syncmer.position = seq.len() - syncmer.position - sp.k;
            syncmer.toggle_orientation();
        }
        push_query_randstrobes(
            RandstrobeIterator::new(syncmers.into_iter(), rp.clone()),
            sp.k,
            &mut randstrobes[1],
        );
    }
    randstrobes
}

/// Drain a randstrobe iterator into a query randstrobe vector.
fn push_query_randstrobes(
    randstrobes: impl Iterator<Item = Randstrobe>,
    k: usize,
    out: &mut Vec<QueryRandstrobe>,
) {
    for randstrobe in randstrobes {
        out.push(QueryRandstrobe {
            hash: randstrobe.hash,
            hash_revcomp: randstrobe.hash_revcomp,
            start: randstrobe.strobe1_pos,
            end: randstrobe.strobe2_pos + k,
        });
    }
}
