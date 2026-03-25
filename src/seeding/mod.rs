mod hash;
pub mod parameters;
pub mod strobes;
pub mod syncmers;

pub use parameters::{InvalidSeedingParameter, SeedingParameters};
pub use strobes::{DEFAULT_AUX_LEN, RandstrobeIterator, RandstrobeParameters, RymerIterator};
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

    if parameters.adna_mode {
        randstrobes_query_adna(seq, parameters, &mut randstrobes);
    } else {
        randstrobes_query_kmer(seq, parameters, &mut randstrobes);
    }
    randstrobes
}

fn randstrobes_query_kmer(
    seq: &[u8],
    parameters: &SeedingParameters,
    randstrobes: &mut [Vec<QueryRandstrobe>; 2],
) {
    // Generate syncmers for the forward sequence
    let syncmer_iter = KmerSyncmerIterator::new(
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
        syncmers[i].toggle_canonical();
    }

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
}

fn randstrobes_query_adna(
    seq: &[u8],
    parameters: &SeedingParameters,
    randstrobes: &mut [Vec<QueryRandstrobe>; 2],
) {
    // Generate rymer syncmers for the forward sequence
    let syncmer_iter = RymerSyncmerIterator::with_ry_len(
        seq,
        parameters.syncmer.k,
        parameters.syncmer.s,
        parameters.syncmer.t,
        parameters.ry_len,
    );
    let mut syncmers: Vec<_> = syncmer_iter.collect();

    // Generate randstrobes for the forward sequence
    let randstrobe_iter =
        RymerIterator::new(syncmers.iter().cloned(), parameters.randstrobe.clone());

    for randstrobe in randstrobe_iter {
        randstrobes[0].push(QueryRandstrobe {
            hash: randstrobe.hash,
            hash_revcomp: randstrobe.hash_revcomp,
            start: randstrobe.strobe1_pos,
            end: randstrobe.strobe2_pos + parameters.syncmer.k,
        });
    }

    // For the reverse complement, re-use syncmers with adjusted coordinates
    syncmers.reverse();
    for i in 0..syncmers.len() {
        syncmers[i].position = seq.len() - syncmers[i].position - parameters.syncmer.k;
        syncmers[i].toggle_canonical();
    }

    let rc_randstrobe_iter =
        RymerIterator::new(syncmers.into_iter(), parameters.randstrobe.clone());
    for randstrobe in rc_randstrobe_iter {
        randstrobes[1].push(QueryRandstrobe {
            hash: randstrobe.hash,
            hash_revcomp: randstrobe.hash_revcomp,
            start: randstrobe.strobe1_pos,
            end: randstrobe.strobe2_pos + parameters.syncmer.k,
        });
    }
}
