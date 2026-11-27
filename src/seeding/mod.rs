mod hash;
pub mod parameters;
pub mod strobes;
pub mod syncmers;

pub use parameters::{InvalidSeedingParameter, Profile, SeedingParameters};
pub use strobes::{DEFAULT_AUX_LEN, RandstrobeIterator, RandstrobeParameters, RymerIterator};

use strobes::Randstrobe;
pub use syncmers::{
    KmerSyncmerIterator, RymerSyncmerIterator, Syncmer, SyncmerIterator, SyncmerParameters,
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
    if seq.len() < parameters.randstrobe.w_max {
        return [vec![], vec![]];
    }
    let sp = &parameters.syncmer;
    let expected = seq.len() / (sp.k - sp.s + 1);

    let mut syncmers = Vec::with_capacity(expected);
    if parameters.adna_mode {
        syncmers.extend(RymerSyncmerIterator::with_ry_len(
            seq,
            sp.k,
            sp.s,
            sp.t,
            parameters.ry_len,
        ));
    } else {
        syncmers.extend(KmerSyncmerIterator::new(seq, sp.k, sp.s, sp.t));
    }

    let push_randstrobes = |syncmers: &[Syncmer], out: &mut Vec<QueryRandstrobe>| {
        let syncmers = syncmers.iter().copied();
        let rp = parameters.randstrobe.clone();
        if parameters.adna_mode {
            push_query_randstrobes(
                RymerIterator::new(syncmers, rp, parameters.ry_len),
                sp.k,
                out,
            );
        } else {
            push_query_randstrobes(RandstrobeIterator::new(syncmers, rp), sp.k, out);
        }
    };

    let mut randstrobes = [
        Vec::with_capacity(syncmers.len()),
        Vec::with_capacity(syncmers.len()),
    ];
    push_randstrobes(&syncmers, &mut randstrobes[0]);

    // The forward syncmers are re-used for the reverse complement because
    // canonical/rymer syncmers are invariant under reverse complementing;
    // only the coordinates and orientation need to be adjusted.
    //
    // The randstrobes themselves cannot be re-used: if syncmer[i] and
    // syncmer[j] were paired up in the forward direction, syncmer[j] is not
    // necessarily paired with syncmer[i] in the reverse direction because i is
    // fixed in the forward direction and j is fixed in the reverse direction.
    syncmers.reverse();
    for syncmer in &mut syncmers {
        syncmer.position = seq.len() - syncmer.position - sp.k;
        syncmer.toggle_orientation();
    }
    push_randstrobes(&syncmers, &mut randstrobes[1]);

    randstrobes
}

/// Drain a randstrobe iterator into a randstrobe vector.
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
