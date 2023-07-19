use crate::index::{IndexParameters, StrobemerIndex};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;

pub fn map_single_end_read(seq: &Vec<u8>, index: &StrobemerIndex) {
    //Details details;
    //Timer strobe_timer;
    let query_randstrobes = randstrobes_query(&seq, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();
/*
    // Find NAMs
    // Timer nam_timer;
    auto [nonrepetitive_fraction, nams] = find_nams(query_randstrobes, index);
    // statistics.tot_find_nams += nam_timer.duration();
    /*
    TODO
    if (map_param.R > 1) {
        Timer rescue_timer;
        if (nams.empty() || nonrepetitive_fraction < 0.7) {
            details.nam_rescue = true;
            nams = find_nams_rescue(query_randstrobes, index, map_param.rescue_cutoff);
        }
        statistics.tot_time_rescue += rescue_timer.duration();
    }
    details.nams = nams.size();
    Timer nam_sort_timer;
    */
    std::sort(nams.begin(), nams.end(), score);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

    // Timer extend_timer;
    align_SE(
        aligner, sam, nams, record, index_parameters.syncmer.k,
        references, details, map_param.dropoff_threshold, map_param.maxTries,
        map_param.max_secondary
    );
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;

 */
    for qr in query_randstrobes {
        println!("qr: {:?}", qr);
    }
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    hash: u64,
    start: usize,
    end: usize,
    is_reverse: bool,
}

/// Generate randstrobes for a query sequence and its reverse complement.
pub fn randstrobes_query(seq: &Vec<u8>, parameters: &IndexParameters) -> Vec<QueryRandstrobe> {
    let mut randstrobes= Vec::<QueryRandstrobe>::new();
    if seq.len() < parameters.randstrobe.w_max {
        return randstrobes;
    }

    // TODO
    // For the reverse complement, we could re-use the syncmers of the forward
    // sequence because canonical syncmers are invariant under reverse
    // complementing. Only the coordinates need to be adjusted.

    let seq_rc = reverse_complement(seq);
    for (s, is_reverse) in [(seq, false), (&seq_rc, true)] {
        // Generate randstrobes for the forward sequence
        let mut syncmer_iter = SyncmerIterator::new(s, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
        let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &parameters.randstrobe);

        for randstrobe in randstrobe_iter {
            randstrobes.push(
                QueryRandstrobe {
                    hash: randstrobe.hash,
                    start: randstrobe.strobe1_pos,
                    end: randstrobe.strobe2_pos + parameters.syncmer.k,
                    is_reverse
                }
            );
        }
    }

    randstrobes
}
