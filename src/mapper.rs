use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::find_nams;
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;

pub fn map_single_end_read(seq: &Vec<u8>, index: &StrobemerIndex) {
    //Details details;
    //Timer strobe_timer;
    let query_randstrobes = randstrobes_query(&seq, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();

    // Timer nam_timer;
    let (nonrepetitive_fraction, mut nams) = find_nams(&query_randstrobes, index);
    println!("nonrepetitive fraction: {}. nams.len(): {}", nonrepetitive_fraction, nams.len());
    for nam in &nams {
        println!("{:?}", nam);
    }

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
    nams.sort_by_key(|&k| k.score);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

    // Timer extend_timer;
    // align_SE(
    //     aligner, sam, nams, record, index.parameters.syncmer.k,
    //     references, details, map_param.dropoff_threshold, map_param.maxTries,
    //     map_param.max_secondary
    // );

    /*
    if nams.is_empty() {
        sam.add_unmapped(record);
        return;
    }

    Read read(record.seq);
    let alignments = std::vector<Alignment> alignments;
    int tries = 0;
    float score_dropoff;
    Nam n_max = all_nams[0];

    int best_align_dist = ~0U >> 1;
    int best_align_sw_score = -1000;

    Alignment best_sam_aln;
    best_sam_aln.sw_score = -100000;
    best_sam_aln.is_unaligned = true;
    int min_mapq_diff = best_align_dist;
    for (auto &n : all_nams) {
        score_dropoff = (float) n.n_hits / n_max.n_hits;

        if (tries >= max_tries || best_align_dist == 0 || score_dropoff < dropoff) { // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
            break;
        }
        bool consistent_nam = reverse_nam_if_needed(n, read, references, k);
        details.nam_inconsistent += !consistent_nam;
        auto sam_aln = get_alignment(aligner, n, references, read, consistent_nam);
        details.tried_alignment++;
        details.gapped += sam_aln.gapped;

        int diff_to_best = std::abs(best_align_sw_score - sam_aln.sw_score);
        min_mapq_diff = std::min(min_mapq_diff, diff_to_best);

        if (max_secondary > 0) {
            alignments.emplace_back(sam_aln);
        }
        if (sam_aln.sw_score > best_align_sw_score) {
            min_mapq_diff = std::max(0, sam_aln.sw_score - best_align_sw_score); // new distance to next best match
            best_align_sw_score = sam_aln.sw_score;
            best_sam_aln = std::move(sam_aln);
            if (max_secondary == 0) {
                best_align_dist = best_sam_aln.global_ed;
            }
        }
        tries++;
    }

    if (max_secondary == 0) {
        best_sam_aln.mapq = std::min(min_mapq_diff, 60);
        sam.add(best_sam_aln, record, read.rc, true, details);
        return;
    }
    // Sort alignments by score, highest first
    std::sort(alignments.begin(), alignments.end(),
        [](const Alignment& a, const Alignment& b) -> bool {
            return a.sw_score > b.sw_score;
        }
    );

    auto max_out = std::min(alignments.size(), static_cast<size_t>(max_secondary) + 1);
    bool is_primary{true};
    for (size_t i = 0; i < max_out; ++i) {
        auto sam_aln = alignments[i];
        if ((sam_aln.sw_score - best_align_sw_score) > (2*aligner.parameters.mismatch + aligner.parameters.gap_open) ){
            break;
        }
        if (!is_primary) {
            sam_aln.mapq = 255;
        } else {
            sam_aln.mapq = std::min(min_mapq_diff, 60);
        }
        sam.add(sam_aln, record, read.rc, is_primary, details);
        is_primary = false;
    }



    */


    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub start: usize,
    pub end: usize,
    pub is_reverse: bool,
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
