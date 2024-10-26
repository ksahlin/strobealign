use std::cmp::{max, min};
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use fastrand::Rng;
use crate::fasta::RefSequence;
use crate::index::StrobemerIndex;
use crate::mapper;
use crate::mapper::QueryRandstrobe;
use crate::read::Read;

/// Non-overlapping approximate match
#[derive(Clone,Debug)]
pub struct Nam {
    pub nam_id: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    query_prev_hit_startpos: usize,
    ref_prev_hit_startpos: usize,
    pub n_hits: usize,
    pub ref_id: usize,
    pub score: u32,
    pub is_revcomp: bool,
}

impl Nam {
    pub fn ref_span(&self) -> usize {
        self.ref_end - self.ref_start
    }

    pub fn query_span(&self) -> usize {
        self.query_end - self.query_start
    }

    pub fn projected_ref_start(&self) -> usize {
        self.ref_start.saturating_sub(self.query_start)
    }
}

struct Hit {
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
}

/// Find a query’s NAMs, ignoring randstrobes that occur too often in the
/// reference (have a count above filter_cutoff).
///
/// Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
///
pub fn find_nams(query_randstrobes: &Vec<QueryRandstrobe>, index: &StrobemerIndex, filter_cutoff: usize) -> (f32, Vec<Nam>) {
    let mut hits_per_ref = [HashMap::with_capacity(100), HashMap::with_capacity(100)];
    let mut nr_good_hits = 0;
    let mut total_hits = 0;
    for randstrobe in query_randstrobes {
        if let Some(position) = index.get(randstrobe.hash) {
            total_hits += 1;
            if index.is_too_frequent(position, filter_cutoff) {
                continue;
            }
            nr_good_hits += 1;
            add_to_hits_per_ref(&mut hits_per_ref[randstrobe.is_revcomp as usize], randstrobe.start, randstrobe.end, index, position);
        }
    }
    let nonrepetitive_fraction = if total_hits > 0 { (nr_good_hits as f32) / (total_hits as f32) } else { 1.0 };
    let nams = merge_hits_into_nams_forward_and_reverse(&mut hits_per_ref, index.parameters.syncmer.k, false);

    (nonrepetitive_fraction, nams)
}

/// Find a query’s NAMs, using also some of the randstrobes that occur more often
/// than the normal (non-rescue) filter cutoff.
pub fn find_nams_rescue(
    query_randstrobes: &[QueryRandstrobe],
    index: &StrobemerIndex,
    rescue_cutoff: usize,
) -> Vec<Nam> {

    struct RescueHit {
        count: usize,
        position: usize,
        query_start: usize,
        query_end: usize,
    }
/*
        bool operator< (const RescueHit& rhs) const {
            return std::tie(count, query_s, query_e, is_revcomp)
                < std::tie(rhs.count, rhs.query_s, rhs.query_e, rhs.is_revcomp);
        }
    }*/

    let mut hits_per_ref = [HashMap::with_capacity(100), HashMap::with_capacity(100)];
    let mut hits_fw = Vec::with_capacity(5000);
    let mut hits_rc = Vec::with_capacity(5000);

    for randstrobe in query_randstrobes {

        if let Some(position) = index.get(randstrobe.hash) {
            let count = index.get_count(position);
            let rh = RescueHit {
                count,
                position,
                query_start: randstrobe.start,
                query_end: randstrobe.end,
            };
            if randstrobe.is_revcomp {
                hits_rc.push(rh);
            } else {
                hits_fw.push(rh);
            }
        }
    }

    let cmp = |a: &RescueHit, b: &RescueHit| (a.count, a.query_start, a.query_end).cmp(&(b.count, b.query_start, b.query_end));
    hits_fw.sort_by(cmp);
    hits_rc.sort_by(cmp);

    for (is_revcomp, rescue_hits) in [(false, hits_fw), (true, hits_rc)] {
        for (i, rh) in rescue_hits.iter().enumerate() {
            if (rh.count > rescue_cutoff && i >= 5) || rh.count > 1000 {
                break;
            }
            add_to_hits_per_ref(&mut hits_per_ref[is_revcomp as usize], rh.query_start, rh.query_end, index, rh.position);
        }
    }

    merge_hits_into_nams_forward_and_reverse(&mut hits_per_ref, index.parameters.syncmer.k, true)
}

fn add_to_hits_per_ref(
    hits_per_ref: &mut HashMap<usize, Vec<Hit>>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    position: usize,
) {
    let query_length = query_end - query_start;
    let mut max_length_diff = usize::MAX;
    for randstrobe in &index.randstrobes[position..] {
        if randstrobe.hash() != index.randstrobes[position].hash() {
            break;
        }
        let ref_start = randstrobe.position();
        let ref_end = ref_start + randstrobe.strobe2_offset() + index.parameters.syncmer.k;
        let ref_length = ref_end - ref_start;
        let length_diff = (query_length as isize - ref_length as isize).unsigned_abs();
        if length_diff <= max_length_diff {
            let hit = Hit{query_start, query_end, ref_start, ref_end};
            let ref_id = randstrobe.reference_index();

            if let Entry::Vacant(e) = hits_per_ref.entry(ref_id) {
                e.insert(vec![hit]);
            } else {
                hits_per_ref.get_mut(&ref_id).unwrap().push(hit);
            }
            max_length_diff = length_diff;
        }
    }
}

// TODO should not be mut
fn merge_hits_into_nams(hits_per_ref: &mut HashMap<usize, Vec<Hit>>, k: usize, sort: bool, is_revcomp: bool, nams: &mut Vec<Nam>) {
    for (ref_id, hits) in hits_per_ref.iter_mut() {
        if sort {
            hits.sort_by_key(|k| (k.query_start, k.ref_start));
        }

        let mut open_nams: Vec<Nam> = Vec::new();
        let mut prev_q_start = 0;
        for h in hits {
            let mut is_added = false;
            for o in &mut open_nams {
                // Extend NAM
                if (o.query_prev_hit_startpos < h.query_start)
                    && (h.query_start <= o.query_end)
                    && (o.ref_prev_hit_startpos < h.ref_start)
                    && (h.ref_start <= o.ref_end)
                {
                    if (h.query_end > o.query_end) && (h.ref_end > o.ref_end) {
                        o.query_end = h.query_end;
                        o.ref_end = h.ref_end;
                        o.query_prev_hit_startpos = h.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_hits += 1;
                        is_added = true;
                        break;
                    } else if (h.query_end <= o.query_end) && (h.ref_end <= o.ref_end) {
                        o.query_prev_hit_startpos = h.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_hits += 1;
                        is_added = true;
                        break;
                    }
                }
            }
            // Add the hit to open matches
            if !is_added {
                open_nams.push(Nam {
                    nam_id: nams.len() + open_nams.len(),
                    query_start: h.query_start,
                    query_end: h.query_end,
                    ref_start: h.ref_start,
                    ref_end: h.ref_end,
                    ref_id: *ref_id,
                    query_prev_hit_startpos: h.query_start,
                    ref_prev_hit_startpos: h.ref_start,
                    n_hits: 1,
                    is_revcomp,
                    score: 0,
                });
            }

            // Only filter if we have advanced at least k nucleotides
            if h.query_start > prev_q_start + k {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for n in &open_nams {
                    if n.query_end < h.query_start {
                        let n_max_span = max(n.query_span(), n.ref_span());
                        let n_min_span = min(n.query_span(), n.ref_span());
                        let n_score =
                            if 2 * n_min_span - n_max_span > 0 {
                                // this is really just n_hits * (min_span - (offset_in_span) ) );
                                n.n_hits * (2 * n_min_span - n_max_span)
                            } else {
                                1
                            } as u32;
//                        n_score = n.n_hits * n.query_span();
                        let mut nam = n.clone();
                        nam.score = n_score;
                        nams.push(nam);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                let c = h.query_start;

                open_nams.retain(|x| x.query_end >= c);

                prev_q_start = h.query_start;
            }
        }

        // Add all current open_matches to final NAMs
        for mut n in open_nams {
            let n_max_span = max(n.query_span(), n.ref_span());
            let n_min_span = min(n.query_span(), n.ref_span());
            n.score =
                if 2 * n_min_span > n_max_span {
                    n.n_hits * (2 * n_min_span - n_max_span)
                } else {
                    1
                } as u32;
            nams.push(n);
        }
    }
}

fn merge_hits_into_nams_forward_and_reverse(hits_per_ref: &mut [HashMap<usize, Vec<Hit>>; 2], k: usize, sort: bool) -> Vec<Nam> {
    let mut nams = Vec::new();
    for is_revcomp in [false, true] {
        merge_hits_into_nams(&mut hits_per_ref[is_revcomp as usize], k, sort, is_revcomp, &mut nams);
    }

    nams
}

/// Determine whether the NAM represents a match to the forward or
/// reverse-complemented sequence by checking in which orientation the
/// first and last strobe in the NAM match
///
/// - If first and last strobe match in forward orientation, return true.
/// - If first and last strobe match in reverse orientation, update the NAM
///   in place and return true.
/// - If first and last strobe do not match consistently, return false.
pub fn reverse_nam_if_needed(nam: &mut Nam, read: &Read, references: &[RefSequence], k: usize) -> bool {
    let ref_start_kmer = &references[nam.ref_id].sequence[nam.ref_start..nam.ref_start + k];
    let ref_end_kmer = &references[nam.ref_id].sequence[nam.ref_end - k..nam.ref_end];

    let (seq, seq_rc) = if nam.is_revcomp {
        (read.rc(), read.seq())
    } else {
        (read.seq(), read.rc())
    };
    let read_start_kmer = &seq[nam.query_start..nam.query_start + k];
    let read_end_kmer = &seq[nam.query_end - k.. nam.query_end];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    // we need two extra checks for this - hopefully this will remove all the false hits we see
    // (true hash collisions should be very few)
    let read_len = read.len();
    let q_start_tmp = read_len - nam.query_end;
    let q_end_tmp = read_len - nam.query_start;
    // false reverse hit, change coordinates in nam to forward
    let read_start_kmer = &seq_rc[q_start_tmp..q_start_tmp + k];
    let read_end_kmer = &seq_rc[q_end_tmp - k..q_end_tmp];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        nam.is_revcomp = !nam.is_revcomp;
        nam.query_start = q_start_tmp;
        nam.query_end = q_end_tmp;
        true
    } else {
        false
    }
}

// TODO rename
pub fn get_nams(sequence: &[u8], index: &StrobemerIndex, rescue_level: usize, rng: &mut Rng) -> (bool, Vec<Nam>) {
    //Timer strobe_timer;
    let query_randstrobes = mapper::randstrobes_query(sequence, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();

    // Timer nam_timer;
    let (nonrepetitive_fraction, mut nams) = find_nams(&query_randstrobes, index, index.filter_cutoff);
    // statistics.tot_find_nams += nam_timer.duration();

    let mut nam_rescue = false;
    if rescue_level > 1 {
        // Timer rescue_timer;
        if nams.is_empty() || nonrepetitive_fraction < 0.7 {
            nam_rescue = true;
            nams = find_nams_rescue(&query_randstrobes, index, index.rescue_cutoff);
        }
        // statistics.tot_time_rescue += rescue_timer.duration();
    }
    // Timer nam_sort_timer;

    nams.sort_by_key(|k| -(k.score as i32));
    shuffle_top_nams(&mut nams, rng);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

    (nam_rescue, nams)
}

/// Shuffle the top-scoring NAMs. Input must be sorted by score.
/// This helps to ensure we pick a random location in case there are multiple
/// equally good ones.
fn shuffle_top_nams(nams: &mut [Nam], rng: &mut Rng) {
    if let Some(best) = nams.first() {
        let best_score = best.score;

        let pos = nams.iter().position(|nam| nam.score != best_score);
        let end = pos.unwrap_or(nams.len());
        if end > 1 {
            rng.shuffle(&mut nams[0..end]);
        }
    }
}
