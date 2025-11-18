use std::cmp::{max, min};
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::hash::{BuildHasherDefault, DefaultHasher};
use std::time::Instant;
use fastrand::Rng;
use log::Level::Trace;
use log::trace;
use crate::chainer::Chainer;
use crate::details::NamDetails;
use crate::fasta::RefSequence;
use crate::hit::Hit;
use crate::index::StrobemerIndex;
use crate::mapper;
use crate::mapper::QueryRandstrobe;
use crate::mcsstrategy::McsStrategy;
use crate::read::Read;

/// Non-overlapping approximate match
#[derive(Clone, Debug)]
pub struct Nam {
    pub nam_id: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub query_prev_match_startpos: usize,
    pub ref_prev_match_startpos: usize,
    pub n_matches: usize,
    pub ref_id: usize,
    pub score: f32,
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

impl Display for Nam {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Nam(ref_id={}, query: {}..{}, ref: {}..{}, rc={}, score={})",
               self.ref_id, self.query_start, self.query_end, self.ref_start, self.ref_end, self.is_revcomp as u8, self.score
        )?;
        Ok(())
    }
}

#[derive(Debug, PartialEq, Eq)]
struct Match {
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
}

// The regular HashMap uses a randomly seeded hash function, but we want reproducible results
type DefaultHashMap<K, V> = HashMap<K, V, BuildHasherDefault<DefaultHasher>>;

/// Find a query’s NAMs, using also some of the randstrobes that occur more often
/// than the normal (non-rescue) filter cutoff.
///
/// Return the number of hits and the matches_map
fn find_matches_rescue(
    query_randstrobes: &[QueryRandstrobe],
    index: &StrobemerIndex,
    rescue_cutoff: usize,
    _mcs_strategy: McsStrategy,
) -> (usize, DefaultHashMap<usize, Vec<Match>>) {

    // TODO we ignore mcs_strategy
    // if we implement mcs_strategy, we also need to return the no. of partial_hits

    struct RescueHit {
        count: usize,
        position: usize,
        query_start: usize,
        query_end: usize,
    }

    let mut matches_map = DefaultHashMap::default();
    matches_map.reserve(100);
    let mut n_hits = 0;

    let mut hits = Vec::with_capacity(5000);

    for randstrobe in query_randstrobes {
        if let Some(position) = index.get_full(randstrobe.hash) {
            let count = index.get_count_full(position);
            let rh = RescueHit {
                count,
                position,
                query_start: randstrobe.start,
                query_end: randstrobe.end,
            };
            hits.push(rh);
        }
    }

    let cmp = |a: &RescueHit, b: &RescueHit| (a.count, a.query_start, a.query_end).cmp(&(b.count, b.query_start, b.query_end));
    hits.sort_by(cmp);

    let rescue_hits = &hits;
    for (i, rh) in rescue_hits.iter().enumerate() {
        if (rh.count > rescue_cutoff && i >= 5) || rh.count > 1000 {
            break;
        }
        add_to_matches_map_full(&mut matches_map, rh.query_start, rh.query_end, index, rh.position);
        n_hits += 1;
    }

    (n_hits, matches_map)
}

fn add_to_matches_map_full(
    matches_map: &mut DefaultHashMap<usize, Vec<Match>>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    position: usize,
) {
    let query_length = query_end - query_start;
    let mut min_length_diff = usize::MAX;
    for randstrobe in &index.randstrobes[position..] {
        if randstrobe.hash() != index.randstrobes[position].hash() {
            break;
        }
        let ref_start = randstrobe.position();
        let ref_end = ref_start + randstrobe.strobe2_offset() + index.k();
        let ref_length = ref_end - ref_start;
        let length_diff = (query_length as isize - ref_length as isize).unsigned_abs();
        if length_diff <= min_length_diff {
            let match_ = Match {query_start, query_end, ref_start, ref_end};
            let ref_id = randstrobe.reference_index();
            matches_map.entry(ref_id).or_default().push(match_);
            min_length_diff = length_diff;
        }
    }
}

fn add_to_matches_map_partial(
    matches_map: &mut DefaultHashMap<usize, Vec<Match>>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    position: usize,
) {
    let hash = index.get_hash_partial(position);
    for pos in position..index.randstrobes.len() {
        if index.get_hash_partial(pos) != hash {
            break;
        }
        let randstrobe = &index.randstrobes[pos];
        let ref_id = randstrobe.reference_index();
        let (ref_start, ref_end) = index.strobe_extent_partial(pos);
        let match_ = Match { query_start, query_end, ref_start, ref_end };
        matches_map.entry(ref_id).or_default().push(match_);
    }
}

// TODO should not be mut
fn merge_matches_into_nams(
    matches_map: &mut DefaultHashMap<usize, Vec<Match>>, k: usize, sort: bool, is_revcomp: bool, nams: &mut Vec<Nam>
) {
    for (ref_id, matches) in matches_map.iter_mut() {
        if sort {
            // -k.query_end to prefer full matches over partial ones
            matches.sort_by_key(|k| (k.query_start, k.ref_start, -(k.query_end as isize)));
        }

        let mut open_nams: Vec<Nam> = Vec::new();
        let mut prev_q_start = 0;
        let prev_match = Match {query_start: 0, query_end: 0, ref_start: 0, ref_end: 0};
        for m in matches {
            assert_ne!(prev_match, *m);

            if m.query_end - m.query_start == k && (prev_match.query_start == m.query_start) && (prev_match.ref_start == m.ref_start) {
                continue; //panic!("redundant partial hit encountered");
            }
            let mut is_added = false;
            for o in &mut open_nams {
                // Extend NAM
                if (o.query_prev_match_startpos <= m.query_start)
                    && (m.query_start <= o.query_end)
                    && (o.ref_prev_match_startpos <= m.ref_start)
                    && (m.ref_start <= o.ref_end)
                {
                    if (m.query_end > o.query_end) && (m.ref_end > o.ref_end) {
                        o.query_end = m.query_end;
                        o.ref_end = m.ref_end;
                        o.query_prev_match_startpos = m.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_match_startpos = m.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_matches += 1;
                        is_added = true;
                        break;
                    } else if (m.query_end <= o.query_end) && (m.ref_end <= o.ref_end) {
                        o.query_prev_match_startpos = m.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_match_startpos = m.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_matches += 1;
                        is_added = true;
                        break;
                    }
                }
            }
            // Add the match to open matches
            if !is_added {
                open_nams.push(Nam {
                    nam_id: nams.len() + open_nams.len(),
                    query_start: m.query_start,
                    query_end: m.query_end,
                    ref_start: m.ref_start,
                    ref_end: m.ref_end,
                    ref_id: *ref_id,
                    query_prev_match_startpos: m.query_start,
                    ref_prev_match_startpos: m.ref_start,
                    n_matches: 1,
                    is_revcomp,
                    score: 0.0,
                });
            }

            // Only filter if we have advanced at least k nucleotides
            if m.query_start > prev_q_start + k {

                // Output all NAMs from open_matches to final_nams that the current match have passed
                for n in &open_nams {
                    if n.query_end < m.query_start {
                        let n_max_span = max(n.query_span(), n.ref_span());
                        let n_min_span = min(n.query_span(), n.ref_span());
                        let n_score =
                            if 2 * n_min_span > n_max_span {
                                // this is really just n_matches * (min_span - (offset_in_span) ) );
                                n.n_matches * (2 * n_min_span - n_max_span)
                            } else {
                                1
                            } as u32;
//                        n_score = n.n_matches * n.query_span();
                        let mut nam = n.clone();
                        nam.score = n_score as f32;
                        nams.push(nam);
                    }
                }

                // Remove all NAMs from open_matches that the current match have passed
                let c = m.query_start;

                open_nams.retain(|x| x.query_end >= c);

                prev_q_start = m.query_start;
            }
        }

        // Add all current open_matches to final NAMs
        for mut n in open_nams {
            let n_max_span = max(n.query_span(), n.ref_span());
            let n_min_span = min(n.query_span(), n.ref_span());
            n.score =
                if 2 * n_min_span > n_max_span {
                    (n.n_matches * (2 * n_min_span - n_max_span)) as f32
                } else {
                    1.0
                };
            nams.push(n);
        }
    }
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
    // we need two extra checks for this - hopefully this will remove all the false matches we see
    // (true hash collisions should be very few)
    let read_len = read.len();
    let q_start_tmp = read_len - nam.query_end;
    let q_end_tmp = read_len - nam.query_start;
    // false reverse match, change coordinates in nam to forward
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

fn hits_to_matches(hits: &[Hit], index: &StrobemerIndex) -> DefaultHashMap<usize, Vec<Match>> {
    let mut matches_map = DefaultHashMap::default();
    matches_map.reserve(100);

    for hit in hits {
        if hit.is_partial {
            add_to_matches_map_partial(&mut matches_map, hit.query_start, hit.query_end, index, hit.position);
        } else {
            add_to_matches_map_full(&mut matches_map, hit.query_start, hit.query_end, index, hit.position);
        }
    }

    matches_map
}

/// Obtain NAMs for a sequence record, doing rescue if needed.
///
/// NAMs are returned sorted by decreasing score
pub fn get_nams_by_chaining(
    sequence: &[u8],
    index: &StrobemerIndex,
    chainer: &Chainer,
    rescue_level: usize,
    mcs_strategy: McsStrategy,
    rng: &mut Rng,
) -> (NamDetails, Vec<Nam>) {
    let timer = Instant::now();
    let query_randstrobes = mapper::randstrobes_query(sequence, &index.parameters);
    let n_randstrobes = query_randstrobes[0].len() + query_randstrobes[1].len();
    let time_randstrobes = timer.elapsed().as_secs_f64();

    trace!("we have {} + {} randstrobes", query_randstrobes[0].len(), query_randstrobes[1].len());
    let (mut nam_details, mut nams) = chainer.get_chains(&query_randstrobes, index, rescue_level, mcs_strategy);

    let timer = Instant::now();
    nams.sort_by_key(|k| -(k.score as i32));
    shuffle_top_nams(&mut nams, rng);
    nam_details.time_sort_nams = timer.elapsed().as_secs_f64();
    nam_details.time_randstrobes = time_randstrobes;

    if log::log_enabled!(Trace) {
        trace!("Found {} NAMs (rescue done: {})", nams.len(), nam_details.nam_rescue);
        let mut printed = 0;
        for nam in &nams {
            if nam.n_matches > 1 || printed < 10 {
                trace!("- {}", nam);
                printed += 1;
            }
        }
        if printed < nams.len() {
            trace!("+ {} single-anchor chains", nams.len() - printed);
        }
    }

    (nam_details, nams)
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
