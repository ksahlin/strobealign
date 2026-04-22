use std::fmt::{Display, Formatter};
use std::time::Instant;

use fastrand::Rng;
use log::Level::Trace;
use log::trace;

use crate::chainer::Anchor;
use crate::chainer::Chainer;
use crate::details::NamDetails;
use crate::index::StrobemerIndex;
use crate::io::fasta::RefSequence;
use crate::mcsstrategy::McsStrategy;
use crate::read::Read;
use crate::seeding::randstrobes_query;
use crate::seeding::ry_equal;

/// Non-overlapping approximate match
#[derive(Clone, Debug)]
pub struct Nam {
    pub nam_id: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub n_matches: usize,
    pub matching_bases: usize,
    pub ref_id: usize,
    pub score: f32,
    pub is_revcomp: bool,
    pub anchors: Vec<Anchor>,
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

    /// Returns whether a NAM represents a consistent match between read and
    /// reference by comparing the nucleotide sequences of the first and last
    /// strobe (taking orientation into account).
    pub fn is_consistent(
        &self,
        read: &Read,
        references: &[RefSequence],
        k: usize,
        adna_mode: bool,
        ry_len: usize,
    ) -> bool {
        let ref_start_kmer = &references[self.ref_id].sequence[self.ref_start..self.ref_start + k];
        let ref_end_kmer = &references[self.ref_id].sequence[self.ref_end - k..self.ref_end];

        let seq = if self.is_revcomp {
            read.rc()
        } else {
            read.seq()
        };
        let read_start_kmer = &seq[self.query_start..self.query_start + k];
        let read_end_kmer = &seq[self.query_end - k..self.query_end];

        if adna_mode {
            ry_equal(ref_start_kmer, read_start_kmer, ry_len)
                && ry_equal(ref_end_kmer, read_end_kmer, ry_len)
        } else {
            ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer
        }
    }
}

impl Display for Nam {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Nam(ref_id={}, query: {}..{}, ref: {}..{}, rc={}, score={})",
            self.ref_id,
            self.query_start,
            self.query_end,
            self.ref_start,
            self.ref_end,
            self.is_revcomp as u8,
            self.score
        )?;
        Ok(())
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
pub fn reverse_nam_if_needed(
    nam: &mut Nam,
    read: &Read,
    references: &[RefSequence],
    k: usize,
    adna_mode: bool,
    ry_len: usize,
) -> bool {
    let ref_start_kmer = &references[nam.ref_id].sequence[nam.ref_start..nam.ref_start + k];
    let ref_end_kmer = &references[nam.ref_id].sequence[nam.ref_end - k..nam.ref_end];

    let kmers_match = |a: &[u8], b: &[u8]| -> bool {
        if adna_mode { ry_equal(a, b, ry_len) } else { a == b }
    };

    let (seq, seq_rc) = if nam.is_revcomp {
        (read.rc(), read.seq())
    } else {
        (read.seq(), read.rc())
    };
    let read_start_kmer = &seq[nam.query_start..nam.query_start + k];
    let read_end_kmer = &seq[nam.query_end - k..nam.query_end];
    if kmers_match(ref_start_kmer, read_start_kmer)
        && kmers_match(ref_end_kmer, read_end_kmer)
    {
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
    if kmers_match(ref_start_kmer, read_start_kmer)
        && kmers_match(ref_end_kmer, read_end_kmer)
    {
        nam.is_revcomp = !nam.is_revcomp;
        nam.query_start = q_start_tmp;
        nam.query_end = q_end_tmp;
        true
    } else {
        false
    }
}

/// Obtain NAMs for a sequence record, doing rescue if needed.
///
/// NAMs are returned sorted by decreasing score
pub fn get_nams_by_chaining(
    sequence: &[u8],
    index: &StrobemerIndex,
    chainer: &Chainer,
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    rng: &mut Rng,
) -> (NamDetails, Vec<Nam>) {
    let timer = Instant::now();
    let query_randstrobes = randstrobes_query(sequence, &index.parameters);
    let time_randstrobes = timer.elapsed().as_secs_f64();

    trace!(
        "we have {} + {} randstrobes",
        query_randstrobes[0].len(),
        query_randstrobes[1].len()
    );

    let (mut nam_details, mut nams) =
        chainer.get_chains(&query_randstrobes, index, rescue_distance, mcs_strategy);

    let timer = Instant::now();

    nams.sort_by(|a, b| b.score.total_cmp(&a.score));
    shuffle_top_nams(&mut nams, rng);
    nam_details.time_sort_nams = timer.elapsed().as_secs_f64();
    nam_details.time_randstrobes = time_randstrobes;

    if log::log_enabled!(Trace) {
        trace!("Found {} NAMs", nams.len());
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::chainer::{Chainer, ChainingParameters};
    use crate::index::StrobemerIndex;
    use crate::io::fasta::RefSequence;
    use crate::seeding::{DEFAULT_AUX_LEN, SeedingParameters};

    #[test]
    fn test_partial_rymer_consistency() {
        // Build a reference containing a known k-mer whose canonical form
        // is the reverse complement. Construct a 200bp reference with unique flanking sequence so
        // the k-mer at the target position is the only hit.
        let target_kmer = b"TTGGCACAACTTCTTCTACAT";
        let k = 21;
        let ry_len = 8;

        let left_flank = b"ACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC";
        let right_flank = b"TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";

        let mut ref_seq = Vec::new();
        ref_seq.extend_from_slice(left_flank);
        ref_seq.extend_from_slice(target_kmer);
        ref_seq.extend_from_slice(right_flank);
        let target_pos = left_flank.len();

        let references = vec![RefSequence {
            name: "test".to_string(),
            sequence: ref_seq.clone(),
        }];

        let parameters = SeedingParameters::from_read_length(
            40,
            Some(k),
            None,
            None,
            None,
            None,
            None,
            DEFAULT_AUX_LEN,
            true,
            Some(ry_len),
        )
        .unwrap();
        assert_eq!(parameters.syncmer.k, k);
        assert_eq!(parameters.ry_len, ry_len);

        let mut index = StrobemerIndex::new(&references, parameters.clone(), None);
        index.populate(0.0000001, 1);

        let chainer = Chainer::new(k, ChainingParameters::default());

        // Create a read that matches the target position but with C→T at
        // position 15 of the k-mer (a deamination-like change).
        // We include some flanking context from the reference so the read
        // is longer than k and produces multiple seeds.
        let deam_pos = 15; // within the target k-mer
        let read_start = target_pos;
        let read_len = 40;
        let mut read_seq: Vec<u8> = ref_seq[read_start..read_start + read_len].to_vec();
        assert_eq!(read_seq[deam_pos], b'C');
        read_seq[deam_pos] = b'T'; // C→T deamination

        let read = Read::new(&read_seq);

        let mut rng = Rng::with_seed(42);
        let (_details, nams) = get_nams_by_chaining(
            &read_seq,
            &index,
            &chainer,
            100,
            McsStrategy::Always,
            &mut rng,
        );

        // Find NAMs that overlap the correct reference position
        let correct_nams: Vec<&Nam> = nams
            .iter()
            .filter(|nam| {
                nam.ref_id == 0
                    && nam.ref_start <= target_pos + deam_pos
                    && nam.ref_end >= target_pos
            })
            .collect();

        assert!(
            !correct_nams.is_empty(),
            "Expected at least one NAM overlapping the correct reference position {}",
            target_pos,
        );

        // Check consistency of matches that contain at least one reverse complement
        for nam in &correct_nams {
            let consistent = nam.is_consistent(
                &read,
                &references,
                k,
                true,
                ry_len,
            );
            if nam.ref_start == left_flank.len() {            
                assert!(consistent, "Found inconsistent NAM {} for ry_len {}", nam, ry_len);
            }
        }
    }
}
