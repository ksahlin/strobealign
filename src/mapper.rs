use std::cmp::{min, Reverse};
use crate::aligner::{AlignmentInfo, hamming_align, hamming_distance};
use crate::cigar::{Cigar, CigarOperation};
use crate::fasta::RefSequence;
use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::{find_nams, find_nams_rescue, Nam, reverse_nam_if_needed};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;
use crate::sam::{SamRecord, REVERSE, SECONDARY, UNMAP};
use crate::read::Read;
use crate::aligner::Aligner;
use crate::details::Details;
use crate::fastq::SequenceRecord;


pub struct MappingParameters {
    r: usize,
    max_secondary: usize,
    dropoff_threshold: f32,
    rescue_level: usize,
    max_tries: usize,
    output_unmapped: bool,
}

impl Default for MappingParameters {
    fn default() -> Self {
        MappingParameters {
            r: 150,
            max_secondary: 0,
            dropoff_threshold: 0.5,
            rescue_level: 2,
            max_tries: 20,
            output_unmapped: true,
        }
    }
}

#[derive(Clone,Default)]
struct Alignment {
    reference_id: usize,
    ref_start: usize,
    cigar: Cigar,
    edit_distance: usize,
    soft_clip_left: usize,
    soft_clip_right: usize,
    score: u32,
    length: usize,
    is_revcomp: bool,
    is_unaligned: bool, // TODO get rid of this
    /// Whether a gapped alignment function was used to obtain this alignment
    /// (even if true, the alignment can still be without gaps)
    gapped: bool,
}

impl Alignment {
    fn global_edit_distance(&self) -> usize {
        self.edit_distance + self.soft_clip_left + self.soft_clip_right
    }
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub start: usize,
    pub end: usize,
    pub is_revcomp: bool,
}

    /// Generate randstrobes for a query sequence and its reverse complement.
pub fn randstrobes_query(seq: &[u8], parameters: &IndexParameters) -> Vec<QueryRandstrobe> {
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
                    is_revcomp: is_reverse
                }
            );
        }
    }

    randstrobes
}

/// Conversion of an Alignment into a SamRecord
#[derive(Default)]
pub struct SamOutput {
    cigar_eqx: bool,
    details: bool,
    rg_id: Option<String>,
}

impl SamOutput {
    pub fn new(details: bool, cigar_eqx: bool, rg_id: Option<String>) -> Self {
        SamOutput {
            cigar_eqx,
            details,
            rg_id,
        }
    }

    /// Convert the alignment into a SamRecord
    fn make_record(&self, alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool, details: Details) -> SamRecord {
        let mut flags = 0;

        if !alignment.is_unaligned && alignment.is_revcomp {
            flags |= REVERSE;
        }
        if !is_primary {
            mapq = 255;
            flags |= SECONDARY;
        }

        let query_sequence = if alignment.is_revcomp {
            reverse_complement(&record.sequence)
        } else {
            record.sequence.clone()
        };
        let query_qualities = if alignment.is_revcomp {
            let mut rev = record.qualities.clone();
            rev.reverse();
            rev
        } else {
            record.qualities.clone()
        };
        let mut cigar = Cigar::new();
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_left);
        cigar.extend(&alignment.cigar);
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_right);

        let reference_name = Some(references[alignment.reference_id].name.clone());
        let details = if self.details { Some(details) } else { None };

        let cigar = if self.cigar_eqx { Some(cigar) } else { Some(cigar.with_m()) };
        SamRecord {
            query_name: record.name.clone(),
            flags,
            reference_name,
            pos: Some(alignment.ref_start as u32),
            mapq,
            cigar,
            query_sequence: Some(query_sequence),
            query_qualities: Some(query_qualities),
            edit_distance: Some(alignment.edit_distance as u32),
            alignment_score: Some(alignment.score),
            details,
            rg_id: self.rg_id.clone(),
            ..SamRecord::default()
        }
    }

    fn make_unmapped_record(&self, record: &SequenceRecord, details: Details) -> SamRecord {
        let details = if self.details { Some(details) } else { None };
        SamRecord {
            query_name: record.name.clone(),
            flags: UNMAP,
            query_sequence: Some(record.sequence.clone()),
            query_qualities: Some(record.qualities.clone()),
            details,
            rg_id: self.rg_id.clone(),
            ..SamRecord::default()
        }
    }
}

pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    sam_output: &SamOutput,
    aligner: &Aligner,
) -> Vec<SamRecord> {
    let mut details = Details::default();
    //Timer strobe_timer;
    let query_randstrobes = randstrobes_query(&record.sequence, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();

    // Timer nam_timer;
    let (nonrepetitive_fraction, mut nams) = find_nams(&query_randstrobes, index, index.filter_cutoff);
    // statistics.tot_find_nams += nam_timer.duration();

    if mapping_parameters.rescue_level > 1 {

        // Timer rescue_timer;
        if nams.is_empty() || nonrepetitive_fraction < 0.7 {
            details.nam_rescue = true;
            nams = find_nams_rescue(&query_randstrobes, index, index.rescue_cutoff);
        }
        // statistics.tot_time_rescue += rescue_timer.duration();
    }
    details.nams = nams.len();
    // Timer nam_sort_timer;

    nams.sort_by_key(|&k| -(k.score as i32));
    // statistics.tot_sort_nams += nam_sort_timer.duration();

    // Timer extend_timer;
    // align_SE(
    //     aligner, sam, nams, record, index.parameters.syncmer.k,
    //     references, details, map_param.dropoff_threshold, map_param.maxTries,
    //     map_param.max_secondary
    // );

    if nams.is_empty() {
        return vec![sam_output.make_unmapped_record(record, details)];
    }
    let mut sam_records = Vec::new();
    let mut alignments = Vec::new();
    let nam_max = nams[0];
    let mut best_edit_distance = usize::MAX;
    let mut best_score = 0;
    let mut second_best_score = 0;
    let mut best_alignment = None;

    let k = index.parameters.syncmer.k;
    let read = Read::new(&record.sequence);

    for (tries, nam) in nams.iter_mut().enumerate() {
        let score_dropoff = nam.n_hits as f32 / nam_max.n_hits as f32;

        // TODO iterate over slice of nams instead of tracking tries
        if tries >= mapping_parameters.max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < mapping_parameters.dropoff_threshold {
            break;
        }
        let consistent_nam = reverse_nam_if_needed(nam, &read, references, k);
        details.nam_inconsistent += (!consistent_nam) as usize;
        let alignment = extend_seed(aligner, nam, references, &read, consistent_nam);
        if alignment.is_none() {
            continue;
        }
        let alignment = alignment.unwrap();
        details.tried_alignment += 1;
        details.gapped += alignment.gapped as usize;

        if alignment.score > best_score {
            second_best_score = best_score;
            best_score = alignment.score;
            best_alignment = Some(alignment.clone());
            if mapping_parameters.max_secondary == 0 {
                best_edit_distance = alignment.global_edit_distance();
            }
        } else if alignment.score > second_best_score {
            second_best_score = alignment.score;
        }
        if mapping_parameters.max_secondary > 0 {
            alignments.push(alignment);
        }
    }
    let mapq = ((60 * (best_score - second_best_score) + best_score - 1) / best_score) as u8;

    if best_alignment.is_none() {
        return vec![sam_output.make_unmapped_record(record, details)];
    }
    let best_alignment = best_alignment.unwrap();
    if mapping_parameters.max_secondary == 0 {
        sam_records.push(
            sam_output.make_record(&best_alignment, references, record, mapq, true, details)
        );
    } else {
        // Highest score first
        alignments.sort_by_key(|k| Reverse(k.score));

        let max_out = min(alignments.len(), mapping_parameters.max_secondary + 1);
        for (i, alignment) in alignments.iter().enumerate().take(max_out) {
            let is_primary = i == 0;
            if alignment.score - best_score > 2 * aligner.scores.mismatch as u32 + aligner.scores.gap_open as u32 {
                break;
            }
            // TODO .clone()
            sam_records.push(sam_output.make_record(alignment, references, record, mapq, is_primary, details.clone()));
        }
    }
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;
    sam_records
}

/// Extend a NAM so that it covers the entire read and return the resulting
/// alignment.
fn extend_seed(
    aligner: &Aligner,
    nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    consistent_nam: bool,
) -> Option<Alignment> {
    let query = if nam.is_revcomp { read.rc() } else { read.seq() };
    let refseq = &references[nam.ref_id].sequence;

    let projected_ref_start = nam.ref_start.saturating_sub(nam.query_start);
    let projected_ref_end = min(nam.ref_end + query.len() - nam.query_end, refseq.len());

    // TODO ugly
    let mut info = AlignmentInfo {
        cigar: Default::default(),
        edit_distance: 0,
        ref_start: 0,
        ref_end: 0,
        query_start: 0,
        query_end: 0,
        score: 0,
    };
    let mut result_ref_start = 0;
    let mut gapped = true;
    if projected_ref_start + query.len() == projected_ref_end && consistent_nam {
        let ref_segm_ham = &refseq[projected_ref_start..projected_ref_end];
        if let Some(hamming_dist) = hamming_distance(query, ref_segm_ham) {
            if (hamming_dist as f32 / query.len() as f32) < 0.05 {
                // ungapped worked fine, no need to do gapped alignment
                info = hamming_align(query, ref_segm_ham, aligner.scores.match_, aligner.scores.mismatch, aligner.scores.end_bonus).expect(
                    "hamming_dist was successful, this should be as well"
                );
                result_ref_start = projected_ref_start + info.ref_start;
                gapped = false;
            }
        }
    }
    if gapped {
        let ref_start = projected_ref_start.saturating_sub(50);
        let ref_end = min(nam.ref_end + 50, refseq.len());
        let segment = &refseq[ref_start..ref_end];
        info = aligner.align(query, segment)?;
        result_ref_start = ref_start + info.ref_start;
    }
    Some(Alignment {
        cigar: info.cigar.clone(),
        edit_distance: info.edit_distance,
        soft_clip_left: info.query_start,
        soft_clip_right: query.len() - info.query_end,
        score: info.score,
        ref_start: result_ref_start,
        length: info.ref_span(),
        is_revcomp: nam.is_revcomp,
        is_unaligned: false,
        reference_id: nam.ref_id,
        gapped,
    })
}
