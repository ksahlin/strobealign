use std::cmp::{max, min, Reverse};
use crate::aligner::{AlignmentInfo, hamming_align, hamming_distance};
use crate::cigar::Cigar;
use crate::fasta::RefSequence;
use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::{find_nams, Nam, reverse_nam_if_needed};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;
use crate::sam::{REVERSE, SamRecord, SECONDARY};
use crate::read::Read;
use crate::aligner::Aligner;
use crate::fastq::SequenceRecord;

enum CigarMode {
    M, Eqx
}

pub struct MappingParameters {
    r: usize, // { 150 },
    max_secondary: usize,
    dropoff_threshold: f32,
    rescue_level: usize,
    max_tries: usize,
    // rescue_cutoff,
    cigar_mode: CigarMode,
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
            // rescue_cutoff: ...,
            cigar_mode: CigarMode::M,
            output_unmapped: true,
            // details: false,
        }
    }
}

#[derive(Clone,Default)]
struct Alignment {
    reference_id: usize,
    ref_start: usize,
    cigar: Cigar,
    edit_distance: u32,
    soft_clipped: u32,
    score: u32,
    length: usize,
    is_revcomp: bool,
    is_unaligned: bool,
    /// Whether a gapped alignment function was used to obtain this alignment
    /// (even if true, the alignment can still be without gaps)
    gapped: bool,
}

impl Alignment {
    /*fn new() -> Self {
        Alignment {

        }
    }*/

    fn global_edit_distance(&self) -> u32 {
        self.edit_distance + self.soft_clipped
    }
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub start: usize,
    pub end: usize,
    pub is_reverse: bool,
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
                    is_reverse
                }
            );
        }
    }

    randstrobes
}

// TODO details
fn make_sam_record(alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool) -> SamRecord {
    let mut flags = 0;

    if !alignment.is_unaligned && alignment.is_revcomp {
        flags |= REVERSE;
    }
    if !is_primary {
        mapq = 255;
        flags |= SECONDARY;
    }

    SamRecord {
        query_name: record.name.clone(),
        flags,
        reference_name: Some(references[alignment.reference_id].name.clone()),
        pos: Some(alignment.ref_start as u32),
        mapq,
        cigar: Some(alignment.cigar.clone()),
        mate_reference_name: None,
        mate_pos: None,
        template_len: None,
        query_sequence: Some(record.sequence.clone()),
        query_qualities: Some(record.qualities.clone()),
        edit_distance: alignment.edit_distance,
        alignment_score: alignment.score,
        // TODO details: details
    }
}

pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    aligner: &Aligner,
) -> Vec<SamRecord> {
    //Details details;
    //Timer strobe_timer;
    let query_randstrobes = randstrobes_query(&record.sequence, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();

    // Timer nam_timer;
    let (nonrepetitive_fraction, mut nams) = find_nams(&query_randstrobes, index);
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

    if nams.is_empty() {
        return Vec::new();
    }
    let mut sam_records = Vec::new();
    let mut alignments = Vec::new();
    let nam_max = nams[0];
    let mut best_edit_distance = u32::MAX;
    let mut best_score = 0;
    let mut second_best_score = 0;
    let mut best_alignment = Alignment::default();

    let k = index.parameters.syncmer.k;
    let read = Read::new(&record.sequence);

    for (tries, nam) in nams.iter_mut().enumerate() {
        let score_dropoff = nam.n_hits as f32 / nam_max.n_hits as f32;

        // TODO iterate over slice of nams instead of tracking tries
        if tries >= mapping_parameters.max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < mapping_parameters.dropoff_threshold {
            break;
        }
        let consistent_nam = reverse_nam_if_needed(nam, &read, references, k);
        // details.nam_inconsistent += !consistent_nam;
        let alignment = get_alignment(aligner, nam, references, &read, consistent_nam);
        // details.tried_alignment += 1;
        // details.gapped += sam_aln.gapped;

        if alignment.score > best_score {
            second_best_score = best_score;
            best_score = alignment.score;
            best_alignment = alignment.clone();
            if mapping_parameters.max_secondary == 0 {
                best_edit_distance = best_alignment.global_edit_distance();
            }
        } else if alignment.score > second_best_score {
            second_best_score = alignment.score;
        }
        if mapping_parameters.max_secondary > 0 {
            alignments.push(alignment);
        }
    }
    let mapq = (60.0 * (best_score - second_best_score) as f32 / best_score as f32) as u8;

    if mapping_parameters.max_secondary == 0 {
        sam_records.push(
            make_sam_record(&best_alignment, references, record, mapq, true) // TODO details
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
            // TODO details
            sam_records.push(make_sam_record(alignment, references, record, mapq, is_primary));
        }
    }
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;
    sam_records
}

/*
 Extend a NAM so that it covers the entire read and return the resulting
 alignment.
*/
fn get_alignment(
    aligner: &Aligner,
    nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    consistent_nam: bool,
) -> Alignment {
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
        info = aligner.align(query, segment).expect("alignment failed");
        result_ref_start = ref_start + info.ref_start;
    }
    // TODO alignment cigar does not contain soft clipping
    let softclipped = info.query_start + (query.len() - info.query_end);
    Alignment {
        cigar: info.cigar.clone(),
        edit_distance: info.edit_distance,
        soft_clipped: softclipped as u32,
        score: info.score,
        ref_start: result_ref_start,
        length: info.ref_span(),
        is_revcomp: nam.is_revcomp,
        is_unaligned: false,
        reference_id: nam.ref_id,
        gapped,
    }
}
