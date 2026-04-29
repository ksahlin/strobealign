use std::cmp::Reverse;

use crate::{
    aligner::Aligner,
    details::{Details, NamDetails},
    io::fasta::RefSequence,
    mapper::{Alignment, extend_seed},
    math::normal_pdf,
    nam::{Nam, reverse_nam_if_needed},
    read::Read,
    shuffle::shuffle_best,
};

use bumpalo::collections::Vec as BumpVec;
use fastrand::Rng;
use memchr::memmem;

/// A scored alignment pair
#[derive(Debug, Clone)]
pub struct PairedAlignments {
    pub score: f64,
    pub alignment1: Alignment,
    pub alignment2: Alignment,
}

/// Align both reads and collect all plausible paired and individual alignments.
/// First tries NAMs that could form proper pairs, then falls back to
/// unpaired chains with mate rescue if the max_tries is not yet exhausted.
pub fn get_paired_alignment<'a>(
    aligner: &Aligner,
    mut nams1: BumpVec<Nam<'a>>,
    mut nams2: BumpVec<Nam<'a>>,
    max_tries: usize,
    dropoff: f32,
    references: &[RefSequence],
    read1: &Read,
    read2: &Read,
    details1: &mut Details,
    details2: &mut Details,
    mu: f32,
    sigma: f32,
    k: usize,
    rng: &mut Rng,
) -> (Vec<PairedAlignments>, Vec<Alignment>, Vec<Alignment>) {
    let mut paired_nams = get_paired_nams(
        &mut nams1,
        &mut nams2,
        mu,
        sigma,
        &details1.nam,
        &details2.nam,
    );
    let max_score = paired_nams.first().map_or(0.0, |p| p.score as f32);

    let mut paired_alignments = vec![];
    let mut alignments1 = vec![];
    let mut alignments2 = vec![];

    for p in paired_nams.iter_mut() {
        if p.score as f32 / max_score < dropoff || paired_alignments.len() == max_tries {
            break;
        }
        let a1 = align_nam(aligner, &mut p.nam1, references, read1, k, details1);
        let a2 = align_nam(aligner, &mut p.nam2, references, read2, k, details2);
        if let (Some(a1), Some(a2)) = (&a1, &a2) {
            paired_alignments.push(PairedAlignments {
                score: compute_combined_score(a1, a2, mu, sigma),
                alignment1: a1.clone(),
                alignment2: a2.clone(),
            });
        }
        alignments1.extend(a1);
        alignments2.extend(a2);
    }

    // Fallback to unpaired NAMs only if we didn't have enough alignment pairs
    if paired_alignments.len() < max_tries {
        let mut unpaired_nams = get_unpaired_nams(nams1, nams2, &paired_nams);
        let max_score = max_score.max(unpaired_nams.first().map_or(0.0, |u| u.nam.score));
        for u in unpaired_nams.iter_mut() {
            if u.nam.score / max_score < dropoff || paired_alignments.len() == max_tries {
                break;
            }
            let (a1, a2) = if u.read1 {
                (
                    align_nam(aligner, &mut u.nam, references, read1, k, details1),
                    rescue_align(aligner, &u.nam, references, read2, mu, sigma, k, details2),
                )
            } else {
                (
                    rescue_align(aligner, &u.nam, references, read1, mu, sigma, k, details1),
                    align_nam(aligner, &mut u.nam, references, read2, k, details2),
                )
            };
            if let (Some(a1), Some(a2)) = (&a1, &a2) {
                paired_alignments.push(PairedAlignments {
                    score: compute_combined_score(a1, a2, mu, sigma),
                    alignment1: a1.clone(),
                    alignment2: a2.clone(),
                });
            }
            alignments1.extend(a1);
            alignments2.extend(a2);
        }
    }
    alignments1.sort_unstable_by_key(|a| Reverse(a.score));
    alignments2.sort_unstable_by_key(|a| Reverse(a.score));
    shuffle_best(&mut alignments1, |a| a.score, rng);
    shuffle_best(&mut alignments2, |a| a.score, rng);

    if let Some(a1) = alignments1.first()
        && let Some(a2) = alignments2.first()
    {
        paired_alignments.push(PairedAlignments {
            score: compute_combined_score(a1, a2, mu, sigma),
            alignment1: a1.clone(),
            alignment2: a2.clone(),
        });
    }

    paired_alignments.sort_unstable_by(|a, b| b.score.total_cmp(&a.score));
    deduplicate_scored_pairs(&mut paired_alignments);
    shuffle_best(&mut paired_alignments, |a| a.score, rng);

    (paired_alignments, alignments1, alignments2)
}

/// Compute the combined score for a read pair.
/// Properly oriented pairs within the expected insert size get a log-normal
/// distance bonus
fn compute_combined_score(a1: &Alignment, a2: &Alignment, mu: f32, sigma: f32) -> f64 {
    let r1_r2 = a2.is_revcomp
        && !a1.is_revcomp
        && a1.ref_start <= a2.ref_start
        && (a2.ref_start - a1.ref_start) < (mu + 10.0 * sigma) as usize;
    let r2_r1 = a1.is_revcomp
        && !a2.is_revcomp
        && a2.ref_start <= a1.ref_start
        && (a1.ref_start - a2.ref_start) < (mu + 10.0 * sigma) as usize;

    if r1_r2 || r2_r1 {
        let x = a1.ref_start.abs_diff(a2.ref_start);
        a1.score as f64
            + a2.score as f64
            + (-20.0f64 + 0.001).max(normal_pdf(x as f32, mu, sigma).ln() as f64)
    } else {
        a1.score as f64 + a2.score as f64 - 20.0
    }
}

/// Remove consecutive identical alignment pairs and leave only the first.
fn deduplicate_scored_pairs(pairs: &mut Vec<PairedAlignments>) {
    pairs.dedup_by(|a, b| {
        a.alignment1.ref_start == b.alignment1.ref_start
            && a.alignment1.reference_id == b.alignment1.reference_id
            && a.alignment2.ref_start == b.alignment2.ref_start
            && a.alignment2.reference_id == b.alignment2.reference_id
    });
}

/// Align a NAM against the reference,
fn align_nam(
    aligner: &Aligner,
    nam: &mut Nam,
    references: &[RefSequence],
    read: &Read,
    k: usize,
    details: &mut Details,
) -> Option<Alignment> {
    let consistent = reverse_nam_if_needed(nam, read, references, k);
    details.inconsistent_nams += !consistent as usize;
    // Can do piecewise here
    let aln = extend_seed(aligner, nam, references, read, consistent, true);
    details.tried_alignment += 1;
    details.gapped += aln.as_ref().map_or(0, |a| a.gapped as usize);
    aln
}

/// Align a read to the reference given the mapping location of its mate.
fn rescue_align(
    aligner: &Aligner,
    mate_nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    mu: f32,
    sigma: f32,
    k: usize,
    details: &mut Details,
) -> Option<Alignment> {
    let read_len = read.len();

    let (r_tmp, ref_start, ref_end) = if mate_nam.is_revcomp {
        (
            read.seq(),
            mate_nam
                .projected_ref_start()
                .saturating_sub((mu + 5.0 * sigma) as usize),
            mate_nam.projected_ref_start() + read_len / 2, // at most half read overlap
        )
    } else {
        (
            read.rc(), // mate is rc since fr orientation
            (mate_nam.ref_end + read_len - mate_nam.query_end).saturating_sub(read_len / 2), // at most half read overlap
            mate_nam.ref_end + read_len - mate_nam.query_end + (mu + 5.0 * sigma) as usize,
        )
    };

    let ref_len = references[mate_nam.ref_id].sequence.len();
    let ref_start = ref_start.min(ref_len);
    let ref_end = ref_end.min(ref_len);

    if ref_end < ref_start + k {
        //        std::cerr << "RESCUE: Caught Bug3! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
        return None;
    }
    let ref_segm = &references[mate_nam.ref_id].sequence[ref_start..ref_end];

    if !has_shared_substring(r_tmp, ref_segm, k) {
        return None;
    }

    let info = aligner.align(r_tmp, ref_segm)?;
    details.mate_rescue += 1;

    Some(Alignment {
        reference_id: mate_nam.ref_id,
        ref_start: ref_start + info.ref_start,
        edit_distance: info.edit_distance,
        soft_clip_left: info.query_start,
        soft_clip_right: read_len - info.query_end,
        score: info.score,
        length: info.ref_span(),
        cigar: info.cigar,
        is_revcomp: !mate_nam.is_revcomp,
        gapped: true,
        rescued: true,
    })
}

/// Determine (roughly) whether the read sequence has some l-mer (with l = k*2/3)
/// in common with the reference sequence
fn has_shared_substring(read_seq: &[u8], ref_seq: &[u8], k: usize) -> bool {
    let sub_size = 2 * k / 3;
    let step_size = k / 3;
    for i in (0..read_seq.len().saturating_sub(sub_size)).step_by(step_size) {
        let submer = &read_seq[i..i + sub_size];
        if memmem::find(ref_seq, submer).is_some() {
            return true;
        }
    }

    false
}

/// Return true if the two alignments form a proper pair:
/// on the same reference, in opposite orientations, within the expected insert size.
pub fn is_proper_pair(a1: &Alignment, a2: &Alignment, mu: f32, sigma: f32) -> bool {
    let dist = a2.ref_start as isize - a1.ref_start as isize;
    let same_reference = a1.reference_id == a2.reference_id;
    let r1_r2 = !a1.is_revcomp && a2.is_revcomp && dist >= 0; // r1 ---> <---- r2
    let r2_r1 = !a2.is_revcomp && a1.is_revcomp && dist <= 0; // r2 ---> <---- r1
    let rel_orientation_good = r1_r2 || r2_r1;
    let insert_good = dist.unsigned_abs() <= (mu + sigma * 6.0) as usize;
    same_reference && insert_good && rel_orientation_good
}

/// Properly paired nams
#[derive(Debug)]
pub struct PairedNams<'a> {
    pub nam1: Nam<'a>,
    pub nam2: Nam<'a>,
    pub score: f64,
}

/// Build all plausible forward/revcomp mapping pairings
pub fn get_paired_nams<'a, 'b>(
    nams1: &'a mut [Nam<'b>],
    nams2: &'a mut [Nam<'b>],
    mu: f32,
    sigma: f32,
    details1: &NamDetails,
    details2: &NamDetails,
) -> Vec<PairedNams<'b>> {
    let mut paired_nams = vec![];
    if nams1.is_empty() || nams2.is_empty() {
        return paired_nams;
    }

    let (fwd1, rev1): (&'a mut [Nam<'b>], &'a mut [Nam<'b>]) =
        split_nams_by_orientation_checked(nams1, details1.both_orientations);
    let (fwd2, rev2): (&'a mut [Nam<'b>], &'a mut [Nam<'b>]) =
        split_nams_by_orientation_checked(nams2, details2.both_orientations);

    if !fwd1.is_empty() && !rev2.is_empty() {
        fwd1.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        rev2.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        paired_nams.extend(find_pairs(fwd1, rev2, mu, sigma, false));
    }

    if !fwd2.is_empty() && !rev1.is_empty() {
        fwd2.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        rev1.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        paired_nams.extend(find_pairs(fwd2, rev1, mu, sigma, true));
    }

    paired_nams.sort_unstable_by(|a, b| b.score.total_cmp(&a.score));
    paired_nams
}

/// Split nams into (forward, revcomp),
/// if only 1 orientation exists, returns it and a empty slice
fn split_nams_by_orientation_checked<'a, 'b>(
    nams: &'a mut [Nam<'b>],
    both: bool,
) -> (&'a mut [Nam<'b>], &'a mut [Nam<'b>]) {
    if both {
        split_nams_by_orientation(nams)
    } else if nams[0].is_revcomp {
        (&mut [], nams)
    } else {
        (nams, &mut [])
    }
}

/// In-place partition of NAMs by orientation:
/// forward on the left, revcomp on the right.
/// Returns two slices separating (forward, revcomp)
fn split_nams_by_orientation<'a, 'b>(
    nams: &'a mut [Nam<'b>],
) -> (&'a mut [Nam<'b>], &'a mut [Nam<'b>]) {
    let mut left = 0;
    let mut right = nams.len();

    while left < right {
        if nams[left].is_revcomp {
            right -= 1;
            nams.swap(left, right);
        } else {
            left += 1;
        }
    }

    nams.split_at_mut(left)
}

/// Find most forward/revcomp pairs using a two-pointer scan.
/// Assumes both slices are sorted by (ref_id, projected_ref_start).
fn find_pairs<'a>(
    fwd: &[Nam<'a>],
    rev: &[Nam<'a>],
    mu: f32,
    sigma: f32,
    swap_order: bool,
) -> Vec<PairedNams<'a>> {
    let mut out = Vec::new();
    let max_dist = (mu + 10.0 * sigma).ceil() as usize; // distance cutoff from insert size distribution
    let mut rev_ptr = 0;
    let mut last_paired = None;

    for f in fwd {
        // Advance revcomp pointer to the first possible candidate
        while rev_ptr < rev.len()
            && (rev[rev_ptr].ref_id < f.ref_id
                || rev[rev_ptr].ref_id == f.ref_id
                    && rev[rev_ptr].projected_ref_start() < f.projected_ref_start())
        {
            rev_ptr += 1;
        }
        if rev_ptr == rev.len() {
            break;
        }
        if rev[rev_ptr].ref_id > f.ref_id {
            continue;
        }

        // Scan window of revcomp nams within distance limit.
        let mut best = None;
        let mut i = rev_ptr;
        while i < rev.len()
            && rev[i].ref_id == f.ref_id
            && (rev[i].projected_ref_start() - f.projected_ref_start()) <= max_dist
        {
            let r = &rev[i];
            // The pairing score gets a bonus based on the reference distance of the two chosen nam
            // paired and from our current knowledge of the reference distance distribution
            let x = f.ref_start.abs_diff(r.ref_start);
            let score = f.score as f64
                + r.score as f64
                + 0.001f64.max((normal_pdf(x as f32, mu, sigma) + 1.0).ln() as f64);

            if best.is_none_or(|(_, highest_score)| score > highest_score) {
                best = Some((i, score));
            }
            i += 1;
        }

        // Highest scoring candidate
        let Some((best_id, score)) = best else {
            continue;
        };
        let r = &rev[best_id];

        // If the same revcomp nam was paired previously, keep only the better scoring pair.
        if let Some((last_id, prev_score)) = last_paired
            && last_id == best_id
        {
            if score <= prev_score {
                continue; // keep the previous better pair
            }
            out.pop(); // replace it 
        }

        out.push(if swap_order {
            PairedNams {
                nam1: r.clone(),
                nam2: f.clone(),
                score,
            }
        } else {
            PairedNams {
                nam1: f.clone(),
                nam2: r.clone(),
                score,
            }
        });

        last_paired = Some((best_id, score));
        rev_ptr = best_id;
    }

    out
}

/// Nam without any pairing partner
#[derive(Debug)]
pub struct UnpairedNam<'a> {
    pub nam: Nam<'a>,
    pub read1: bool,
}

/// Returns the nams that did not get paired
pub fn get_unpaired_nams<'a>(
    nams1: BumpVec<Nam<'a>>,
    nams2: BumpVec<Nam<'a>>,
    paired_nams: &[PairedNams],
) -> Vec<UnpairedNam<'a>> {
    let mut paired1 = vec![false; nams1.len()];
    let mut paired2 = vec![false; nams2.len()];

    for pair in paired_nams {
        paired1[pair.nam1.nam_id] = true;
        paired2[pair.nam2.nam_id] = true;
    }

    let mut unpaired_nams =
        Vec::with_capacity((nams1.len() + nams2.len()).saturating_sub(paired_nams.len()));

    for nam in nams1 {
        if !paired1[nam.nam_id] {
            unpaired_nams.push(UnpairedNam { nam, read1: true });
        }
    }
    for nam in nams2 {
        if !paired2[nam.nam_id] {
            unpaired_nams.push(UnpairedNam { nam, read1: false });
        }
    }

    unpaired_nams.sort_unstable_by(|a, b| b.nam.score.total_cmp(&a.nam.score));
    unpaired_nams
}

#[cfg(test)]
mod tests {
    use bumpalo::Bump;
    use bumpalo::collections::Vec as BumpVec;

    use super::*;
    use crate::{chainer::Anchor, mapper::Alignment, nam::Nam};

    fn make_nam<'a>(
        nam_id: usize,
        ref_id: usize,
        ref_start: usize,
        ref_end: usize,
        score: f32,
        is_revcomp: bool,
        arena: &'a mut Bump,
    ) -> Nam<'a> {
        let mut anchors = BumpVec::new_in(arena);
        anchors.push(Anchor {
            ref_id: 0,
            ref_start: 0,
            query_start: 0,
        });

        Nam {
            nam_id,
            ref_id,
            ref_start,
            ref_end,
            query_start: 0,
            query_end: ref_end - ref_start,
            anchors,
            matching_bases: ref_end - ref_start,
            score,
            is_revcomp,
        }
    }

    fn make_alignment(
        reference_id: usize,
        ref_start: usize,
        score: u32,
        is_revcomp: bool,
    ) -> Alignment {
        Alignment {
            reference_id,
            ref_start,
            score,
            is_revcomp,
            edit_distance: 0,
            soft_clip_left: 0,
            soft_clip_right: 0,
            length: 50,
            cigar: Default::default(),
            gapped: false,
            rescued: false,
        }
    }

    #[test]
    fn shared_substring_found() {
        assert!(has_shared_substring(b"ACGTACGTACGT", b"NNNNACGTACNNNN", 9));
    }

    #[test]
    fn shared_substring_not_found() {
        assert!(!has_shared_substring(b"ACGTACGTACGT", b"TTTTTTTTTTTTTT", 9));
    }

    #[test]
    fn shared_substring_empty() {
        assert!(!has_shared_substring(b"", b"ACGT", 9));
    }

    #[test]
    fn shared_substring_shorter_than_submer() {
        assert!(!has_shared_substring(b"ACG", b"ACGACGACG", 9));
    }

    #[test]
    fn get_paired_nams_empty() {
        let mut nams1: Vec<Nam> = vec![];
        let mut nams2: Vec<Nam> = vec![];
        let d1 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let d2 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let pairs = get_paired_nams(&mut nams1, &mut nams2, 300.0, 50.0, &d1, &d2);
        assert!(pairs.is_empty());
    }
}
