use std::marker::PhantomData;
use crate::aligner::AlignmentInfo;
use crate::cigar::Cigar;

mod raw;


fn translate(query: &[u8]) -> Vec<i8> {
    static TABLE: [i8; 128] = [
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    ];

    query.iter().map(|c| TABLE[*c as usize]).collect()
}

#[derive(Debug)]
struct Profile {
    profile: *mut raw::s_profile,
}

/// Wrapper for the s_align struct that ensures memory is freed on drop()
struct SswAlignment<'a> {
    raw: *mut raw::s_align,
    _align: PhantomData<&'a raw::s_align>,
}

impl<'a> SswAlignment<'a> {
    fn is_valid(&self) -> bool {
        unsafe { (*self.raw).flag == 0 }
    }
}

impl<'a> Drop for SswAlignment<'a> {
    fn drop(&mut self) {
        unsafe { raw::align_destroy(self.raw); }
    }
}

impl<'a> From<SswAlignment<'a>> for AlignmentInfo {
    fn from(alignment: SswAlignment) -> AlignmentInfo {
        let raw = unsafe { alignment.raw.as_ref().unwrap() };
        let cigar_slice = unsafe { std::slice::from_raw_parts(raw.cigar, raw.cigar_length as usize) };
        let cigar = Cigar::try_from(cigar_slice).expect("Invalid CIGAR");

        AlignmentInfo {
            score: raw.score1 as u32,
            ref_start: raw.ref_begin1 as usize,
            ref_end: (raw.ref_end1 + 1) as usize,
            query_start: raw.read_begin1 as usize,
            query_end: (raw.read_end1 + 1) as usize,
            cigar,
            edit_distance: 0, // TODO
        }
    }
}

// Wrapper for s_profile that frees memory on drop()
impl Profile {
    fn new(translated_query: &[i8], score_matrix: &[i8]) -> Self {
        // TODO should return an error if query.is_empty()
        let score_size = 2;
        let profile = unsafe {
            // TODO hardcoded 5
            raw::ssw_init(translated_query.as_ptr(), translated_query.len() as i32, score_matrix.as_ptr(), 5i32, score_size)
        };
        Profile { profile }
    }

    #[allow(clippy::too_many_arguments)]
    fn align(&self, translated_refseq: &[i8], gap_open_penalty: u8, gap_extend_penalty: u8, flag: u8, score_filter: u16, distance_filter: i32, mask_len: i32) -> Option<SswAlignment> {
        let alignment;
        unsafe {
            alignment = raw::ssw_align(self.profile, translated_refseq.as_ptr(), translated_refseq.len() as i32, gap_open_penalty, gap_extend_penalty, flag, score_filter, distance_filter, mask_len);
            if (*alignment).ref_begin1 == -1 || (*alignment).read_begin1 == -1 {
                return None;
            }
        };
        Some(SswAlignment { raw: alignment, _align: PhantomData })
    }
}

impl Drop for Profile {
    fn drop(&mut self) {
        unsafe { raw::init_destroy(self.profile); }
    }
}

#[derive(Clone, Debug)]
pub struct SswAligner {
    score_matrix: Vec<i8>,
    gap_open_penalty: u8,
    gap_extend_penalty: u8,
}

impl SswAligner {
    pub fn new(match_score: u8, mismatch_penalty: u8, gap_open_penalty: u8, gap_extend_penalty: u8) -> Self {
        let mat = match_score as i8;
        let mis = -(mismatch_penalty as i8);
        let score_matrix = vec![
            mat, mis, mis, mis, mis,
            mis, mat, mis, mis, mis,
            mis, mis, mat, mis, mis,
            mis, mis, mis, mat, mis,
            mis, mis, mis, mis, mis,
        ];
        SswAligner { score_matrix, gap_open_penalty, gap_extend_penalty }
    }

    pub fn align(&self, query: &[u8], refseq: &[u8]) -> Option<AlignmentInfo> {
        if query.is_empty() {
            return None;
        }

        let translated_query = translate(query);
        let translated_refseq = translate(refseq);
        let profile = Profile::new(&translated_query, &self.score_matrix);
        let flag = 0x0f;
        let score_filter = 0;
        let distance_filter = i32::MAX;
        let mask_len = std::cmp::max(translated_query.len() / 2, 15);

        let alignment = profile.align(&translated_refseq, self.gap_open_penalty, self.gap_extend_penalty, flag, score_filter, distance_filter, mask_len as i32)?;
        if !alignment.is_valid() {
            return None;
        }
        let mut alignment = AlignmentInfo::from(alignment);
        alignment.cigar = alignment.cigar.with_eqx(&query[alignment.query_start..alignment.query_end], &refseq[alignment.ref_start..alignment.ref_end]);
        alignment.edit_distance = alignment.cigar.edit_distance();
        Some(alignment)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_ssw_align_no_result() {
        let aligner = SswAligner::new(2, 8, 12, 1);
        let query = b"TCTCTCCCTCTCTCTCTCTCCCTCCCTCTCTCTCCCTCTCTCTCTCTCTCTCCCTCCCTT";
        let refseq = b"GAGGGAGAGAGAGAGAGGGAGAGAGAGAGAGAG";
        let info = aligner.align(query, refseq);
        assert!(info.is_none());
    }
}
