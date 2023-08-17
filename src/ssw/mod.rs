// partially based on ssw_cpp.h/.cpp by Wan-Ping Lee and Mengyao Zhao

use std::marker::PhantomData;
use crate::cigar::Cigar;

mod raw;

  // @function Align the query againt the reference.
  //           [NOTICE] The reference won't replace the reference
  //                      set by SetReferenceSequence.
  // @param    query     The query sequence.
  // @param    ref       The reference sequence.
  //                     [NOTICE] It is not necessary null terminated.
  // @param    ref_len   The length of the reference sequence.
  // @param    filter    The filter for the alignment.
  // @param    alignment The container contains the result.
  // @param    maskLen   The distance between the optimal and suboptimal alignment ending position will >= maskLen. We suggest to
  //                     use readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function
  //                     will NOT return the suboptimal alignment information.
  // @return   If the alignment path is accurate (or has missing part). 0: accurate; 1: banded_sw is totally failed; 2: banded_sw returned path has missing part
  // =========
  // uint16_t Align(const char* query, const char* ref, const int& ref_len,
  //            const Filter& filter, Alignment* alignment, const int32_t maskLen) const;


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
struct Profile<'a> {
    profile: *mut raw::s_profile,
    query: &'a [i8],
    score_matrix: &'a [i8]
}

#[derive(Debug)]
pub struct SswAlignment {
    pub score: u16,
    //pub sw_score_next_best: u16, // The next best alignment score
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    //pub ref_end_next_best: usize,
    //pub mismatches: i32,
    pub cigar: Cigar,
}

struct Alignment<'a> {
    raw: *mut raw::s_align,
    _align: PhantomData<&'a raw::s_align>,
}

impl<'a> Alignment<'a> {
    fn is_valid(&self) -> bool {
        unsafe { (*self.raw).flag == 0 }
    }
}

impl<'a> Drop for Alignment<'a> {
    fn drop(&mut self) {
        unsafe { raw::align_destroy(self.raw); }
    }
}

impl<'a> From<Alignment<'a>> for SswAlignment {
    fn from(alignment: Alignment) -> SswAlignment {
        let raw = unsafe { alignment.raw.as_ref().expect("hm") };
        let cigar_slice = unsafe { std::slice::from_raw_parts(raw.cigar, raw.cigar_length as usize) };
        let cigar = Cigar::try_from(cigar_slice).expect("Invalid CIGAR");

            /*al->cigar.clear();
  al->cigar_string.clear();

  if (s_al.cigarLen > 0) {
    std::ostringstream cigar_string;
    if (al->query_begin > 0) {
      uint32_t cigar = to_cigar_int(al->query_begin, 'S');
      al->cigar.push_back(cigar);
      cigar_string << al->query_begin << 'S';
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      al->cigar.push_back(s_al.cigar[i]);
      cigar_string << cigar_int_to_len(s_al.cigar[i]) << cigar_int_to_op(s_al.cigar[i]);
    }

    int end = query_len - al->query_end - 1;
    if (end > 0) {
      uint32_t cigar = to_cigar_int(end, 'S');
      al->cigar.push_back(cigar);
      cigar_string << end << 'S';
    }

  } // end if
}
*/
/*  ConvertAlignment(*s_al, query_len, alignment);
  alignment->mismatches = CalculateNumberMismatch(&*alignment, translated_ref, translated_query, query_len);
*/

        SswAlignment {
            score: raw.score1,
            ref_start: raw.ref_begin1 as usize,
            ref_end: raw.ref_end1 as usize + 1,
            query_start: raw.read_begin1 as usize,
            query_end: raw.read_end1 as usize + 1,
            cigar,
        }
    }
}

impl<'a> Profile<'a> {
    fn new(translated_query: &'a [i8], score_matrix: &'a [i8]) -> Self {
        // TODO should return an error if query.is_empty()
        let score_size = 2;
        let profile = unsafe {
            // TODO hardcoded 5
            raw::ssw_init(translated_query.as_ptr(), translated_query.len() as i32, score_matrix.as_ptr(), 5i32, score_size)
        };
        Profile { profile, query: translated_query, score_matrix }
    }

    fn align(&self, translated_refseq: &[i8], gap_open_penalty: u8, gap_extend_penalty: u8, flag: u8, score_filter: u16, distance_filter: i32, mask_len: i32) -> Alignment {
        let alignment = unsafe {
            raw::ssw_align(self.profile, translated_refseq.as_ptr(), translated_refseq.len() as i32, gap_open_penalty, gap_extend_penalty, flag, score_filter, distance_filter, mask_len)
        };
        Alignment { raw: alignment, _align: PhantomData }
    }
}


impl<'a> Drop for Profile<'a> {
    fn drop(&mut self) {
        unsafe { raw::init_destroy(self.profile); }
    }
}

pub struct Aligner {
    score_matrix: Vec<i8>,
    match_score: u8,
    mismatch_penalty: u8,
    gap_open_penalty: u8,
    gap_extend_penalty: u8,
}

impl Aligner {
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
        Aligner { score_matrix, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty }
    }

    pub fn align(&self, query: &[u8], refseq: &[u8]) -> Option<SswAlignment> {
        if query.is_empty() {
            return None;
        }

        let query = translate(query);
        let refseq = translate(refseq);
        let profile = Profile::new(&query, &self.score_matrix);
        let flag = 0x0f;
        let score_filter = 0;
        let distance_filter = i32::MAX;
        let mask_len = std::cmp::max(query.len() / 2, 15);

        let alignment = profile.align(&refseq, self.gap_open_penalty, self.gap_extend_penalty, flag, score_filter, distance_filter, mask_len as i32);
        if !alignment.is_valid() {
            return None;
        }

        Some(SswAlignment::from(alignment))
    }
}
