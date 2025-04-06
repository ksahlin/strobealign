//! Bindings for StripedSmithWaterman (SSW) library

#[repr(C)]
#[allow(non_camel_case_types)]
#[derive(Debug)]
pub struct s_align {
    pub score1: u16,
    pub score2: u16,
    pub ref_begin1: i32,  // is -1 when not available
    pub ref_end1: i32,
    pub read_begin1: i32,  // is -1 when not available
    pub read_end1: i32,
    pub ref_end2: i32,
    pub cigar: *const u32,
    pub cigar_length: i32,
    pub flag: u16,
}

#[repr(C)] // TODO possibly not needed as this is an opaque struct
#[allow(non_camel_case_types)]
pub struct s_profile {
    _data: [u8; 0],
    _marker:
        core::marker::PhantomData<(*mut u8, core::marker::PhantomPinned)>,
}

#[link(name = "ssw")]
extern "C" {
    pub fn ssw_init(read: *const i8, read_length: i32, mat: *const i8, n: i32, score_size: i8) -> *mut s_profile;

    pub fn init_destroy(p: *mut s_profile);

    pub fn ssw_align(
        prof: *const s_profile,
        refseq: *const i8,
        ref_length: i32,
        weight_gap_open: u8,
        weight_gap_extend: u8,
        flag: u8,
        filters: u16,
        filterd: i32,
        mask_length: i32
    ) -> *mut s_align;

    pub fn align_destroy(a: *mut s_align);
}
