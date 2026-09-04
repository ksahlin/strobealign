//! x86_64 backends: [`Avx2`] at 32 lanes and [`Sse41`] at 16.
//!
//! Both are written out longhand rather than through a macro, so that a reviewer checking a port
//! against the intrinsics reference reads the intrinsic itself.
//!
//! # Why SSE4.1 exists
//!
//! Not for pre-2013 hardware, but because it is **16 lanes on a machine that also has the
//! 32-lane kernel**, which makes it the only way to differentially test the `LANES != 32`
//! bookkeeping (the chunk counts, the `nm` packing, the band's edge cases) on hardware we can
//! actually run. NEON and simd128 are the same 16-lane shape, so a disagreement there is much
//! cheaper to diagnose once 16 lanes as such is ruled out. It also removes the hard panic on a
//! CPU without AVX2.
//!
//! # Where `#[target_feature]` is, and is not
//!
//! On [`Backend::fill`] alone, and deliberately **not** on the 21 operations; see
//! [`Backend::fill`] for why. The intrinsics are still `unsafe` to execute on a CPU without the
//! feature, and `select_backend` in `aligner.rs`, the only thing that names one of these types,
//! detects first.

use std::arch::x86_64::*;

use super::Backend;
use crate::simdaligner::kernel::U8Probe;

/// 32 `u8` lanes. The fastest path on any modern x86_64.
#[derive(Clone, Copy)]
pub struct Avx2;

unsafe impl Backend for Avx2 {
    const NAME: &'static str = "AVX2";

    fn is_available() -> bool {
        is_x86_feature_detected!("avx2")
    }
    const LANES: usize = 32;
    const LANES32: usize = 8;

    type V8 = __m256i;
    type V32 = __m256i;

    #[inline(always)]
    unsafe fn load8(p: *const u8) -> Self::V8 {
        unsafe { _mm256_loadu_si256(p as *const __m256i) }
    }

    #[inline(always)]
    unsafe fn store8(p: *mut u8, v: Self::V8) {
        unsafe { _mm256_storeu_si256(p as *mut __m256i, v) }
    }

    #[inline(always)]
    unsafe fn splat8(x: u8) -> Self::V8 {
        unsafe { _mm256_set1_epi8(x as i8) }
    }

    #[inline(always)]
    unsafe fn max8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm256_max_epu8(a, b) }
    }

    #[inline(always)]
    unsafe fn add8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm256_add_epi8(a, b) }
    }

    #[inline(always)]
    unsafe fn sub8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm256_sub_epi8(a, b) }
    }

    #[inline(always)]
    unsafe fn eq8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm256_cmpeq_epi8(a, b) }
    }

    #[inline(always)]
    unsafe fn and8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm256_and_si256(a, b) }
    }

    #[inline(always)]
    unsafe fn or8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm256_or_si256(a, b) }
    }

    #[inline(always)]
    unsafe fn blend8(if_false: Self::V8, if_true: Self::V8, mask: Self::V8) -> Self::V8 {
        unsafe { _mm256_blendv_epi8(if_false, if_true, mask) }
    }

    #[inline(always)]
    unsafe fn movemask8(v: Self::V8) -> u32 {
        unsafe { _mm256_movemask_epi8(v) as u32 }
    }

    #[inline(always)]
    unsafe fn load32(p: *const i32) -> Self::V32 {
        unsafe { _mm256_loadu_si256(p as *const __m256i) }
    }

    #[inline(always)]
    unsafe fn store32(p: *mut i32, v: Self::V32) {
        unsafe { _mm256_storeu_si256(p as *mut __m256i, v) }
    }

    #[inline(always)]
    unsafe fn splat32(x: i32) -> Self::V32 {
        unsafe { _mm256_set1_epi32(x) }
    }

    #[inline(always)]
    unsafe fn widen8to32(p: *const u8) -> Self::V32 {
        unsafe { _mm256_cvtepu8_epi32(_mm_loadl_epi64(p as *const __m128i)) }
    }

    #[inline(always)]
    unsafe fn add32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm256_add_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn sub32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm256_sub_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn max32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm256_max_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn eq32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm256_cmpeq_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn gt32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm256_cmpgt_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn movemask32(v: Self::V32) -> u32 {
        unsafe { _mm256_movemask_ps(_mm256_castsi256_ps(v)) as u32 }
    }

    #[target_feature(enable = "avx2")]
    unsafe fn fill<
        const TRACE: bool,
        const TRACK_LAST_ROW: bool,
        const TRACK_ROW_MAX: bool,
        const ARG_LARGEST: bool,
        const BANDED: bool,
    >(
        probe: &mut U8Probe<Self>,
        qlen: usize,
        rlen: usize,
        match_score: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
        w: usize,
    ) -> (i32, usize) {
        unsafe {
            probe.fill_generic::<TRACE, TRACK_LAST_ROW, TRACK_ROW_MAX, ARG_LARGEST, BANDED>(
                qlen,
                rlen,
                match_score,
                mismatch,
                gap_open,
                gap_extend,
                w,
            )
        }
    }
}

/// 16 `u8` lanes. See the module docs for why this is here.
#[derive(Clone, Copy)]
pub struct Sse41;

unsafe impl Backend for Sse41 {
    const NAME: &'static str = "SSE4.1";

    fn is_available() -> bool {
        is_x86_feature_detected!("sse4.1")
    }
    const LANES: usize = 16;
    const LANES32: usize = 4;

    type V8 = __m128i;
    type V32 = __m128i;

    #[inline(always)]
    unsafe fn load8(p: *const u8) -> Self::V8 {
        unsafe { _mm_loadu_si128(p as *const __m128i) }
    }

    #[inline(always)]
    unsafe fn store8(p: *mut u8, v: Self::V8) {
        unsafe { _mm_storeu_si128(p as *mut __m128i, v) }
    }

    #[inline(always)]
    unsafe fn splat8(x: u8) -> Self::V8 {
        unsafe { _mm_set1_epi8(x as i8) }
    }

    #[inline(always)]
    unsafe fn max8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm_max_epu8(a, b) }
    }

    #[inline(always)]
    unsafe fn add8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm_add_epi8(a, b) }
    }

    #[inline(always)]
    unsafe fn sub8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm_sub_epi8(a, b) }
    }

    #[inline(always)]
    unsafe fn eq8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm_cmpeq_epi8(a, b) }
    }

    #[inline(always)]
    unsafe fn and8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm_and_si128(a, b) }
    }

    #[inline(always)]
    unsafe fn or8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { _mm_or_si128(a, b) }
    }

    #[inline(always)]
    unsafe fn blend8(if_false: Self::V8, if_true: Self::V8, mask: Self::V8) -> Self::V8 {
        unsafe { _mm_blendv_epi8(if_false, if_true, mask) }
    }

    #[inline(always)]
    unsafe fn movemask8(v: Self::V8) -> u32 {
        unsafe { _mm_movemask_epi8(v) as u32 }
    }

    #[inline(always)]
    unsafe fn load32(p: *const i32) -> Self::V32 {
        unsafe { _mm_loadu_si128(p as *const __m128i) }
    }

    #[inline(always)]
    unsafe fn store32(p: *mut i32, v: Self::V32) {
        unsafe { _mm_storeu_si128(p as *mut __m128i, v) }
    }

    #[inline(always)]
    unsafe fn splat32(x: i32) -> Self::V32 {
        unsafe { _mm_set1_epi32(x) }
    }

    #[inline(always)]
    unsafe fn widen8to32(p: *const u8) -> Self::V32 {
        // Four bytes, not eight: `_mm_cvtepu8_epi32` widens the low four lanes. Read them as a
        // `u32` so this is one unaligned scalar load rather than a 16-byte load that could run
        // off the end of the buffer's LANES-byte padding.
        let four = unsafe { std::ptr::read_unaligned(p as *const u32) };
        unsafe { _mm_cvtepu8_epi32(_mm_cvtsi32_si128(four as i32)) }
    }

    #[inline(always)]
    unsafe fn add32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm_add_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn sub32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm_sub_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn max32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm_max_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn eq32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm_cmpeq_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn gt32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { _mm_cmpgt_epi32(a, b) }
    }

    #[inline(always)]
    unsafe fn movemask32(v: Self::V32) -> u32 {
        unsafe { _mm_movemask_ps(_mm_castsi128_ps(v)) as u32 }
    }

    #[target_feature(enable = "sse4.1")]
    unsafe fn fill<
        const TRACE: bool,
        const TRACK_LAST_ROW: bool,
        const TRACK_ROW_MAX: bool,
        const ARG_LARGEST: bool,
        const BANDED: bool,
    >(
        probe: &mut U8Probe<Self>,
        qlen: usize,
        rlen: usize,
        match_score: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
        w: usize,
    ) -> (i32, usize) {
        unsafe {
            probe.fill_generic::<TRACE, TRACK_LAST_ROW, TRACK_ROW_MAX, ARG_LARGEST, BANDED>(
                qlen,
                rlen,
                match_score,
                mismatch,
                gap_open,
                gap_extend,
                w,
            )
        }
    }
}
