//! The portable fallback: the same kernel, with arrays standing in for vector registers.
//!
//! This is **not** an independently-written reference implementation of the DP. It is the exact
//! same kernel, with `[u8; 16]` where a real backend has a register. So it is a *portability*
//! fallback, not a correctness oracle: if the shared kernel had a bug, this would reproduce it
//! faithfully.
//!
//! 16 lanes rather than 1, matching the 128-bit backends: the width is what the kernel's chunk
//! bookkeeping and `nm` packing are parameterized on, so reusing a real width keeps this path on
//! the same code paths as SSE4.1 and NEON instead of on a degenerate one of its own. LLVM
//! auto-vectorizes most of these loops anyway, so "scalar" describes the source, not necessarily
//! the object code.

use super::Backend;
use crate::simdaligner::kernel::U8Probe;

/// The portable fallback backend. See the [module docs](self).
#[derive(Clone, Copy)]
pub struct Scalar;

unsafe impl Backend for Scalar {
    const NAME: &'static str = "scalar";

    fn is_available() -> bool {
        // Always. That is what makes it the last resort.
        true
    }
    const LANES: usize = 16;
    const LANES32: usize = 4;

    type V8 = [u8; 16];
    type V32 = [i32; 4];

    #[inline(always)]
    unsafe fn load8(p: *const u8) -> Self::V8 {
        unsafe { std::ptr::read_unaligned(p as *const [u8; 16]) }
    }

    #[inline(always)]
    unsafe fn store8(p: *mut u8, v: Self::V8) {
        unsafe { std::ptr::write_unaligned(p as *mut [u8; 16], v) }
    }

    #[inline(always)]
    unsafe fn splat8(x: u8) -> Self::V8 {
        [x; 16]
    }

    #[inline(always)]
    unsafe fn max8(a: Self::V8, b: Self::V8) -> Self::V8 {
        std::array::from_fn(|i| a[i].max(b[i]))
    }

    #[inline(always)]
    unsafe fn add8(a: Self::V8, b: Self::V8) -> Self::V8 {
        std::array::from_fn(|i| a[i].wrapping_add(b[i]))
    }

    #[inline(always)]
    unsafe fn sub8(a: Self::V8, b: Self::V8) -> Self::V8 {
        std::array::from_fn(|i| a[i].wrapping_sub(b[i]))
    }

    #[inline(always)]
    unsafe fn eq8(a: Self::V8, b: Self::V8) -> Self::V8 {
        std::array::from_fn(|i| if a[i] == b[i] { 0xFF } else { 0 })
    }

    #[inline(always)]
    unsafe fn and8(a: Self::V8, b: Self::V8) -> Self::V8 {
        std::array::from_fn(|i| a[i] & b[i])
    }

    #[inline(always)]
    unsafe fn or8(a: Self::V8, b: Self::V8) -> Self::V8 {
        std::array::from_fn(|i| a[i] | b[i])
    }

    #[inline(always)]
    unsafe fn blend8(if_false: Self::V8, if_true: Self::V8, mask: Self::V8) -> Self::V8 {
        // MSB-select, matching `_mm256_blendv_epi8`. Every call site passes a saturated mask, so
        // this agrees with a bit-select too.
        std::array::from_fn(|i| {
            if mask[i] & 0x80 != 0 {
                if_true[i]
            } else {
                if_false[i]
            }
        })
    }

    #[inline(always)]
    unsafe fn movemask8(v: Self::V8) -> u32 {
        let mut m = 0u32;
        for (i, lane) in v.iter().enumerate() {
            m |= ((lane >> 7) as u32) << i;
        }
        m
    }

    #[inline(always)]
    unsafe fn load32(p: *const i32) -> Self::V32 {
        unsafe { std::ptr::read_unaligned(p as *const [i32; 4]) }
    }

    #[inline(always)]
    unsafe fn store32(p: *mut i32, v: Self::V32) {
        unsafe { std::ptr::write_unaligned(p as *mut [i32; 4], v) }
    }

    #[inline(always)]
    unsafe fn splat32(x: i32) -> Self::V32 {
        [x; 4]
    }

    #[inline(always)]
    unsafe fn widen8to32(p: *const u8) -> Self::V32 {
        let b = unsafe { std::ptr::read_unaligned(p as *const [u8; 4]) };
        std::array::from_fn(|i| b[i] as i32)
    }

    #[inline(always)]
    unsafe fn add32(a: Self::V32, b: Self::V32) -> Self::V32 {
        std::array::from_fn(|i| a[i].wrapping_add(b[i]))
    }

    #[inline(always)]
    unsafe fn sub32(a: Self::V32, b: Self::V32) -> Self::V32 {
        std::array::from_fn(|i| a[i].wrapping_sub(b[i]))
    }

    #[inline(always)]
    unsafe fn max32(a: Self::V32, b: Self::V32) -> Self::V32 {
        std::array::from_fn(|i| a[i].max(b[i]))
    }

    #[inline(always)]
    unsafe fn eq32(a: Self::V32, b: Self::V32) -> Self::V32 {
        std::array::from_fn(|i| if a[i] == b[i] { -1 } else { 0 })
    }

    #[inline(always)]
    unsafe fn gt32(a: Self::V32, b: Self::V32) -> Self::V32 {
        std::array::from_fn(|i| if a[i] > b[i] { -1 } else { 0 })
    }

    #[inline(always)]
    unsafe fn movemask32(v: Self::V32) -> u32 {
        let mut m = 0u32;
        for (i, lane) in v.iter().enumerate() {
            m |= (((*lane as u32) >> 31) & 1) << i;
        }
        m
    }

    // No `#[target_feature]`: this backend uses no instruction the base target lacks, which is
    // the entire point of it. The trampoline is still here because the trait requires it, and
    // it costs nothing: it is one inlinable call per alignment.
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
