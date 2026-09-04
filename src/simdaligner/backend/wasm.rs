//! wasm32 simd128, 16 lanes.
//!
//! The easiest of the three ports: simd128 was specified with the x86 and NEON backends of real
//! codebases in view, and it has a native bitmask (`i8x16_bitmask` and `i32x4_bitmask`, lane
//! `i` in bit `i`), so unlike NEON there is nothing to synthesize. Every operation is one
//! intrinsic.
//!
//! # This backend is compile-time only
//!
//! wasm has no runtime feature detection, so there is no equivalent of the `is_x86_feature_detected!`
//! that picks between AVX2 and SSE4.1. Either the module is built with simd128 enabled and this
//! backend exists, or it is not and the aligner falls back to the portable kernel with a
//! warning. Enable it with:
//!
//! ```text
//! RUSTFLAGS="-C target-feature=+simd128" cargo build --target wasm32-unknown-unknown --release
//! ```
//!
//! The `cfg` on this module *is* the check: see `select_backend` in `aligner.rs`.

use std::arch::wasm32::*;

use super::Backend;
use crate::simdaligner::kernel::U8Probe;

/// 16 `u8` lanes on wasm simd128.
#[derive(Clone, Copy)]
pub struct Simd128;

unsafe impl Backend for Simd128 {
    const NAME: &'static str = "simd128";

    fn is_available() -> bool {
        // wasm has no runtime detection; the module cfg IS the check.
        true
    }
    const LANES: usize = 16;
    const LANES32: usize = 4;

    type V8 = v128;
    type V32 = v128;

    #[inline(always)]
    unsafe fn load8(p: *const u8) -> Self::V8 {
        unsafe { v128_load(p as *const v128) }
    }

    #[inline(always)]
    unsafe fn store8(p: *mut u8, v: Self::V8) {
        unsafe { v128_store(p as *mut v128, v) }
    }

    #[inline(always)]
    unsafe fn splat8(x: u8) -> Self::V8 {
        u8x16_splat(x)
    }

    #[inline(always)]
    unsafe fn max8(a: Self::V8, b: Self::V8) -> Self::V8 {
        u8x16_max(a, b)
    }

    #[inline(always)]
    unsafe fn add8(a: Self::V8, b: Self::V8) -> Self::V8 {
        // Wrapping. NOT `u8x16_add_sat`.
        u8x16_add(a, b)
    }

    #[inline(always)]
    unsafe fn sub8(a: Self::V8, b: Self::V8) -> Self::V8 {
        // Wrapping. NOT `u8x16_sub_sat`.
        u8x16_sub(a, b)
    }

    #[inline(always)]
    unsafe fn eq8(a: Self::V8, b: Self::V8) -> Self::V8 {
        u8x16_eq(a, b)
    }

    #[inline(always)]
    unsafe fn and8(a: Self::V8, b: Self::V8) -> Self::V8 {
        v128_and(a, b)
    }

    #[inline(always)]
    unsafe fn or8(a: Self::V8, b: Self::V8) -> Self::V8 {
        v128_or(a, b)
    }

    #[inline(always)]
    unsafe fn blend8(if_false: Self::V8, if_true: Self::V8, mask: Self::V8) -> Self::V8 {
        // Bit-select, like NEON's `vbslq_u8`; agrees with x86's MSB-select because every mask
        // the kernel passes is a saturated comparison result.
        v128_bitselect(if_true, if_false, mask)
    }

    #[inline(always)]
    unsafe fn movemask8(v: Self::V8) -> u32 {
        // Native, and already lane `i` -> bit `i`.
        i8x16_bitmask(v) as u32
    }

    #[inline(always)]
    unsafe fn load32(p: *const i32) -> Self::V32 {
        unsafe { v128_load(p as *const v128) }
    }

    #[inline(always)]
    unsafe fn store32(p: *mut i32, v: Self::V32) {
        unsafe { v128_store(p as *mut v128, v) }
    }

    #[inline(always)]
    unsafe fn splat32(x: i32) -> Self::V32 {
        i32x4_splat(x)
    }

    #[inline(always)]
    unsafe fn widen8to32(p: *const u8) -> Self::V32 {
        // Four bytes only: see the note in the NEON backend for why a wider load would run past
        // the delta buffers' padding on the group loop's last group.
        let four = unsafe { std::ptr::read_unaligned(p as *const u32) };
        // The splat puts the four bytes in every lane; only the low four are extended.
        u32x4_extend_low_u16x8(u16x8_extend_low_u8x16(u32x4_splat(four)))
    }

    #[inline(always)]
    unsafe fn add32(a: Self::V32, b: Self::V32) -> Self::V32 {
        i32x4_add(a, b)
    }

    #[inline(always)]
    unsafe fn sub32(a: Self::V32, b: Self::V32) -> Self::V32 {
        i32x4_sub(a, b)
    }

    #[inline(always)]
    unsafe fn max32(a: Self::V32, b: Self::V32) -> Self::V32 {
        i32x4_max(a, b)
    }

    #[inline(always)]
    unsafe fn eq32(a: Self::V32, b: Self::V32) -> Self::V32 {
        i32x4_eq(a, b)
    }

    #[inline(always)]
    unsafe fn gt32(a: Self::V32, b: Self::V32) -> Self::V32 {
        i32x4_gt(a, b)
    }

    #[inline(always)]
    unsafe fn movemask32(v: Self::V32) -> u32 {
        i32x4_bitmask(v) as u32
    }

    #[target_feature(enable = "simd128")]
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
