//! aarch64 NEON, 16 lanes.
//!
//! Most of the twenty-one operations are a one-to-one transcription of the AVX2 backend at half
//! the width. Three are not, and they are the whole of the interesting content of this file: the
//! two movemasks, which NEON does not have and which are synthesized below, and `widen8to32`,
//! which must not over-read the delta buffers' padding. `blend8` is also a bit-select rather than
//! an MSB-select, and agrees only because every mask the kernel passes it is saturated.
//!
//! # NEON has no movemask
//!
//! `_mm256_movemask_epi8` extracts one bit per lane in a single instruction. NEON has no
//! equivalent, and the kernel needs it on the hot path: with `TRACE` on (which every alignment
//! mode compiles), it is four extractions per chunk, plus four more per chunk under
//! `TRACK_ROW_MAX`.
//!
//! [`Neon::movemask8`] synthesizes it in four instructions: an arithmetic shift right by 7 to
//! turn each lane's MSB into a full `0x00`/`0xFF` mask, an `and` against per-lane bit weights
//! `1, 2, 4 ... 128` repeated across both halves, and one `addv` per half to sum each half's
//! weights into a byte. The two bytes are the low and high halves of the 16-bit result.
//!
//! The weights repeat every 8 lanes rather than running `1 << i` across all 16 because `addv`
//! reduces 8 lanes into a `u8`, and `1 + 2 + ... + 128` is exactly 255, the widest sum that
//! still fits. Splitting into halves is what keeps the reduction in 8-bit lanes; a 16-lane
//! version would need 16-bit lanes and a wider reduction to hold the same 16 bits.
//!
//! The cheaper `shrn`-by-4 trick is not used: it produces a nibble mask, where `arg_for_row`,
//! `flag_at` and `replay` all decode one bit per lane, so taking it would mean a per-arch flag
//! layout threaded through all four decoders.
//!
//! # Little-endian only
//!
//! This backend is gated on `target_endian = "little"`, which excludes `aarch64_be`. Two things
//! here depend on byte order and it is not obvious they cancel: `widen8to32` reads four bytes as
//! a scalar `u32` and then reinterprets the vector back to bytes, and big-endian NEON bitcasts
//! carry an implicit lane reversal. There is no big-endian machine here to test on, and both
//! reference backends are little-endian, so `aarch64_be` takes the portable kernel instead.

use std::arch::aarch64::*;

use super::Backend;
use crate::simdaligner::kernel::U8Probe;

/// 16 `u8` lanes on NEON. Baseline on aarch64, so there is nothing to detect.
#[derive(Clone, Copy)]
pub struct Neon;

/// Per-lane bit weights for [`Neon::movemask8`], repeating every 8 lanes so that each half sums
/// to at most 255. See the module docs.
const LANE_BITS: [u8; 16] = [1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128];

/// Per-lane bit weights for [`Neon::movemask32`].
const LANE_BITS32: [u32; 4] = [1, 2, 4, 8];

unsafe impl Backend for Neon {
    const NAME: &'static str = "NEON";

    fn is_available() -> bool {
        // The module cfg already required NEON; nothing is left to detect.
        true
    }
    const LANES: usize = 16;
    const LANES32: usize = 4;

    type V8 = uint8x16_t;
    type V32 = int32x4_t;

    #[inline(always)]
    unsafe fn load8(p: *const u8) -> Self::V8 {
        unsafe { vld1q_u8(p) }
    }

    #[inline(always)]
    unsafe fn store8(p: *mut u8, v: Self::V8) {
        unsafe { vst1q_u8(p, v) }
    }

    #[inline(always)]
    unsafe fn splat8(x: u8) -> Self::V8 {
        unsafe { vdupq_n_u8(x) }
    }

    #[inline(always)]
    unsafe fn max8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { vmaxq_u8(a, b) }
    }

    #[inline(always)]
    unsafe fn add8(a: Self::V8, b: Self::V8) -> Self::V8 {
        // Wrapping, matching `_mm256_add_epi8`. NOT `vqaddq_u8`, which saturates.
        unsafe { vaddq_u8(a, b) }
    }

    #[inline(always)]
    unsafe fn sub8(a: Self::V8, b: Self::V8) -> Self::V8 {
        // Wrapping, matching `_mm256_sub_epi8`. NOT `vqsubq_u8`.
        unsafe { vsubq_u8(a, b) }
    }

    #[inline(always)]
    unsafe fn eq8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { vceqq_u8(a, b) }
    }

    #[inline(always)]
    unsafe fn and8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { vandq_u8(a, b) }
    }

    #[inline(always)]
    unsafe fn or8(a: Self::V8, b: Self::V8) -> Self::V8 {
        unsafe { vorrq_u8(a, b) }
    }

    #[inline(always)]
    unsafe fn blend8(if_false: Self::V8, if_true: Self::V8, mask: Self::V8) -> Self::V8 {
        // `vbslq_u8` selects per *bit*, where `_mm256_blendv_epi8` selects per lane on the MSB.
        // Every mask the kernel passes here is a saturated `0x00`/`0xFF` comparison result, so
        // the two agree. A caller passing a partially-set mask would see them differ.
        unsafe { vbslq_u8(mask, if_true, if_false) }
    }

    #[inline(always)]
    unsafe fn movemask8(v: Self::V8) -> u32 {
        unsafe {
            // MSB of each lane -> 0x00 / 0xFF, by arithmetic shift right on signed lanes.
            let mask = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(v), 7));
            let bits = vandq_u8(mask, vld1q_u8(LANE_BITS.as_ptr()));
            // One `addv` per half: each sums to at most 255, so each lands in one byte.
            let lo = vaddv_u8(vget_low_u8(bits)) as u32;
            let hi = vaddv_u8(vget_high_u8(bits)) as u32;
            lo | (hi << 8)
        }
    }

    #[inline(always)]
    unsafe fn load32(p: *const i32) -> Self::V32 {
        unsafe { vld1q_s32(p) }
    }

    #[inline(always)]
    unsafe fn store32(p: *mut i32, v: Self::V32) {
        unsafe { vst1q_s32(p, v) }
    }

    #[inline(always)]
    unsafe fn splat32(x: i32) -> Self::V32 {
        unsafe { vdupq_n_s32(x) }
    }

    #[inline(always)]
    unsafe fn widen8to32(p: *const u8) -> Self::V32 {
        // Read exactly four bytes as a scalar `u32` rather than `vld1_u8`, which would read
        // eight. The delta buffers are padded by `LANES` (16) past `qlen`, and the group loop's
        // last group starts at `hi + 3 * LANES32`; an eight-byte read there would reach
        // `hi + 20`, which is past the padding. Four keeps it inside.
        let four = unsafe { std::ptr::read_unaligned(p as *const u32) };
        unsafe {
            // Broadcast puts the four bytes in both halves; only the low half is widened.
            let bytes = vreinterpret_u8_u32(vdup_n_u32(four));
            vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(vmovl_u8(bytes))))
        }
    }

    #[inline(always)]
    unsafe fn add32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { vaddq_s32(a, b) }
    }

    #[inline(always)]
    unsafe fn sub32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { vsubq_s32(a, b) }
    }

    #[inline(always)]
    unsafe fn max32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { vmaxq_s32(a, b) }
    }

    #[inline(always)]
    unsafe fn eq32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { vreinterpretq_s32_u32(vceqq_s32(a, b)) }
    }

    #[inline(always)]
    unsafe fn gt32(a: Self::V32, b: Self::V32) -> Self::V32 {
        unsafe { vreinterpretq_s32_u32(vcgtq_s32(a, b)) }
    }

    #[inline(always)]
    unsafe fn movemask32(v: Self::V32) -> u32 {
        unsafe {
            let mask = vreinterpretq_u32_s32(vshrq_n_s32(v, 31));
            let bits = vandq_u32(mask, vld1q_u32(LANE_BITS32.as_ptr()));
            // Four weights summing to at most 15: one reduction, no halving needed.
            vaddvq_u32(bits)
        }
    }

    // NEON is baseline on aarch64, so this attribute enables nothing that was not already on.
    // It is here so the selection stays honest if this is ever built for an aarch64 target
    // without it, in which case `select_backend`'s `cfg` excludes this backend entirely.
    #[target_feature(enable = "neon")]
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
