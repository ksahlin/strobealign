//! The SIMD operations the kernel needs, and one implementation per instruction set.
//!
//! The kernel in `kernel.rs` is written once, generic over [`Backend`]. A backend supplies two
//! vector types and 21 operations; everything else (the difference recurrence, the band geometry,
//! the traceback and the replay) is shared and arch-independent.
//!
//! # What a backend has to provide
//!
//! Two lane widths, because the kernel runs two things at once:
//!
//! - **`V8`, `LANES` lanes of `u8`**: the DP itself. This is the hot loop.
//! - **`V32`, `LANES32` lanes of `i32`**: the absolute-score reconstruction that rides
//!   alongside it for the modes that argmax over absolute `H[i][j]`. See `TRACK_ROW_MAX`.
//!
//! [`Backend::LANES`] must be exactly `4 * LANES32`, which every real instruction set satisfies
//! (AVX2: 32 and 8; SSE4.1/NEON/simd128: 16 and 4). The kernel's `nm` packing writes four
//! `LANES32`-lane groups into one `u32` per chunk and relies on it; [`AssertLaneRatio`] pins
//! it at compile time.
//!
//! # The two contracts that are easy to get wrong
//!
//! **Comparisons return an all-ones/all-zeros lane mask** (`0xFF` per byte for [`Backend::eq8`],
//! `0xFFFF_FFFF` per lane for [`Backend::eq32`]/[`Backend::gt32`]), not `1`. The kernel feeds
//! these to [`Backend::and8`]/[`Backend::or8`] and to the movemasks, both of which assume a
//! saturated mask.
//!
//! **The movemasks put lane `i` in bit `i`**, and take each lane's *most significant* bit.
//! `arg_for_row` and `replay` decode these bits by lane index, so a backend that packed them in
//! the other order would produce a silently wrong CIGAR under a right score.
//!
//! # `LANES` is not fixed at 32
//!
//! The traceback words stay `u32` on every backend: a 16-lane backend leaves the upper 16 bits of
//! each word unused rather than repacking. That wastes half of `tb`/`nm` on 128-bit targets and is
//! deliberate: it keeps one flag format across every arch, so `replay` and `arg_for_row` are
//! genuinely shared code rather than four near-copies.

use std::marker::PhantomData;

use super::kernel::U8Probe;

#[cfg(test)]
mod tests;

pub mod scalar;
pub use scalar::Scalar;

#[cfg(target_arch = "x86_64")]
pub mod x86;
#[cfg(target_arch = "x86_64")]
pub use x86::{Avx2, Sse41};

#[cfg(all(
    target_arch = "aarch64",
    target_feature = "neon",
    target_endian = "little"
))]
pub mod neon;
#[cfg(all(
    target_arch = "aarch64",
    target_feature = "neon",
    target_endian = "little"
))]
pub use neon::Neon;

#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
pub mod wasm;
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
pub use wasm::Simd128;

/// Compile-time check that a backend's two lane counts are in the ratio the kernel's `nm`
/// packing assumes. Instantiate it once per backend.
///
/// Being a post-monomorphization const error, it surfaces at codegen: `cargo build` and
/// `cargo test` report it, `cargo check` does not. The runtime assertions in `U8Probe::new` and
/// in the backend conformance tests are the backstop.
pub struct AssertLaneRatio<B: Backend>(PhantomData<B>);

impl<B: Backend> AssertLaneRatio<B> {
    pub const CHECK: () = assert!(
        B::LANES == 4 * B::LANES32,
        "the kernel packs four LANES32-lane groups into one u32 per chunk, so LANES must be \
         exactly 4 * LANES32"
    );
}

/// One SIMD instruction set, as the kernel uses it.
///
/// Every operation is `unsafe`: the loads and stores genuinely are, and the arithmetic is legal
/// only once the CPU has been checked. The check is carried by [`Backend::fill`], which holds the
/// `#[target_feature]`, not by the operations themselves.
///
/// # Safety
///
/// The contract is the same for all 21 operations, so it is stated once here rather than
/// repeated on each: the CPU must support this backend's instruction set, and pointers must be
/// valid for `LANES` (or `LANES32`) elements. `fill` states its own, which differs.
///
/// Implementations must also uphold the lane-mask and movemask contracts described in the
/// [module docs](self). A backend that violates them compiles and runs, and produces wrong
/// alignments.
#[allow(clippy::missing_safety_doc)]
pub unsafe trait Backend: Copy + 'static {
    /// Human-readable name, for the dispatch log line.
    const NAME: &'static str;
    /// `u8` lanes per vector. Must equal `4 * LANES32`.
    const LANES: usize;
    /// `i32` lanes per vector.
    const LANES32: usize;

    /// Whether this CPU can actually execute this backend.
    ///
    /// Safe to call anywhere, and the single source of truth for selection. `U8Probe::new`
    /// asserts it, so that naming a kernel for an instruction set the CPU lacks panics rather
    /// than being undefined behaviour, which it otherwise would be, since the alignment methods
    /// are safe fns and `U8Probe` is re-exported through `simdaligner::internal`.
    fn is_available() -> bool;

    /// `LANES` lanes of `u8`.
    type V8: Copy;
    /// `LANES32` lanes of `i32`.
    type V32: Copy;

    // --- u8 lanes: the DP itself ---

    /// Unaligned load of `LANES` bytes. The kernel never aligns its buffers.
    unsafe fn load8(p: *const u8) -> Self::V8;
    /// Unaligned store of `LANES` bytes.
    ///
    /// Always writes all `LANES` lanes, including on the partial last chunk: the kernel's
    /// write-ordering hazards are built on that. Never narrow this to a masked store.
    unsafe fn store8(p: *mut u8, v: Self::V8);
    /// Broadcast one byte to every lane.
    unsafe fn splat8(x: u8) -> Self::V8;
    /// Lanewise **unsigned** max. The whole u8 kernel rests on this being unsigned.
    unsafe fn max8(a: Self::V8, b: Self::V8) -> Self::V8;
    /// Lanewise wrapping add. The kernel's values are provably in range, so wrapping vs
    /// saturating never differs; do not substitute a saturating add, which would change the
    /// answer if the range invariant were ever broken, and hide it.
    unsafe fn add8(a: Self::V8, b: Self::V8) -> Self::V8;
    /// Lanewise wrapping subtract. Same note as [`Backend::add8`].
    unsafe fn sub8(a: Self::V8, b: Self::V8) -> Self::V8;
    /// Lanewise equality, returning `0xFF` per equal byte and `0x00` otherwise.
    unsafe fn eq8(a: Self::V8, b: Self::V8) -> Self::V8;
    /// Bitwise and.
    unsafe fn and8(a: Self::V8, b: Self::V8) -> Self::V8;
    /// Bitwise or.
    unsafe fn or8(a: Self::V8, b: Self::V8) -> Self::V8;
    /// Select per lane: `if_true` where `mask` is set, `if_false` where it is clear.
    ///
    /// `mask` is always a saturated comparison mask at every call site, so a bit-select and an
    /// MSB-select are equivalent here.
    unsafe fn blend8(if_false: Self::V8, if_true: Self::V8, mask: Self::V8) -> Self::V8;
    /// The MSB of each lane, lane `i` in bit `i`. Bits above `LANES` are zero.
    unsafe fn movemask8(v: Self::V8) -> u32;

    // --- i32 lanes: the absolute-score reconstruction ---

    /// Unaligned load of `LANES32` `i32`s.
    unsafe fn load32(p: *const i32) -> Self::V32;
    /// Unaligned store of `LANES32` `i32`s.
    unsafe fn store32(p: *mut i32, v: Self::V32);
    /// Broadcast one `i32` to every lane.
    unsafe fn splat32(x: i32) -> Self::V32;
    /// Load `LANES32` **bytes** and zero-extend each to an `i32` lane.
    ///
    /// Zero-extension, not sign-extension: the deltas are unsigned and the kernel un-offsets
    /// them afterwards by subtracting `gap_open`.
    ///
    /// **Must not read more than `LANES32` bytes.** The kernel's last group sits flush against
    /// the end of the delta buffers, so a wider load (`vld1_u8`'s 8 bytes where 4 were asked
    /// for, say) reads out of bounds. The padding leaves no slack to spare.
    unsafe fn widen8to32(p: *const u8) -> Self::V32;
    /// Lanewise wrapping add.
    unsafe fn add32(a: Self::V32, b: Self::V32) -> Self::V32;
    /// Lanewise wrapping subtract.
    unsafe fn sub32(a: Self::V32, b: Self::V32) -> Self::V32;
    /// Lanewise **signed** max. These are absolute scores and are routinely negative.
    unsafe fn max32(a: Self::V32, b: Self::V32) -> Self::V32;
    /// Lanewise equality, returning all-ones per equal lane.
    unsafe fn eq32(a: Self::V32, b: Self::V32) -> Self::V32;
    /// Lanewise **signed** greater-than, returning all-ones per lane where `a > b`.
    unsafe fn gt32(a: Self::V32, b: Self::V32) -> Self::V32;
    /// The MSB of each lane, lane `i` in bit `i`. Bits above `LANES32` are zero.
    unsafe fn movemask32(v: Self::V32) -> u32;

    /// Run the kernel's fill with this backend's instruction set enabled.
    ///
    /// # Why this exists, and why it is not just a call to `U8Probe::fill_generic`
    ///
    /// On x86 the intrinsics need `#[target_feature(enable = "avx2")]`, and a
    /// `#[target_feature]` function **cannot be inlined into a caller that lacks the feature**.
    /// Putting the attribute on the 21 operations above would therefore leave 21 real function
    /// calls per cell. (The obvious fix, `#[inline(always)]` alongside `#[target_feature]`, is
    /// still unstable: see rust-lang/rust#145574. This crate builds on stable.)
    ///
    /// So the feature goes **here** instead, on a trampoline wrapping the entire fill. The
    /// operations stay plain `#[inline(always)]` functions with no feature of their own, so they
    /// inline freely into this one and end up inside a function that does have the feature. That
    /// is one `call` per alignment and none per cell, and it also puts the surrounding scalar
    /// work under the same target.
    ///
    /// # Safety
    ///
    /// The CPU must actually support this backend's instruction set. `select_backend` in
    /// `aligner.rs` is the only caller's source of `Self`, and it detects before it selects.
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
    ) -> (i32, usize)
    where
        Self: Sized;
}
