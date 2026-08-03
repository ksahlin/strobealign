//! Exact pairwise sequence alignment with affine gap costs.
//!
//! The entry point is [`SimdAligner`]: build one with a [`Scores`] scheme (or
//! [`Scores::default()`]) and call an alignment method on it repeatedly. Construct
//! one per thread and reuse it; the scoring matrix and every scratch buffer live
//! inside it and are reused across calls.
//!
//! Five alignment types are named on a single axis, *which side is local*, i.e. free
//! to stop or start wherever it scores best:
//!
//! - [`SimdAligner::global_alignment`]: nothing is local; both sequences are aligned
//!   end to end.
//! - [`SimdAligner::local_reference_end_alignment`]: the query is still spanned in
//!   full, but the reference may stop anywhere and its trailing tail is free.
//! - [`SimdAligner::local_reference_start_alignment`]: the mirror. The reference may
//!   *begin* anywhere and its leading prefix is free.
//! - [`SimdAligner::local_end_alignment`]: both sequences start at 0 and the alignment
//!   ends wherever it scores best; both trailing tails are free, and
//!   [`Scores::end_bonus`] nudges it towards covering the whole query.
//! - [`SimdAligner::local_start_alignment`]: the mirror. The alignment *begins*
//!   wherever it scores best and both leading prefixes are free.
//!
//! And one on a second axis, *which side is split*, taking **two** references:
//!
//! - [`SimdAligner::split_reference_alignment`]: one query aligned across a
//!   `left_reference` and a `right_reference` with a single jump between them,
//!   returning two [`AlignmentResult`]s that partition the query exactly. Built for
//!   spotting a tandem duplication that a linear aligner reported as an insertion.
//!
//! ```
//! use strobealign::simdaligner::{SimdAligner, Scores};
//! let mut aligner = SimdAligner::new(Scores::default());
//! let result = aligner.global_alignment(b"ACGTACGT", b"ACGTTACGT", None);
//! assert_eq!(result.cigar.to_string(), "3=1D5=");
//! ```
//!
//! Every alignment method takes a trailing `bandwidth: Option<usize>`: `None` aligns
//! exactly, `Some(w)` confines the fill to `w` diagonals either side of the anchor
//! diagonal. Each method's own `# bandwidth` section says where its band is anchored and
//! whether `w` can be widened.
//!
//! # Alphabet
//! The alphabet is **A, C, G, T and U** (case-insensitive), and nothing else. `U` is RNA's
//! T-equivalent and scores identically to `T`.
//!
//! Any other byte (`N`, an IUPAC ambiguity code, a `-` from a gapped FASTA, a stray `\n`) is
//! **unknown**, and an unknown byte never scores as a match against anything, *including
//! another unknown byte*.
//!
//! `=`/`X` in the CIGAR are **literal byte identity**, deliberately decoupled from how a
//! position scored. So `U` against `T` scores as a match but reports `X`, and an unknown byte
//! against a byte-identical copy of itself scores as a mismatch but reports `=`. Scoring is
//! about biology; `=`/`X` is about the bytes the caller handed us.
//!
//! # Guarantees
//! With `bandwidth: None` the alignment is *exact* at every length: the reported score is the
//! optimal score under the given scheme, and the CIGAR is an alignment achieving it. There is
//! no x-drop and no early exit. With `Some(w)` it is optimal among alignments staying inside
//! the band. Either way every input pair produces an alignment, including ones whose optimal
//! score is negative.
//!
//! `piecewisealigner.rs` passes `None` between anchors and `Some(w)` for the end extensions, so
//! "exact" is not a property of the aligner as strobealign uses it.
//!
//! # The kernel
//! One kernel (u8 lanes on a difference recurrence) and every call takes it. Written once,
//! generic over [`backend::Backend`], and instantiated per instruction set: AVX2 at 32 lanes,
//! SSE4.1 / NEON / wasm simd128 at 16, and a portable [`backend::Scalar`] fallback at 16 that
//! uses arrays instead of registers. [`SimdAligner::new`] picks the fastest one this machine can
//! run and **never panics for want of an instruction set**; [`SimdAligner::is_supported`] reports
//! whether that pick was a vectorized one.
//!
//! The **scoring scheme** is bounded, though (`match + 3*gap_open - gap_extend <= 255`, i.e.
//! `gap_open <= ~84` at `match = 2`), and [`SimdAligner::new`] rejects anything past it rather
//! than silently wrap a lane. Sequence *lengths* are not: every stored value is bounded by the
//! scheme alone, with no length term.

mod aligner;
pub mod backend;
mod kernel;

#[cfg(test)]
mod tests;

pub use aligner::SimdAligner;

// Re-exported so the kernel can reach them as `super::Cigar` without a conversion at the
// boundary.
pub use super::{
    aligner::Scores,
    cigar::{Cigar, CigarOperation},
};

/// Kernel internals, for tests. **Not part of the alignment API.**
#[doc(hidden)]
pub mod internal {
    pub use super::backend::{Backend, Scalar};
    pub use super::kernel::{U8Probe, fits_u8_lanes};

    #[cfg(all(
        target_arch = "aarch64",
        target_feature = "neon",
        target_endian = "little"
    ))]
    pub use super::backend::Neon;
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    pub use super::backend::Simd128;
    #[cfg(target_arch = "x86_64")]
    pub use super::backend::{Avx2, Sse41};

    /// The lane count of the widest kernel this build has. Anything backend-aware should read
    /// `B::LANES` instead; this stays 32 on x86_64 so that an out-of-tree oracle written against
    /// the AVX2 kernel does not silently start testing a different width.
    #[cfg(target_arch = "x86_64")]
    pub const U8_LANES: usize = <Avx2 as Backend>::LANES;
    #[cfg(not(target_arch = "x86_64"))]
    pub const U8_LANES: usize = 16;
}

/// Result of aligning a query against a reference.
///
/// `query_start`/`query_end`/`ref_start`/`ref_end` are half-open ranges into the
/// original slices, and the `cigar` covers exactly those ranges and nothing else: a
/// free end gap is **not** represented in the CIGAR.
///
/// Which of the four coordinates can move depends on the alignment type:
///
/// - [`SimdAligner::global_alignment`]: both spans cover their inputs in full.
/// - [`SimdAligner::local_reference_end_alignment`]: the query span is still the whole
///   query, but `ref_end` is wherever the alignment chose to stop.
/// - [`SimdAligner::local_reference_start_alignment`]: likewise, but it is `ref_start`
///   that moves and `ref_end` is always `reference.len()`.
/// - [`SimdAligner::local_end_alignment`]: `query_end` and `ref_end` both move.
/// - [`SimdAligner::local_start_alignment`]: `query_start` and `ref_start` both move.
///
/// `Default` is the **empty alignment**: score 0, all four coordinates at 0, and an empty
/// CIGAR. It is a real result, returned by [`SimdAligner::local_end_alignment`] when no
/// extension scores better than stopping immediately.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignmentResult {
    pub score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub cigar: Cigar,
}

/// Result of [`SimdAligner::split_reference_alignment`]: one query aligned across two
/// references, with a single jump point between them.
///
/// The two arms **partition the query exactly**: `right.query_start` is 0,
/// `left.query_end` is `query.len()`, and `right.query_end == left.query_start` (the
/// jump point). Each arm's CIGAR covers its own half of the query against its own
/// reference span, and the two CIGARs together consume every query base once.
///
/// `right` is anchored at the start of `right_reference` (`ref_start == 0`) and may
/// stop anywhere; `left` is anchored at the end of `left_reference`
/// (`ref_end == left_reference.len()`) and may begin anywhere.
///
/// `score` is `left.score + right.score`, and unlike the local types it **can be
/// negative**: the query must be covered whether or not either reference explains it.
/// An arm may be the empty alignment, or may consume no reference at all and report its
/// query segment as one pure insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitReferenceAlignment {
    /// The arm against `left_reference`: covers the query *suffix*, and ends at the
    /// end of the left reference.
    pub left: AlignmentResult,
    /// The arm against `right_reference`: covers the query *prefix*, and starts at the
    /// start of the right reference.
    pub right: AlignmentResult,
    /// `left.score + right.score`.
    pub score: i32,
}
