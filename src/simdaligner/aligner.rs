use super::Scores;

use super::backend::{self, Backend};
use super::kernel::{U8Probe, fits_u8_lanes};
use super::{AlignmentResult, SplitReferenceAlignment};

/// The selected kernel, one variant per instruction set this build can reach.
///
/// Which variants exist is decided by `cfg`; which one is built is decided at runtime only on
/// x86_64, where AVX2 is not baseline.
enum Workspace {
    #[cfg(target_arch = "x86_64")]
    Avx2(U8Probe<backend::Avx2>),
    #[cfg(target_arch = "x86_64")]
    Sse41(U8Probe<backend::Sse41>),
    #[cfg(all(
        target_arch = "aarch64",
        target_feature = "neon",
        target_endian = "little"
    ))]
    Neon(U8Probe<backend::Neon>),
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    Simd128(U8Probe<backend::Simd128>),
    Scalar(U8Probe<backend::Scalar>),
}

/// Pick the fastest kernel this machine can run.
///
/// Ordered fastest-first per architecture, with [`backend::Scalar`] as the last resort.
#[allow(unreachable_code)]
fn select_backend() -> Workspace {
    #[cfg(target_arch = "x86_64")]
    {
        if backend::Avx2::is_available() {
            return Workspace::Avx2(U8Probe::new());
        }
        if backend::Sse41::is_available() {
            return Workspace::Sse41(U8Probe::new());
        }
    }
    #[cfg(all(
        target_arch = "aarch64",
        target_feature = "neon",
        target_endian = "little"
    ))]
    {
        return Workspace::Neon(U8Probe::new());
    }
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        return Workspace::Simd128(U8Probe::new());
    }

    Workspace::Scalar(U8Probe::new())
}

/// Say so, once per process, when the aligner drops to the portable kernel.
///
/// Worth logging at all because nothing else would reveal it: the portable kernel gives
/// identical alignments, just several times slower, which is indistinguishable from a large
/// workload.
fn warn_scalar_fallback() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        log::warn!(
            "no SIMD kernel is available for this CPU ({}); falling back to the portable \
             scalar kernel, which gives identical alignments several times slower. On \
             x86_64 this means the CPU lacks both AVX2 and SSE4.1; on wasm32 it means the \
             build did not enable simd128 (`-C target-feature=+simd128`).",
            std::env::consts::ARCH
        );
    });
}

/// The two preconditions the kernel cannot check for itself, asserted once per scheme so the
/// alignment path carries no validation. Shared by every constructor.
fn validate(scores: Scores) {
    assert!(
        scores.gap_open >= scores.gap_extend,
        "gap_open ({}) must be >= gap_extend ({}): the kernel collapses the open-vs-extend \
         recurrence on that assumption",
        scores.gap_open,
        scores.gap_extend
    );
    assert!(
        fits_u8_lanes(scores.match_, scores.gap_open, scores.gap_extend),
        "scores do not fit the u8 kernel: match ({}) + 3*gap_open ({}) - gap_extend ({}) \
         must be <= 255, i.e. gap_open <= ~84 at match=2. Past this an intermediate wraps \
         and the score comes back better than optimal, silently.",
        scores.match_,
        scores.gap_open,
        scores.gap_extend
    );
}

/// Call the same method on whichever kernel was selected.
///
/// The body is written once and monomorphized per backend by inference on `$ws`, so dispatch
/// costs one well-predicted branch per alignment call, not per cell.
macro_rules! dispatch {
    ($self:ident, $ws:ident => $body:expr) => {
        match &mut $self.workspace {
            #[cfg(target_arch = "x86_64")]
            Workspace::Avx2($ws) => $body,
            #[cfg(target_arch = "x86_64")]
            Workspace::Sse41($ws) => $body,
            #[cfg(all(
                target_arch = "aarch64",
                target_feature = "neon",
                target_endian = "little"
            ))]
            Workspace::Neon($ws) => $body,
            #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
            Workspace::Simd128($ws) => $body,
            Workspace::Scalar($ws) => $body,
        }
    };
}

/// A reusable aligner configured with a fixed [`Scores`] scheme.
pub struct SimdAligner {
    workspace: Workspace,
    scores: Scores,
}

impl Default for SimdAligner {
    fn default() -> Self {
        SimdAligner::new(Scores::default())
    }
}

impl SimdAligner {
    /// Builds an aligner for the given scoring scheme.
    ///
    /// The scheme is validated and the backend selected here, once per scheme, so that the
    /// alignment path carries no feature detection and no validation.
    ///
    /// # Panics
    ///
    /// - If `scores.gap_open < scores.gap_extend`. The kernel collapses the gap-open and
    ///   gap-extend cases into a single recurrence on that assumption.
    /// - If `match_ + 3 * gap_open - gap_extend > 255`, i.e. roughly `gap_open > 84` at
    ///   `match_ = 2`. Beyond that bound an intermediate value exceeds the kernel's 8-bit
    ///   lanes.
    pub fn new(scores: Scores) -> Self {
        validate(scores);
        let workspace = select_backend();
        if matches!(workspace, Workspace::Scalar(_)) {
            warn_scalar_fallback();
        }
        SimdAligner { workspace, scores }
    }

    /// An aligner pinned to the portable kernel, bypassing backend selection. **Tests only.**
    #[cfg(test)]
    pub(crate) fn with_scalar_kernel(scores: Scores) -> Self {
        validate(scores);
        SimdAligner {
            workspace: Workspace::Scalar(U8Probe::new()),
            scores,
        }
    }

    /// Whether this build has a **vectorized** kernel for this machine.
    ///
    /// A performance question, not an availability one: `false` means [`SimdAligner::new`] falls
    /// back to the portable kernel, which is correct and slow. It never panics for want of an
    /// instruction set.
    pub fn is_supported() -> bool {
        !matches!(select_backend(), Workspace::Scalar(_))
    }

    /// The name of the kernel this aligner is using, for diagnostics: `"AVX2"`, `"SSE4.1"`,
    /// `"NEON"`, `"simd128"` or `"scalar"`.
    pub fn backend_name(&self) -> &'static str {
        match &self.workspace {
            #[cfg(target_arch = "x86_64")]
            Workspace::Avx2(_) => backend::Avx2::NAME,
            #[cfg(target_arch = "x86_64")]
            Workspace::Sse41(_) => backend::Sse41::NAME,
            #[cfg(all(
                target_arch = "aarch64",
                target_feature = "neon",
                target_endian = "little"
            ))]
            Workspace::Neon(_) => backend::Neon::NAME,
            #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
            Workspace::Simd128(_) => backend::Simd128::NAME,
            Workspace::Scalar(_) => backend::Scalar::NAME,
        }
    }

    /// The scoring scheme this aligner was built with.
    pub fn scores(&self) -> Scores {
        self.scores
    }

    /// Global (end-to-end) alignment: both sequences are aligned in full, and terminal gaps are
    /// penalised.
    ///
    /// `query_start`/`query_end` and `ref_start`/`ref_end` always span the entire inputs.
    ///
    /// Empty inputs follow affine-gap semantics: aligning a length-`n` sequence against an empty
    /// one yields a single length-`n` gap scoring `-(gap_open + (n - 1) * gap_extend)`, and two
    /// empty inputs yield an empty alignment scoring 0.
    ///
    /// # `bandwidth`
    ///
    /// Global is the one mode pinned at **both** corners, and `(qlen, rlen)` sits `|qlen - rlen|`
    /// off the diagonal a band is centred on, so a narrower band contains no path at all and `w`
    /// is widened to reach it.
    pub fn global_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        dispatch!(self, ws => ws.global_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        ))
    }

    /// Alignment with a **local reference end**: the query is spanned in full, but the reference
    /// may end anywhere, and the unaligned reference tail costs nothing.
    ///
    /// Formally, the optimum over every end cell `(query.len(), j)` for `j` in
    /// `0..=reference.len()`. `query_start`/`query_end` span the whole query, `ref_start` is 0,
    /// and `ref_end` is where the alignment chose to stop: the CIGAR covers exactly
    /// `reference[..ref_end]` and says nothing about the rest. `j = 0` is a legal end cell, in
    /// which case the whole reference is the free end gap. Among equally-scoring end cells the
    /// *smallest* `ref_end` wins.
    ///
    /// # `bandwidth`
    ///
    /// Anchored at `(0, 0)`, like [`local_end_alignment`](SimdAligner::local_end_alignment), so
    /// `Some(w)` bounds `ref_end` to `[qlen - w, qlen + w]`. Spanning the query needs
    /// `w >= qlen - rlen`, and `w` is widened to that where necessary.
    pub fn local_reference_end_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        dispatch!(self, ws => ws.local_reference_end_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        ))
    }

    /// Alignment with a **local reference start**: the mirror of
    /// [`local_reference_end_alignment`](SimdAligner::local_reference_end_alignment). The query
    /// is spanned in full, but the reference may *begin* anywhere, and the unaligned reference
    /// prefix costs nothing.
    ///
    /// Formally, the optimum over every start cell `(0, s)` for `s` in `0..=reference.len()`,
    /// ending at `(query.len(), reference.len())`. `query_start`/`query_end` span the whole
    /// query, `ref_end` is always `reference.len()`, and `ref_start` is where the alignment
    /// chose to begin: the CIGAR covers exactly `reference[ref_start..]`. Among equally-scoring
    /// start cells the *largest* `ref_start` wins: the shortest reference span, mirroring the
    /// `_end` variant's rule.
    ///
    /// # `bandwidth`
    ///
    /// As [`local_reference_end_alignment`](SimdAligner::local_reference_end_alignment),
    /// **mirrored with the mode**: the band is anchored at the *ends* of the two inputs, where
    /// this alignment is pinned, and bounds `ref_start` to within `w` of `rlen - qlen`.
    pub fn local_reference_start_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        dispatch!(self, ws => ws.local_reference_start_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        ))
    }

    /// **Local-end alignment**: both sequences start at 0, and the alignment ends wherever it
    /// scores best. Both trailing tails (query and reference) are free.
    ///
    /// Formally, the optimum over every end cell `(i, j)` of
    /// `H[i][j] + (if i == query.len() { scores.end_bonus } else { 0 })`. `query_end` and
    /// `ref_end` are wherever it stopped, and the CIGAR covers exactly `query[..query_end]`
    /// against `reference[..ref_end]`. Ties go to the longest extension: the largest
    /// `query_end` wins, then the largest `ref_end`.
    ///
    /// `end_bonus` is added only if the alignment covers the whole query, and never enters the
    /// recurrence; it is applied once, at the end, to choose between the best cell anywhere and
    /// the best cell that finishes the query. That nudge is the difference between this and
    /// [`local_reference_end_alignment`](SimdAligner::local_reference_end_alignment).
    ///
    /// The CIGAR ends with an insertion only if `end_bonus >= gap_open`, and never with a
    /// deletion when `gap_open > 0`: a local maximum cannot end in a gap. `ref_end` does not
    /// advance with a trailing insertion.
    ///
    /// # The score is never negative
    ///
    /// Cell `(0, 0)` is a legal end cell, so a hopeless extension returns the empty alignment:
    /// spans `0..0` on both sequences, score 0, empty CIGAR. That is a normal result.
    ///
    /// One exception: an **empty query** trivially covers itself, so `(0, 0)` earns `end_bonus`
    /// and the empty alignment comes back scoring `end_bonus` rather than 0. Don't use
    /// `score == 0` alone as a "gave up" sentinel unless the query is known non-empty.
    ///
    /// # `bandwidth`
    ///
    /// Anchored at `(0, 0)`, which is in every band, so `Some(w)` is exactly `w`: this and
    /// [`local_start_alignment`](SimdAligner::local_start_alignment) are the only two modes
    /// never widened.
    pub fn local_end_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        dispatch!(self, ws => ws.local_end_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            self.scores.end_bonus,
            bandwidth,
        ))
    }

    /// **Local-start alignment**: the mirror of
    /// [`local_end_alignment`](SimdAligner::local_end_alignment). Both sequences end at their
    /// ends, and the alignment *begins* wherever it scores best; both leading prefixes are free.
    ///
    /// Formally, the optimum over every start cell `(a, b)` of the best alignment of
    /// `query[a..]` against `reference[b..]`, plus `end_bonus` if `a == 0`. `query_start` and
    /// `ref_start` are wherever it began; `query_end` and `ref_end` are always the ends of the
    /// inputs. Ties go to the longest extension, mirrored: the *smallest* `query_start` wins,
    /// then the smallest `ref_start`. Likewise a *leading* insertion requires
    /// `end_bonus >= gap_open`, and with `gap_open > 0` the CIGAR never begins with a deletion.
    ///
    /// The score is never negative here either (the start cell
    /// `(query.len(), reference.len())` is legal and scores 0, at the *end* of both sequences
    /// rather than the start), and the empty-query exception applies:
    /// see [`local_end_alignment`](SimdAligner::local_end_alignment).
    ///
    /// # `bandwidth`
    ///
    /// As [`local_end_alignment`](SimdAligner::local_end_alignment), **mirrored with the mode**:
    /// the band is anchored at the *ends* of the two inputs, where this alignment is pinned, so
    /// it constrains start cells to `|(qlen - a) - (rlen - b)| <= w`.
    pub fn local_start_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        dispatch!(self, ws => ws.local_start_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            self.scores.end_bonus,
            bandwidth,
        ))
    }

    /// **Split-reference alignment**: one query aligned across two references, with a single
    /// jump between them. Built for detecting a tandem duplication that a linear aligner reports
    /// as an insertion.
    ///
    /// The alignment starts at `query[0]` against `right_reference[0]` and ends at the end of
    /// `left_reference`, jumping once from the right reference to the left. It returns two
    /// [`AlignmentResult`]s that **partition the query exactly**:
    ///
    /// - `right` covers `query[..k]` against a prefix of `right_reference`;
    /// - `left` covers `query[k..]` against a suffix of `left_reference`;
    ///
    /// where `k`, the jump point, is chosen to maximise the total score:
    ///
    /// ```text
    /// score = max over k in 0..=query.len() of  F[k] + G[k]
    ///   F[k] = the best alignment of query[..k] against any prefix of right_reference
    ///   G[k] = the best alignment of query[k..] against any suffix of left_reference
    /// ```
    ///
    /// `end_bonus` is ignored: the query is spanned in full by construction.
    ///
    /// Total score decides between splits. Among equal totals the **most balanced** wins (the
    /// smallest `|left.score - right.score|`, since `5 + 5` describes a real duplication and
    /// `-20 + 30` does not), then the smallest `k`. Each arm's own reference endpoint breaks
    /// ties the way its one-reference counterpart does.
    ///
    /// # The score can be negative, and an arm can be empty
    ///
    /// Unlike the local modes there is no "give up" option: the two arms must cover the query,
    /// so a query that neither reference explains still gets covered, at a cost. `k = 0` and
    /// `k = query.len()` are legal, leaving one arm spanning the whole query and the other the
    /// empty alignment; and an arm may consume no reference at all, reporting its query segment
    /// as one pure insertion with `ref_start == ref_end`.
    ///
    /// # `bandwidth`
    ///
    /// `Some(w)` is **one band, applied to both arms**, each measured from that arm's own anchor:
    /// the right arm from `(query[0], right_reference[0])`, the left arm from the ends of `query`
    /// and `left_reference`. **It does not constrain the jump point**: `k` is still chosen over
    /// every value the two bands admit.
    ///
    /// Both arms must together cover the query, and an arm reaches at most `reflen + w` query
    /// bases in band, so `w` is widened where necessary to
    /// `(qlen - left_reference.len() - right_reference.len()) / 2`.
    pub fn split_reference_alignment(
        &mut self,
        query: &[u8],
        left_reference: &[u8],
        right_reference: &[u8],
        bandwidth: Option<usize>,
    ) -> SplitReferenceAlignment {
        dispatch!(self, ws => ws.split_reference_alignment(
            query,
            left_reference,
            right_reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        ))
    }
}
