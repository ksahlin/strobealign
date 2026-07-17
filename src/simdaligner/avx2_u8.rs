//! AVX2 `u8` kernel: the difference recurrence on anti-diagonals, 32 lanes.
//!
//! This is the only kernel; every alignment type routes at it, differing only in which cell the
//! traceback starts from and which maxima the fill tracks.
//!
//! # Where things are
//!
//! A public mode on `U8Probe` picks a band strategy, then goes through one of three wrappers to
//! the fill and the traceback:
//!
//! ```text
//! pub fn <mode>  ->  align | local | row_pass | split_reference
//!                        -> prepare     encode both sequences, size the delta buffers
//!                        -> size_flags  size the tb/nm buffers
//!                        -> run         the fill: deltas forward, flags out
//!                        -> arg_for_row[_banded]   which column won a row (from the nm bits)
//!                        -> replay      walk the flags back, emit the CIGAR
//! ```
//!
//! See `run`'s "Which mode compiles which flags" table for what each mode switches on.
//!
//! # The difference recurrence
//!
//! Rather than storing the DP cells, this stores the **differences between adjacent cells**
//! (Suzuki & Kasahara 2018, BMC Bioinformatics 19(Suppl 1):45). The differences are bounded by
//! the **scoring scheme alone**, with no length term, so they fit `u8` lanes at any sequence
//! length - which is what makes an exact, unbanded aligner affordable, and what removes any need
//! for an overflow envelope or a wider fallback kernel.
//!
//! The recurrence is the paper's Eq 3, the *offset* form, in this crate's gap convention
//! (`go = gap_open - gap_extend`, `ge = gap_extend`, `C = 2*gap_open`, `s_G(x,y) = s(x,y) + C`).
//! The offset variables fold the constants in so that `ΔE'_G` absorbs `ΔV` and `ΔF'_G` absorbs
//! `ΔH`, which leaves nothing to add before the maxes and gives a critical path of 4:
//!
//! ```text
//! A_G  [i,j] = max( s_G(a,b), ΔE'_G[i-1,j], ΔF'_G[i,j-1] )
//! ΔH_G [i,j] = A_G[i,j] - ΔV_G[i-1,j]
//! ΔV_G [i,j] = A_G[i,j] - ΔH_G[i,j-1]
//! ΔE'_G[i,j] = max( A_G[i,j], ΔE'_G[i-1,j] + go ) - ΔH_G[i,j-1]
//! ΔF'_G[i,j] = max( A_G[i,j], ΔF'_G[i,j-1] + go ) - ΔV_G[i-1,j]
//! ```
//!
//! with the boundary `ΔV_G[0,j] = ΔE'_G[0,j] = ΔH_G[i,0] = ΔF'_G[i,0] = 0` at index 1 and `go`
//! after. Note the paper's Eq 2 is the same idea in a *naive* difference form; its critical path
//! is 8, so Eq 3 is the one to build on. (Two transcription traps, if comparing against the
//! paper: its gap convention is `G(k) = G_o + k*G_e` against this crate's
//! `gap_open + (k-1)*gap_extend`, so `G_o = gap_open - gap_extend`; and its supplement's Eq 19
//! is sign-flipped against its own Eq 7.)
//!
//! Everything is unsigned: `0 <= ΔH_G, ΔV_G, ΔE'_G, ΔF'_G <= match + 2*gap_open`. The hard
//! **zero** lower bound is the load-bearing half - it is what permits `u8` lanes with
//! `max_epu8`, and it means the kernel needs **no sentinel at all**. Every value is provably in
//! `[0, upper]`, so the subtractions cannot wrap and plain `_mm256_sub_epi8` suffices.
//!
//! # Why anti-diagonals, and why there are no shifts
//!
//! `ΔE'_G[i,j]` depends on `ΔE'_G[i-1,j]` - same column, previous row - so a column-vectorised
//! layout has an intra-vector serial dependency. On an anti-diagonal `p = i + j`, `A_G`'s inputs
//! are:
//!
//! ```text
//! (i,j) needs  ΔV_G[i-1, j]  and ΔE'_G[i-1, j]     <- anti-diagonal p-1, row i-1
//!              ΔH_G[i, j-1]  and ΔF'_G[i, j-1]     <- anti-diagonal p-1, row i
//! ```
//!
//! **Everything comes from `p-1`, and there is no `p-2` term**: `A_G` *is* the diagonal step, so
//! the diagonal predecessor is already folded in. Holding each anti-diagonal in a buffer indexed
//! by **row** makes "row `i-1`" and "row `i`" two `loadu`s one byte apart, so the shifts the
//! paper needs for its fixed-width sliding vector are not needed here - the addressing does it.
//!
//! The sequences fall out the same way. Cell `(i, p-i)` reads `query[i-1]`, contiguous in `i`,
//! and `refseq[p-i-1]`, which runs *backwards* in `i` - so the reference is stored **reversed**
//! and becomes `ref_rev[rlen + i - p]`, contiguous too. (Group it `rlen + i - p`, never
//! `rlen - p + i`: the latter underflows a `usize` once the anti-diagonal starts sliding right.)
//!
//! # Two hazards
//!
//! - **Boundary cells are written after the chunk loop, never before.** A 32-lane store is
//!   unconditional, so the partial last chunk scribbles over the row just past the interior -
//!   which for `p <= qlen` is the `(p, 0)` boundary cell.
//! - **Garbage lanes need no mask.** Anti-diagonal cells have no cross-lane dependency, so
//!   garbage stays in its own lane; and anti-diagonal `p` provably reads only rows
//!   `[lo(p-1), hi(p-1)]`, which `p-1`'s interior and boundary writes cover exactly.

use std::arch::x86_64::*;

use super::{AlignmentResult, Cigar, CigarOperation, SplitReferenceAlignment};

/// Lanes per vector: 32 `u8`s in a `__m256i`.
pub const LANES: usize = 32;

/// Traceback flag words per 32-lane chunk: `from_diag`, `from_e`, `e_open`, `f_open`.
///
/// **The alignment state is already in the deltas** - this is the observation the whole
/// traceback rests on, and it is why the difference form is *easier* to trace back than the
/// absolute one, not harder.
///
/// `A_G = max(s_G, ΔE'_G[i-1,j], ΔF'_G[i,j-1])`, so *which operand achieved the max* is the
/// traceback direction, and the fill already computed it. Likewise
/// `ΔE'_G[i,j] = max(A_G, ΔE'_G[i-1,j] + go) - ΔH_G[i,j-1]`: which operand won says whether the
/// gap opened or extended - our open/extend bits, for free.
///
/// So every flag is a `_mm256_cmpeq_epi8` against values already in registers, extracted with
/// `_mm256_movemask_epi8`: **4 bits per cell**, 8 ops per 32 cells.
///
/// Every flag is an *equality*, deliberately. AVX2 has no unsigned byte compare
/// (`_mm256_cmpgt_epi8` is signed, and these values can exceed 127 with a large `gap_open`), so
/// "did `ΔE'_G[i-1,j] + go` win the max" is phrased as `max(A_G, ΔE'_G+go) == A_G` - reusing the
/// max the fill already needs - rather than as a `>`. Equality does not care about signedness.
///
/// # Only three of the four directions are stored
///
/// `A_G`'s max has three operands and at least one of them achieves it, so one direction can be
/// left unstored and recovered as the fallthrough. **Which one is not free to choose: the
/// inferred direction is the one the replay checks last, and therefore the one it can never
/// prefer on a tie.** Since the replay prefers the diagonal (see `replay`), `from_diag` is
/// stored and `from_f` is inferred. (`e_open`/`f_open` are a two-way max, so neither is
/// inferable from the other; both stay.)
///
/// `from_diag` cannot be `cmpeq(A_G, s_G)` on its own - see [`TB_FROM_DIAG`].
const TB_WORDS: usize = 4;

/// One `u32` of "this cell took its row's max" per 32-lane chunk - a bit per lane, in the same
/// `(p, chunk)` coordinates as the traceback flags, so `tb_base[p]` indexes both. See
/// `U8Probe::nm`.
const NM_WORDS: usize = 1;

const TB_FROM_E: usize = 0;
const TB_E_OPEN: usize = 1;
const TB_F_OPEN: usize = 2;
/// Did the **diagonal** achieve `A_G`'s max?
///
/// `cmpeq(A_G, s_G)` is the obvious answer and it is **wrong for a clamped `s_G`**. Any scheme
/// with `mismatch > 2*gap_open` makes `s_G(mismatch) = C - mismatch` negative, and a u8 lane
/// cannot hold that, so the fill clamps it to zero. Clamping is lossless *for the value* -
/// `A_G` only ever maxes `s_G` against `ΔE'_G`/`ΔF'_G`, which are `>= 0`, so a negative
/// diagonal loses every max it enters and zero loses just as thoroughly. It is **not** lossless
/// for the flag: `A_G == s_G == 0` then holds while the true diagonal is strictly negative and
/// the real max came from a gap. Believing that bit takes a path worth `mismatch` points less
/// than the fill reported - a wrong CIGAR under a right score.
///
/// The fix costs two ops and no extra state: on a clamped scheme the diagonal is masked off at
/// **mismatch cells**, which is exactly where the clamp bites, using the `q == r` mask the
/// blend already computed. It is sound rather than conservative: where `s_G` clamps, the true
/// diagonal is negative while `ΔE'_G`/`ΔF'_G` are `>= 0`, so at a mismatch cell the diagonal
/// **provably never achieves the max** and masking it off discards nothing.
///
/// Only a test that re-scores the CIGAR can catch a regression here; one that trusts the
/// reported score cannot, since the score stays right while the CIGAR goes wrong.
const TB_FROM_DIAG: usize = 3;

/// The interior row range of anti-diagonal `p`, band included.
///
/// # Band geometry on the anti-diagonal
///
/// A band is the constraint `|i - j| <= w`. On anti-diagonal `p` the cell in row `i` sits at
/// column `j = p - i`, so `i - j = 2i - p` and the constraint reads
///
/// ```text
/// |2i - p| <= w   <=>   ceil((p - w) / 2) <= i <= floor((p + w) / 2)
/// ```
///
/// That is **two more clamps on a range the fill already clamps twice**. The matrix corners
/// already force `i >= max(1, p - rlen)` and `i <= min(qlen, p - 1)`; the band is `max` and
/// `min` against one more term each, and every loop below inherits it for free. Nothing else
/// about the fill changes shape: same recurrence, same lanes, same stores, fewer of them.
///
/// The width is `w + 1` rows when `p + w` is even and `w` when it is odd - the two bounds
/// advance on alternating anti-diagonals, since each moves at half of `p`'s rate. So the fill
/// costs `(qlen + rlen) * (w + 1) / LANES` chunks instead of `qlen * rlen / LANES`, which at
/// the shape this is for (`q = 20000`, `r = 24000`, `w = 100`) is 4.4M cells against 480M.
///
/// Returns `(lo, hi)`, inclusive. `lo > hi` means the anti-diagonal has no interior, which
/// happens at `p = 1` and nowhere else for `w >= 1`.
#[inline(always)]
fn band_rows<const BANDED: bool>(p: usize, qlen: usize, rlen: usize, w: usize) -> (usize, usize) {
    let lo = 1.max(p.saturating_sub(rlen));
    let hi = qlen.min(p - 1);
    if BANDED {
        // `saturating_sub` on both: `p - w` goes negative for every anti-diagonal before the
        // band's lower edge lifts off row 0, and `p - rlen` for every one before the band
        // starts sliding right. Same trap, same fix, and the same one the unbanded bound
        // already needed.
        (lo.max(p.saturating_sub(w).div_ceil(2)), hi.min((p + w) / 2))
    } else {
        (lo, hi)
    }
}

/// Does a band of `w` reach every cell of a `qlen x rlen` matrix?
///
/// The furthest any cell sits from the main diagonal is `max(qlen, rlen)` - the corners
/// `(qlen, 0)` and `(0, rlen)`. At or past that the band excludes nothing, so the banded optimum
/// *is* the exact optimum: same candidate set, same tie-breaks, same CIGAR.
///
/// # Why the callers take the exact path instead
///
/// Not to change the answer - it cannot - but because the exact path is **cheaper**, and most so
/// exactly where this fires. `global` and the `local_reference_*` modes get their exact score
/// from the cheap `anchor` walk, while every banded path runs `TRACK_ROW_MAX`; and the band's own
/// per-anti-diagonal bookkeeping is pure waste when the strip excludes nothing.
///
/// This is constant folding, not a policy: `Some(w)` still means exactly what it says, and a
/// band wider than the problem simply stops costing anything.
#[inline]
fn band_reaches_everything(qlen: usize, rlen: usize, w: usize) -> bool {
    w >= qlen.max(rlen)
}

/// The query counterpart of [`ref_span`]. Only the `local_*` modes - the ones that let the
/// query end early - ever pass `consumed != qlen`.
#[inline]
fn query_span<const REVERSED: bool>(qlen: usize, consumed: usize) -> (usize, usize) {
    if REVERSED {
        (qlen - consumed, qlen)
    } else {
        (0, consumed)
    }
}

/// Map "how much reference the CIGAR consumed" onto a span of the caller's reference.
///
/// Under `REVERSED` the kernel ran on back-to-front codes, so a run of `consumed` bases
/// starting at the reversed origin is a run *ending* at the forward end. One subtraction - the
/// only coordinate work `REVERSED` costs.
#[inline]
fn ref_span<const REVERSED: bool>(rlen: usize, consumed: usize) -> (usize, usize) {
    if REVERSED {
        (rlen - consumed, rlen)
    } else {
        (0, consumed)
    }
}

/// Does this scoring scheme survive `u8` lanes?
///
/// The bound is `M + 3*gap_open - gap_extend`, **not** the `M + 2*gap_open` that bounds the
/// stored values. The looser rule is a trap: `ΔE'_G`'s recurrence forms `ΔE'_G[i-1,j] + go`
/// **before** subtracting, and that intermediate is not a stored value, so it is not covered.
/// It reaches `upper + go`.
///
/// Saturating instead of rejecting does not help: `_mm256_adds_epu8` clamping to 255 makes the
/// following `max` pick the wrong operand - a wrong answer rather than a caught one.
///
/// At default scores this is `2 + 36 - 1 = 37` against a ceiling of 255, and it holds out to
/// `gap_open ~ 84`, so nothing real is near it.
///
/// Note this takes **no lengths**: the bound depends on the scoring scheme alone, so it is
/// decided once per scheme rather than per call, and no sequence is ever too long for it.
pub fn fits_u8_lanes(match_score: u8, gap_open: u8, gap_extend: u8) -> bool {
    let upper = match_score as i64 + 2 * gap_open as i64;
    let go = gap_open as i64 - gap_extend as i64;
    upper + go <= u8::MAX as i64
}

/// One anti-diagonal's four difference vectors, indexed by row `i`.
///
/// Grow-only, and padded by `LANES` so the partial last chunk can store all 32 lanes without a
/// bounds check.
#[derive(Default)]
struct Deltas {
    dh: Vec<u8>,
    dv: Vec<u8>,
    de: Vec<u8>,
    df: Vec<u8>,
}

impl Deltas {
    fn grow(&mut self, rows: usize) {
        for b in [&mut self.dh, &mut self.dv, &mut self.de, &mut self.df] {
            if b.len() < rows {
                b.resize(rows, 0);
            }
        }
    }
}

/// Reusable scratch, owned by the aligner and reused across calls.
///
/// Every buffer here is **grow-only**: sized up when a call needs more and never shrunk, so a
/// steady stream of similar-sized alignments allocates once. Small alignments dominate the
/// intended workload, where a per-call allocation would dominate the call.
#[derive(Default)]
pub struct U8Probe {
    q_codes: Vec<u8>,
    r_rev: Vec<u8>,
    prev: Deltas,
    cur: Deltas,
    /// Traceback flags, `TB_WORDS` `u32`s per 32-lane chunk: **4 bits per cell**, laid out in
    /// `(p, lane)` coordinates rather than `(i, j)`.
    ///
    /// Under a band this stores only the cells actually visited, not the full matrix.
    tb: Vec<u32>,
    /// Chunk index at which anti-diagonal `p` starts within `tb`. Anti-diagonals have unequal
    /// lengths - the ramp at both matrix corners - so a fixed stride would waste ~50%.
    tb_base: Vec<u32>,

    /// The per-row running **absolute** score: `abs[i] == H[i][p - i]` for the anti-diagonal
    /// being filled. See `TRACK_ROW_MAX` on `run` for why this is exact and needs no bound.
    ///
    /// `i32`, not `u8`: the DP itself stays in u8 lanes: this is a reconstruction carried
    /// alongside it for the modes that argmax over absolute scores, so `fits_u8_lanes` is
    /// untouched and still has no length term.
    abs: Vec<i32>,
    /// `row_max[i] == max over j of H[i][j]`.
    row_max: Vec<i32>,

    /// **Where** each row's maximum was attained, as one bit per cell: set iff that cell took
    /// its row's running max. One `u32` per 32-lane chunk, in the same `(p, chunk)` coordinates
    /// as `tb`, so `tb_base[p]` indexes both.
    ///
    /// The argmax is only ever asked about one row, so carrying a winning `j` per row would
    /// compute one for every row on every anti-diagonal and throw all but one away. The bit is
    /// enough: `arg_for_row` recovers the row that matters afterwards, scalar.
    ///
    /// Recovery is exact, and `take` is why: it is the update condition itself, so "the last
    /// cell whose bit is set" is by construction the cell a blend would have left, under either
    /// tie-break (`ARG_LARGEST` makes `take` non-strict so later ties win; otherwise it is
    /// strict, so the first winner sticks).
    nm: Vec<u32>,
    /// Row 0's argmax, which the bit stream cannot carry: row 0 is never *interior* to an
    /// anti-diagonal (`i_from >= 1` always), so it has no lane and no bit. It is maintained
    /// scalar in the boundary, exactly as it always was, and it is load-bearing - with
    /// `gap_open == gap_extend == 0` the whole of row 0 ties at 0 and the longest-extension
    /// tie-break has to take the furthest cell.
    row0_arg: i32,

    /// A second set of flags and row maxima, for `split_reference_alignment` and nothing else.
    ///
    /// That mode runs the kernel **twice** - once forward against the right reference, once
    /// REVERSED against the left - and then replays **both** arms, so both tracebacks have to
    /// be alive at the same time. Empty for every other mode, and grow-only like the rest, so
    /// a program that never splits never pays for these.
    tb_b: Vec<u32>,
    tb_base_b: Vec<u32>,
    row_max_b: Vec<i32>,
    nm_b: Vec<u32>,
    row0_arg_b: i32,
}

/// Recover `argmax over j of H[row][j]` from the new-max bit stream.
///
/// Takes its buffers explicitly rather than reaching through `self`, because `split_reference`
/// holds two passes' worth at once and has to say which.
///
/// The stream's bit is the fill's `take` mask, so the cells with the bit set are every cell that
/// *became* row `row`'s max, and the last of them is the one still holding it when the fill
/// ended. Walking backwards therefore needs no knowledge of the tie-break: `ARG_LARGEST` makes
/// `take` non-strict, so later ties re-fire the bit and the last wins; otherwise `take` is
/// strict, so ties do not re-fire and the first winner stays.
///
/// Row `row` is interior to anti-diagonal `p` exactly for `p` in `[row+1, row+rlen]`, so those
/// are the only `p` scanned - which is also why the garbage the partial chunk tail writes above
/// `i_to` is never read: it lands at `p <= row`, below where this starts.
///
/// Falling off the bottom means no interior cell ever took the max, so the row's max is still
/// its column-0 seed: `j = 0`.
///
/// The scan is two phases because `lane = row - max(1, p - rlen)` is constant at `row - 1` while
/// `p <= rlen + 1`, which is most of it; hoisting the lane out of that phase leaves a load, a
/// shift, a load and a test per step.
fn arg_for_row(nm: &[u32], tb_base: &[u32], row0_arg: i32, row: usize, rlen: usize) -> usize {
    // Row 0 has no lane and no bit - `i_from >= 1` always - so it is carried scalar.
    if row == 0 {
        return row0_arg as usize;
    }
    let word_at = |p: usize, lane: usize| -> bool {
        // SAFETY: `p <= row + rlen <= qlen + rlen`, which is what `size_flags` sized `tb_base`
        // to, and `tb_base[p] + lane/LANES` is a chunk the fill wrote - `row` is interior to
        // every `p` this scans (see below), so its lane is one the group loop covered.
        unsafe {
            let w = *nm.get_unchecked(*tb_base.get_unchecked(p) as usize + lane / LANES);
            w >> (lane % LANES) & 1 != 0
        }
    };

    let mut p = row + rlen;
    // Phase 2: `p > rlen + 1`, so `i_from = p - rlen` and the lane slides with `p`. At most
    // `qlen` iterations - the ramp off the end of the reference.
    while p > rlen + 1 && p > row {
        if word_at(p, row + rlen - p) {
            return p - row;
        }
        p -= 1;
    }
    // Phase 1: `i_from == 1`, so the lane is fixed and so is everything derived from it.
    let lane = row - 1;
    let off = lane / LANES;
    let bit = 1u32 << (lane % LANES);
    while p > row {
        // SAFETY: as above.
        if unsafe { *nm.get_unchecked(*tb_base.get_unchecked(p) as usize + off) } & bit != 0 {
            return p - row;
        }
        p -= 1;
    }
    0
}

/// [`arg_for_row`] under a band: `argmax over j of H[row][j]`, from the same bit stream.
///
/// A separate scan rather than a clamp on the other one, because bounding the range here is a
/// **correctness** requirement, not an optimisation. Row `row` is in band only for `p` in
/// `[2*row - w, 2*row + w]`. Above that the row sits above the band, `row - lo(p)` goes negative
/// and the unchecked read lands anywhere; below it the row sits inside the garbage the partial
/// chunk's tail wrote, and `lane / LANES` can walk into the next anti-diagonal's chunks and read
/// a bit that means something else. The unbanded scan needs neither bound because the unbanded
/// matrix has neither edge.
///
/// Falling off the bottom returns the row's seed column, which under a band is `row - w`, not 0:
/// row `row > w` has no column-0 cell and is seeded at its first in-band column by the fill's
/// band-entry seed. `saturating_sub` covers both cases.
///
/// That also disposes of the one garbage bit in the stream. At `p = 2*row - w` the row's `take`
/// came from a `ΔV` whose left neighbour is outside the band, so the bit means nothing - but the
/// scan only reaches it after finding no set bit above, which is exactly the case where the
/// answer is the seed column anyway. Set or clear, it returns `row - w`, so the fill leaves it
/// alone rather than paying a read-modify-write to correct it.
fn arg_for_row_banded(
    nm: &[u32],
    tb_base: &[u32],
    row0_arg: i32,
    row: usize,
    qlen: usize,
    rlen: usize,
    w: usize,
) -> usize {
    if row == 0 {
        return row0_arg as usize;
    }
    // Row `row` is interior to `p` for `j = p - row` in `[1, rlen]` and `|2*row - p| <= w`.
    let p_hi = (row + rlen).min(2 * row + w);
    let p_lo = (row + 1).max((2 * row).saturating_sub(w));

    let mut p = p_hi;
    while p >= p_lo {
        let (lo, _) = band_rows::<true>(p, qlen, rlen, w);
        let lane = row - lo;
        // SAFETY: `p` is in `[p_lo, p_hi]`, so `row` is interior to it and `lane` is a lane the
        // group loop covered; `p <= row + rlen <= qlen + rlen` is what `size_flags` sized
        // `tb_base` to.
        let set = unsafe {
            let word = *nm.get_unchecked(*tb_base.get_unchecked(p) as usize + lane / LANES);
            word >> (lane % LANES) & 1 != 0
        };
        if set {
            return p - row;
        }
        p -= 1;
    }
    row.saturating_sub(w)
}

/// The alphabet: **A, C, G, T, U**, and one code for everything else.
///
/// Codes 0-3 are A/C/G/T, case-insensitive, with `U` folded onto `T`. Code 4 is *unknown* -
/// every other byte - and it must never match anything, **not even another unknown byte**.
///
/// There is no substitution matrix here, and that is the point. The scheme is `match` on the
/// diagonal and `-mismatch` everywhere else: a binary *predicate*, not a table. So scoring a
/// cell is one `cmpeq` against the codes already in registers - no gather, no `pshufb`, no
/// 25-entry lookup.
///
/// The one thing a plain compare cannot express is `unknown != unknown`, because `x == x` is
/// unavoidably true. Hence the **asymmetry**: [`U8Probe::prepare`] re-encodes unknown to code 5
/// on the *reference* side only, so the compare says "not equal" and the blend picks
/// `-mismatch`. One `cmpeq` and one `blendv`, exactly what a 2-bit ACGT-only encoding would
/// cost - so the fail-safe is free, and a stray `-` or IUPAC `R` scores a mismatch rather than
/// silently becoming an `A` and scoring a false match. See `mod.rs`, "Alphabet".
fn code(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

impl U8Probe {
    pub fn new() -> Self {
        Self::default()
    }

    /// Check the preconditions, encode both sequences, and size the delta buffers.
    ///
    /// Assumes both sequences are non-empty - a zero-length side
    /// is one pure gap and never reaches the kernel.
    fn prepare<const REVERSED: bool>(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        match_score: u8,
        gap_open: u8,
        gap_extend: u8,
    ) {
        assert!(
            fits_u8_lanes(match_score, gap_open, gap_extend),
            "scores do not fit u8 lanes"
        );
        assert!(gap_open >= gap_extend, "gap_open must be >= gap_extend");
        assert!(is_x86_feature_detected!("avx2"), "AVX2 is required");

        let (qlen, rlen) = (query.len(), refseq.len());

        // Encoding is O(n) and happens every call regardless, which is why the asymmetric
        // unknown-base encoding costs nothing.
        //
        // Grow-then-overwrite, not `clear()` + `extend()`: the latter re-checks capacity per
        // element.
        // REVERSED is free, and this is the whole reason: the sequences have to be encoded
        // every call anyway, so a `_start` mode just encodes them **back-to-front**. Same pass,
        // same cost, no buffer, and the kernel never learns the mode exists.
        //
        // Note the two axes compose: the kernel always wants the reference laid out backwards
        // relative to the query (that is what makes `refseq[p-i-1]` a contiguous load), so
        // under REVERSED the reference is reversed *twice* and comes out forward.
        if self.q_codes.len() < qlen + LANES {
            self.q_codes.resize(qlen + LANES, 4);
        }
        if REVERSED {
            for (dst, &b) in self.q_codes[..qlen].iter_mut().zip(query.iter().rev()) {
                *dst = code(b);
            }
        } else {
            for (dst, &b) in self.q_codes[..qlen].iter_mut().zip(query) {
                *dst = code(b);
            }
        }
        self.q_codes[qlen..qlen + LANES].fill(4);

        // Code 5 is N on the reference side - the asymmetry that makes `N != N` fall out of
        // the compare.
        if self.r_rev.len() < rlen + LANES {
            self.r_rev.resize(rlen + LANES, 5);
        }
        let enc_r = |b: u8| -> u8 {
            let c = code(b);
            if c == 4 { 5 } else { c }
        };
        if REVERSED {
            for (dst, &b) in self.r_rev[..rlen].iter_mut().zip(refseq.iter()) {
                *dst = enc_r(b);
            }
        } else {
            for (dst, &b) in self.r_rev[..rlen].iter_mut().zip(refseq.iter().rev()) {
                *dst = enc_r(b);
            }
        }
        self.r_rev[rlen..rlen + LANES].fill(5);

        let rows = qlen + 1 + LANES;
        self.prev.grow(rows);
        self.cur.grow(rows);
        // Each buffer gets its OWN check. Gating all three on `abs.len()` looks equivalent and
        // is not: `split_reference` swaps `row_max` with its `_b` twin between its two passes,
        // so `abs` can already be grown while the row buffers are the empty ones - the check
        // would not fire and the kernel would write off the end.
        if self.abs.len() < rows {
            self.abs.resize(rows, 0);
        }
        if self.row_max.len() < rows {
            self.row_max.resize(rows, 0);
        }
    }

    /// Global alignment: both sequences spanned end to end.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: the scheme plus the bandwidth
    pub fn global_alignment(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        match bandwidth {
            // A band that reaches every cell is the exact mode, and the exact mode is cheaper -
            // see `band_reaches_everything`.
            Some(w) if !band_reaches_everything(query.len(), refseq.len(), w) => {
                self.global_banded(query, refseq, m, mm, go, ge, w)
            }
            _ => self.align::<false, false>(query, refseq, m, mm, go, ge),
        }
    }

    /// The banded `global` path: `S[qlen][rlen]`, confined to `|i - j| <= w`.
    ///
    /// # Widening: global is pinned at both corners
    ///
    /// Every other mode here has exactly one anchor and lets the far end float, so `|i-j| <= w`
    /// is a statement about drift away from that anchor and any `w` is meaningful. Global must
    /// also *arrive* at `(qlen, rlen)`, which sits `|qlen - rlen|` off the diagonal the band is
    /// centred on. A narrower band does not contain a bad alignment, it contains **no path at
    /// all**: `band_rows` returns `lo > hi` for the final anti-diagonal and the corner is never
    /// computed. Hence the widening to `|qlen - rlen|` - the narrowest diagonal band that still
    /// joins the two corners.
    ///
    /// That widening is free when `qlen ~ rlen`, which is the shape global is normally handed
    /// (two segments already believed to correspond end to end), and it is what makes the band
    /// worth having here at all. It bites on a lopsided pair, where it degrades to
    /// correct-but-slow rather than to wrong.
    ///
    /// # The score comes from `abs`, not `anchor`
    ///
    /// See [`U8Probe::local_reference`]: `align`'s `anchor` walk cannot survive a band. But
    /// `TRACK_ROW_MAX` already maintains `abs[i]`, the running **absolute** score of row `i`'s
    /// cell on the current anti-diagonal, seeded at the row's band entry and advanced by one
    /// `ΔV` per cell. Row `qlen` is live for `p` up to `qlen + rlen`, and the widening guarantees
    /// it is still in band there - so when the fill returns, `abs[qlen]` *is* `S[qlen][rlen]`,
    /// with nothing to add. The end cell is read, not recomputed.
    #[allow(clippy::too_many_arguments)] // ditto, plus the bandwidth
    fn global_banded(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        w: usize,
    ) -> AlignmentResult {
        let (qlen, rlen) = (query.len(), refseq.len());
        // Band-independent given the widening below: with one side empty the only path is the
        // single gap `degenerate` returns, and `w >= |qlen - rlen|` is exactly wide enough to
        // hold it.
        if qlen == 0 || rlen == 0 {
            return self.degenerate::<false, false>(qlen, rlen, go, ge);
        }
        let w = w.max(1).max(qlen.abs_diff(rlen));

        self.row_pass::<false, true>(query, refseq, m, mm, go, ge, w);

        // `abs[qlen]` after the fill is `S[qlen][rlen]` - see this function's docs. The widening
        // is what puts the corner in band, and so what makes this read meaningful.
        debug_assert!(
            qlen.abs_diff(rlen) <= w,
            "the widening should have put the corner in band"
        );
        let score = self.abs[qlen];

        let cigar = self.replay::<false, true>(query, refseq, qlen, rlen, w);
        AlignmentResult {
            score,
            query_start: 0,
            query_end: qlen,
            ref_start: 0,
            ref_end: rlen,
            cigar,
        }
    }

    /// The query is spanned in full; the **reference may stop anywhere**, and its trailing
    /// tail is free and does not appear in the CIGAR.
    ///
    /// `j = 0` is a legal end cell and sometimes the winning one. Ties resolve to the
    /// **smallest** `ref_end`.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: the scheme plus the bandwidth
    pub fn local_reference_end_alignment(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        match bandwidth {
            Some(w) if !band_reaches_everything(query.len(), refseq.len(), w) => {
                self.local_reference::<false>(query, refseq, m, mm, go, ge, w)
            }
            _ => self.align::<true, false>(query, refseq, m, mm, go, ge),
        }
    }

    /// The mirror: the reference may **begin** anywhere and its leading prefix is free.
    ///
    /// This is `local_reference_end_alignment` on reversed inputs, and that is exact rather
    /// than an approximation: affine gap costs are symmetric under reversal (a gap of length
    /// `L` costs `open + (L-1)*extend` whichever way you read it) and `=`/`X` is literal byte
    /// identity, also symmetric. Ties resolve to the **largest** `ref_start`.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: the scheme plus the bandwidth
    pub fn local_reference_start_alignment(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        match bandwidth {
            Some(w) if !band_reaches_everything(query.len(), refseq.len(), w) => {
                self.local_reference::<true>(query, refseq, m, mm, go, ge, w)
            }
            _ => self.align::<true, true>(query, refseq, m, mm, go, ge),
        }
    }

    /// The banded `local_reference_*` path: `max over j of H[qlen][j]`, confined to `|i-j| <= w`.
    ///
    /// # Why not `align` with `BANDED` switched on
    ///
    /// `align` does not read its score out of the fill. It walks `anchor`, a closed-form
    /// recurrence along row `qlen` that reads `ΔV[qlen]` on **every** anti-diagonal past
    /// `qlen`, and row `qlen` is a horizontal line while the band is a diagonal strip, so the
    /// two meet on only `O(w)` of those. Under a band the rest of that walk reads lanes the
    /// fill never wrote. `BANDED` and `TRACK_LAST_ROW` were never compiled together for
    /// exactly this reason.
    ///
    /// `TRACK_ROW_MAX` is the way in: it carries **absolute** per-cell scores, it is already
    /// banded (it is what `local` runs), and `row_max[qlen]` *is* this mode's objective. So the
    /// banded path is a `row_pass` and a reduce that is one array read, and the anchor never
    /// enters it.
    ///
    /// It costs more per cell than `align` - a `row_max` store and a bit of the new-max stream
    /// per chunk - which is why `None` still takes `align`: the two are separate
    /// monomorphizations, so the exact path never pays for this one. The band's own win (the
    /// fill drops from `qlen*rlen` to `(qlen+rlen)*(w+1)` cells) buys that back many times.
    #[allow(clippy::too_many_arguments)] // ditto, plus the bandwidth
    fn local_reference<const REVERSED: bool>(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        w: usize,
    ) -> AlignmentResult {
        let (qlen, rlen) = (query.len(), refseq.len());
        // Band-independent: see the widening below, which always admits the whole of the one
        // live row or column these cases consist of.
        if qlen == 0 || rlen == 0 {
            return self.degenerate::<true, REVERSED>(qlen, rlen, go, ge);
        }

        // The widening. This mode must span the query, so it must reach a cell `(qlen, j)` with
        // `j <= rlen`; in band that needs `qlen - w <= rlen`. Below that the band contains no
        // alignment at all - not a bad one, *none* - so `Some(w)` is honoured as "at least w",
        // and widened to the narrowest band that has an answer. It only bites when the query is
        // longer than the reference, which this mode is not normally handed.
        let w = w.max(1).max(qlen.saturating_sub(rlen));

        self.row_pass::<REVERSED, true>(query, refseq, m, mm, go, ge, w);

        // `last_row == min(qlen, rlen + w) == qlen` after the widening, so row `qlen` is filled
        // and this read is in bounds - that is what the widening is for.
        debug_assert!(
            qlen <= rlen + w,
            "the widening should have made row qlen reachable"
        );
        let score = self.row_max[qlen];
        // `ARG_LARGEST = false` in `row_pass`, so ties give the smallest `j` - this mode's rule.
        let j = arg_for_row_banded(&self.nm, &self.tb_base, self.row0_arg, qlen, qlen, rlen, w);
        debug_assert!(
            qlen.abs_diff(j) <= w,
            "the reduce picked ({qlen}, {j}), which is outside the band of {w}"
        );

        let cigar = self.replay::<REVERSED, true>(query, refseq, qlen, j, w);
        let (ref_start, ref_end) = ref_span::<REVERSED>(rlen, j);
        AlignmentResult {
            score,
            query_start: 0,
            query_end: qlen,
            ref_start,
            ref_end,
            cigar,
        }
    }

    /// Ends wherever it scores best: `max over (i, j) of H[i][j] + (end_bonus if i == qlen)`.
    ///
    /// Cell `(0,0)` is a legal end cell, so **the score is never negative** - a hopeless
    /// extension returns the empty alignment (spans `0..0`, score 0, empty CIGAR).
    ///
    /// `end_bonus` never enters the DP. It is applied once, at the end, to choose between "the
    /// best cell anywhere" and "the best cell that finishes the query". Ties resolve to the
    /// **longest extension**: largest `query_end`, then largest `ref_end`.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: the scheme is five scalars
    pub fn local_end_alignment(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        end_bonus: u32,
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        match bandwidth {
            Some(w) if !band_reaches_everything(query.len(), refseq.len(), w.max(1)) => {
                self.local::<false, true>(query, refseq, m, mm, go, ge, end_bonus, w.max(1))
            }
            _ => self.local::<false, false>(query, refseq, m, mm, go, ge, end_bonus, 0),
        }
    }

    /// The mirror: begins wherever it scores best, both leading prefixes free.
    ///
    /// Its empty alignment sits at the *end* of both sequences (`qlen..qlen`, `rlen..rlen`) -
    /// the mirror of `local_end`'s, which sits at the start. Tie-break mirrored too: the
    /// **smallest** `query_start`, then the smallest `ref_start`.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: the scheme is five scalars
    pub fn local_start_alignment(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        end_bonus: u32,
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        match bandwidth {
            Some(w) if !band_reaches_everything(query.len(), refseq.len(), w.max(1)) => {
                self.local::<true, true>(query, refseq, m, mm, go, ge, end_bonus, w.max(1))
            }
            _ => self.local::<true, false>(query, refseq, m, mm, go, ge, end_bonus, 0),
        }
    }

    /// `w` is read only when `BANDED`, and is the caller's bandwidth already clamped to
    /// `>= 1`: a zero-width band is not meaningful, so callers clamp rather than pass it on.
    #[allow(clippy::too_many_arguments)] // ditto, plus the bandwidth
    fn local<const REVERSED: bool, const BANDED: bool>(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        match_score: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
        end_bonus: u32,
        w: usize,
    ) -> AlignmentResult {
        let (qlen, rlen) = (query.len(), refseq.len());
        // A zero-length side skips the kernel - there is no anti-diagonal to walk - but it is
        // NOT simply "the empty alignment", and assuming so was a real bug:
        //
        //  * with `qlen == 0`, EVERY cell has `i == qlen`, so `end_bonus` applies to all of
        //    them. The empty alignment scores `end_bonus`, not 0.
        //  * with free gaps (`gap_open == gap_extend == 0`) the whole boundary row/column ties
        //    at 0, and the longest-extension tie-break must then take the FURTHEST cell, not
        //    the origin.
        //
        // So scan the one live row or column in closed form. O(qlen + rlen), and only here.
        if qlen == 0 || rlen == 0 {
            let gap = |k: usize| -> i32 {
                if k == 0 {
                    0
                } else {
                    -(gap_open as i32) - (k as i32 - 1) * gap_extend as i32
                }
            };
            // The band reaches the live row/column too, and it bites here for the same reason
            // it bites anywhere: `H[0][j]` sits `j` off the diagonal, so a band of `w` can only
            // see the first `w` of it. Without this the free-gap case (`gap_open == gap_extend
            // == 0`, where the whole row ties at 0) would report the furthest cell in the row
            // rather than the furthest cell in the *band*.
            let reach = |k: usize| if BANDED { k.min(w) } else { k };
            let (mut best, mut bi, mut bj) = (i32::MIN, 0usize, 0usize);
            if qlen == 0 {
                // Row 0: H[0][j] is j deletions, and i == qlen == 0 throughout, so the bonus
                // is on every candidate. `>=` takes the largest j on ties.
                for j in 0..=reach(rlen) {
                    let v = gap(j) + end_bonus as i32;
                    if v >= best {
                        best = v;
                        bj = j;
                    }
                }
            } else {
                // Column 0: H[i][0] is i insertions; only i == qlen earns the bonus. `>=`
                // takes the largest i on ties.
                //
                // The bonus needs the band's permission: it lands on `i == qlen`, and under a
                // band that cell only exists if the band reaches it. Where it does not, the
                // loop simply never visits `i == qlen` and no bonus is ever added, which is
                // right - the query cannot be spanned inside the band at all.
                for i in 0..=reach(qlen) {
                    let v = gap(i) + if i == qlen { end_bonus as i32 } else { 0 };
                    if v >= best {
                        best = v;
                        bi = i;
                    }
                }
            }
            let mut cigar = Cigar::default();
            if bi > 0 {
                cigar.push(CigarOperation::Insertion, bi);
            }
            if bj > 0 {
                cigar.push(CigarOperation::Deletion, bj);
            }
            let (query_start, query_end) = query_span::<REVERSED>(qlen, bi);
            let (ref_start, ref_end) = ref_span::<REVERSED>(rlen, bj);
            return AlignmentResult {
                score: best,
                query_start,
                query_end,
                ref_start,
                ref_end,
                cigar,
            };
        }

        self.prepare::<REVERSED>(query, refseq, match_score, gap_open, gap_extend);
        self.size_flags::<true, BANDED>(qlen, rlen, w);

        unsafe {
            self.run::<true, false, true, true, BANDED>(
                qlen,
                rlen,
                match_score,
                mismatch,
                gap_open,
                gap_extend,
                w,
            )
        };

        // How far down the query the band still has cells to offer.
        //
        // Row `i`'s in-band columns are `[i - w, i + w]` intersected with `[0, rlen]`, which is
        // empty once `i > rlen + w`. Those rows are not merely bad candidates, they were never
        // filled: `row_max[i]` there is whatever the previous call left, or the garbage the
        // partial chunk's tail wrote. So the reduce must not look at them - a wrong-but-huge
        // value would win, and would win *silently*.
        //
        // Only reachable when the query is longer than the reference by more than the band, so
        // never in the intended use (an end extension is handed more reference than query) -
        // which is exactly why the bound is enforced rather than assumed.
        let last_row = if BANDED { qlen.min(rlen + w) } else { qlen };

        // Reduce the per-row maxima to a single cell. Scalar, and cheap: qlen+1 values against
        // the fill's qlen*rlen cells.
        //
        // Ties resolve to the longest extension - largest i, then largest j. `>=` on i gives the
        // largest i directly; the largest j comes from `arg_for_row`, whose bit stream was
        // written under `ARG_LARGEST`'s non-strict `take` and so records the last tie.
        let (mut best, mut bi) = (0i32, 0usize);
        for i in 0..=last_row {
            let v = self.row_max[i];
            if v >= best {
                best = v;
                bi = i;
            }
        }
        // The end_bonus candidate: the best cell that finishes the query. Applied once, here,
        // never in the DP.
        //
        // Under a band it needs `row_max[qlen]` to exist at all. Where the band cannot reach
        // row `qlen` there is no cell that finishes the query, so there is nothing for the
        // bonus to apply to - and `row_max[qlen]` is unfilled, so asking anyway would read
        // garbage and hand it the answer.
        if last_row == qlen {
            let finish = self.row_max[qlen] + end_bonus as i32;
            if finish > best || (finish == best && qlen >= bi) {
                best = finish;
                bi = qlen;
            }
        }
        // The argmax `j`, recovered for the ONE row that won - which is the whole reason the
        // fill stores a bit per cell instead of a `j` per row. `O(rlen)` unbanded, `O(w)` under
        // a band.
        let bj = if BANDED {
            arg_for_row_banded(&self.nm, &self.tb_base, self.row0_arg, bi, qlen, rlen, w)
        } else {
            arg_for_row(&self.nm, &self.tb_base, self.row0_arg, bi, rlen)
        };
        // No early return for `best == 0`, and that is not an oversight. (0,0) is the reduce's
        // seed candidate, so the score can never go negative - the "never negative" guarantee
        // falls out rather than being clamped. But bailing out at `best <= 0` would ALSO
        // discard a *longer* extension that happens to score exactly 0, which the
        // longest-extension tie-break says must win over (0,0) - a right score with a wrong
        // endpoint, which is what the endpoint assertion below exists to catch.
        debug_assert!(
            best >= 0,
            "(0,0) seeds the reduce, so the score cannot go negative"
        );
        debug_assert!(
            !BANDED || bi.abs_diff(bj) <= w,
            "the reduce picked ({bi}, {bj}), which is outside the band of {w}"
        );

        let cigar = self.replay::<REVERSED, BANDED>(query, refseq, bi, bj, w);
        let (query_start, query_end) = query_span::<REVERSED>(qlen, bi);
        let (ref_start, ref_end) = ref_span::<REVERSED>(rlen, bj);
        AlignmentResult {
            score: best,
            query_start,
            query_end,
            ref_start,
            ref_end,
            cigar,
        }
    }

    /// An upper bound on the chunks a **banded** fill writes, and therefore on the `tb`/`nm` it
    /// needs.
    ///
    /// # Derived from the fill's own accounting, not from a cell count
    ///
    /// `run` advances its chunk counter exactly once per live anti-diagonal, by
    /// `(hi - lo) / LANES + 1` (see the `if TRACE && lo <= hi` in the fill). So the total is
    /// `sum over live p of ((hi_p - lo_p) / LANES + 1)`, and bounding it needs only a bound on
    /// `hi - lo`:
    ///
    /// `band_rows::<true>` returns `ceil((p - w)/2) ..= floor((p + w)/2)` before the matrix
    /// clamps, which spans `w + 1` rows when `p + w` is even and `w` when it is odd - so
    /// `hi - lo <= w` always, and the matrix clamps only shrink it. There are at most
    /// `qlen + rlen` live anti-diagonals (`p = 0` has no interior). Hence
    ///
    /// ```text
    /// chunks <= (qlen + rlen) * (w / LANES + 1)
    /// ```
    ///
    /// **`w / LANES + 1`, not `(w + 1) / LANES`.** A strip of `w + 1` rows can straddle a chunk
    /// boundary, and the fill starts each anti-diagonal's chunks at the *band's* top row rather
    /// than the matrix's, so the `+ 1` is the partial chunk and is not slack to be tuned away.
    /// The `+ 2` on top is the same margin the unbanded bound carries.
    ///
    /// **This is a bound the fill must not exceed, and exceeding it is heap corruption rather
    /// than a wrong answer** - the fill's stores are unchecked. `run` debug-asserts the final
    /// chunk count against the buffer it was actually handed, which is what makes this derivation
    /// checked rather than argued - so run the banded paths in a debug build after touching it.
    ///
    /// At `q=20000 r=24000 w=64` this is 132,002 chunks against the full matrix's 15,044,002:
    /// **2.6 MB instead of 300 MB**, for a fill that only ever visits 2.9M cells.
    #[inline]
    fn banded_chunks_bound(qlen: usize, rlen: usize, w: usize) -> usize {
        (qlen + rlen) * (w / LANES + 1) + 2
    }

    /// Size the traceback flag buffer and its per-anti-diagonal index. Shared by every mode.
    ///
    /// `NEED_NM` additionally sizes the new-max bit stream, which only the `TRACK_ROW_MAX` modes
    /// write. A const generic, not a runtime flag, for the usual reason: `global` must not carry
    /// so much as a branch for a buffer it never touches, and grow-only means a program that
    /// only ever aligns globally never allocates it at all.
    ///
    /// `BANDED`/`w` are the same deal, and they are what stops a band allocating for a matrix it
    /// will never visit - see [`U8Probe::banded_chunks_bound`].
    fn size_flags<const NEED_NM: bool, const BANDED: bool>(
        &mut self,
        qlen: usize,
        rlen: usize,
        w: usize,
    ) {
        // Size the flag buffer. Anti-diagonals have unequal lengths - the ramp at both matrix
        // corners - so a fixed stride would waste ~50%; `tb_base` records where each starts,
        // and the fill writes it as it goes.
        //
        // An exact count would need its own O(qlen + rlen) pass, which at small n is
        // comparable to the entire fill. So bound it instead: every cell lives on exactly one
        // anti-diagonal, giving `qlen*rlen/LANES` chunks, plus at most one partial chunk per
        // anti-diagonal. Over-allocates by ~(qlen+rlen) chunks, and the buffers are grow-only
        // anyway.
        let full_bound = qlen * rlen / LANES + qlen + rlen + 2;
        // Under a band, take whichever bound is smaller. **Both are valid** - a banded fill
        // visits a subset of the cells an unbanded one does, so `full_bound` never stops being
        // an upper bound - and the `min` matters in both directions: `banded_chunks_bound`
        // *exceeds* `full_bound` once `w` approaches the matrix, because it charges a partial
        // chunk per anti-diagonal against a strip that is no longer narrow.
        let chunks_bound = if BANDED {
            full_bound.min(Self::banded_chunks_bound(qlen, rlen, w))
        } else {
            full_bound
        };
        if self.tb_base.len() < qlen + rlen + 1 {
            self.tb_base.resize(qlen + rlen + 1, 0);
        }
        if self.tb.len() < chunks_bound * TB_WORDS {
            self.tb.resize(chunks_bound * TB_WORDS, 0);
        }
        if NEED_NM && self.nm.len() < chunks_bound * NM_WORDS {
            self.nm.resize(chunks_bound * NM_WORDS, 0);
        }
    }

    /// One query aligned across **two** references, with a single jump point between them.
    ///
    /// Returns two [`AlignmentResult`]s that **partition the query exactly**: `right` covers
    /// `query[..k]` from `right_reference[0]`, `left` covers `query[k..]` to the end of
    /// `left_reference`, and `k` is chosen to maximize the total.
    ///
    /// # Reduction to two ordinary passes
    ///
    /// Over the same `H` matrix every other mode computes:
    ///
    /// ```text
    /// F[k] = max over j of H(query, right)[k][j]                  - the best right arm handing over at k
    /// G[k] = max over l of (best alignment of query[k..] vs left[l..])  - the best left arm taking over
    /// score = max over k in 0..=qlen of F[k] + G[k]
    /// ```
    ///
    /// `F` is exactly `TRACK_ROW_MAX`'s output. And `G` is `F` on **reversed** inputs, so it is
    /// the same pass with `REVERSED` - free, as always. The arms interact *only* through `k`,
    /// so the whole "zig-zag" decision is a one-dimensional scan over `qlen + 1` values.
    /// Everything before it is the ordinary kernel, run twice.
    ///
    /// # Ties
    ///
    /// Resolved to the smallest `|left.score - right.score|`, then the smallest `k`. So if the
    /// other arm can pick up a base for free, it will.
    ///
    /// # Panics
    /// If `fits_u8_lanes` rejects the scheme, if `gap_open < gap_extend`, or without AVX2.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: two references and the scheme
    pub fn split_reference_alignment(
        &mut self,
        query: &[u8],
        left_reference: &[u8],
        right_reference: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        bandwidth: Option<usize>,
    ) -> SplitReferenceAlignment {
        match bandwidth {
            Some(w) => {
                // One band, both arms - each measured from its own arm's anchor. See
                // `SimdAligner::split_reference_alignment`.
                //
                // The widening: the two arms must between them cover the whole query, and an arm
                // can reach at most `reflen + w` query bases inside the band, so a `w` below
                // `(qlen - llen - rrlen) / 2` admits no split at all. Widening to it is the
                // narrowest band that still has an answer.
                let deficit = query
                    .len()
                    .saturating_sub(left_reference.len() + right_reference.len());
                let w = w.max(1).max(deficit.div_ceil(2));
                // A band that reaches every cell is the exact mode, and the exact mode is
                // cheaper - see `band_reaches_everything`. Both arms have to be reached, and
                // they align against different references, so the wider of the two decides.
                let widest = left_reference.len().max(right_reference.len());
                if band_reaches_everything(query.len(), widest, w) {
                    return self.split_reference::<false>(
                        query,
                        left_reference,
                        right_reference,
                        m,
                        mm,
                        go,
                        ge,
                        0,
                    );
                }
                self.split_reference::<true>(
                    query,
                    left_reference,
                    right_reference,
                    m,
                    mm,
                    go,
                    ge,
                    w,
                )
            }
            None => self.split_reference::<false>(
                query,
                left_reference,
                right_reference,
                m,
                mm,
                go,
                ge,
                0,
            ),
        }
    }

    #[allow(clippy::too_many_arguments)] // ditto, plus the bandwidth
    fn split_reference<const BANDED: bool>(
        &mut self,
        query: &[u8],
        left_reference: &[u8],
        right_reference: &[u8],
        m: u8,
        mm: u8,
        go: u8,
        ge: u8,
        w: usize,
    ) -> SplitReferenceAlignment {
        let qlen = query.len();

        // --- pass A: forward against the right reference. row_max[k] == F[k]. ---
        self.row_pass::<false, BANDED>(query, right_reference, m, mm, go, ge, w);
        // Park A's results in the `_b` slots and let pass B use the primary ones. A swap of Vec
        // headers, once per call - not per anti-diagonal, where it would matter.
        std::mem::swap(&mut self.tb, &mut self.tb_b);
        std::mem::swap(&mut self.tb_base, &mut self.tb_base_b);
        std::mem::swap(&mut self.row_max, &mut self.row_max_b);
        // `nm` and `row0_arg` travel with `tb_base`, and they have to: `arg_for_row` reads all
        // three, and reading one pass's bit stream against the other pass's chunk index would
        // silently return a `j` from the wrong matrix.
        std::mem::swap(&mut self.nm, &mut self.nm_b);
        std::mem::swap(&mut self.row0_arg, &mut self.row0_arg_b);
        // From here: `*_b` holds pass A (right arm), `*` holds pass B (left arm).

        // --- pass B: REVERSED against the left reference. row_max[qlen - k] == G[k]. ---
        self.row_pass::<true, BANDED>(query, left_reference, m, mm, go, ge, w);

        // --- the jump point: a 1-D scan over qlen+1 values, against the fill's qlen*rlen ---
        //
        // Under a band the scan must not look at every `k`. Pass A filled rows `0..=min(qlen,
        // rrlen + w)` and pass B rows `0..=min(qlen, llen + w)`, and `k` reads `row_max_b[k]`
        // from A and `row_max[qlen - k]` from B, so both bounds bite - from opposite ends.
        // Outside them the value is not a bad candidate, it is the previous call's leftovers,
        // and being unbounded it would win. `split_reference_alignment` widened `w` so this
        // range is never empty.
        let (k_lo, k_hi) = if BANDED {
            (
                qlen.saturating_sub(left_reference.len() + w),
                qlen.min(right_reference.len() + w),
            )
        } else {
            (0, qlen)
        };
        debug_assert!(
            k_lo <= k_hi,
            "no jump point is in band: k_lo={k_lo} > k_hi={k_hi} at w={w}, qlen={qlen}"
        );
        let (mut best, mut best_diff, mut best_k) = (i64::MIN, i64::MAX, k_lo);
        for k in k_lo..=k_hi {
            let f = self.row_max_b[k] as i64;
            let g = self.row_max[qlen - k] as i64;
            let total = f + g;
            let diff = (f - g).abs();
            if total > best || (total == best && diff < best_diff) {
                best = total;
                best_diff = diff;
                best_k = k;
            }
        }
        let k = best_k;

        // --- the right arm: query[..k] from right_reference[0], stopping wherever it scored ---
        //
        // Pass A's triple, explicitly - this reads before the `tb`/`tb_base` swap below, so the
        // `_b` slots are the ones holding A.
        let rj = if BANDED {
            arg_for_row_banded(
                &self.nm_b,
                &self.tb_base_b,
                self.row0_arg_b,
                k,
                qlen,
                right_reference.len(),
                w,
            )
        } else {
            arg_for_row(
                &self.nm_b,
                &self.tb_base_b,
                self.row0_arg_b,
                k,
                right_reference.len(),
            )
        };
        std::mem::swap(&mut self.tb, &mut self.tb_b);
        std::mem::swap(&mut self.tb_base, &mut self.tb_base_b);
        let right = AlignmentResult {
            score: self.row_max_b[k],
            query_start: 0,
            query_end: k,
            ref_start: 0,
            ref_end: rj,
            cigar: self.replay::<false, BANDED>(query, right_reference, k, rj, w),
        };
        std::mem::swap(&mut self.tb, &mut self.tb_b);
        std::mem::swap(&mut self.tb_base, &mut self.tb_base_b);

        // --- the left arm: query[k..] ending at the end of left_reference ---
        //
        // The replay runs in the REVERSED frame, walking from (qlen-k, lj) back to the origin.
        // Row `i` there is `query[qlen-i]`, so starting at `qlen-k` covers exactly `query[k..]`,
        // and the CIGAR comes out in forward order without a flip.
        let llen = left_reference.len();
        // Pass B's triple: the swaps above have been undone, so the primary slots hold B.
        let lj = if BANDED {
            arg_for_row_banded(
                &self.nm,
                &self.tb_base,
                self.row0_arg,
                qlen - k,
                qlen,
                llen,
                w,
            )
        } else {
            arg_for_row(&self.nm, &self.tb_base, self.row0_arg, qlen - k, llen)
        };
        let left = AlignmentResult {
            score: self.row_max[qlen - k],
            query_start: k,
            query_end: qlen,
            ref_start: llen - lj,
            ref_end: llen,
            cigar: self.replay::<true, BANDED>(query, left_reference, qlen - k, lj, w),
        };

        SplitReferenceAlignment {
            score: left.score + right.score,
            left,
            right,
        }
    }

    /// One `TRACK_ROW_MAX` pass, leaving `row_max`/`nm`/`tb` populated.
    ///
    /// `ARG_LARGEST = false`: the row maxima pin the **smallest** `j`.
    ///
    /// `w` is read only when `BANDED`. Under a band **only rows `0..=min(qlen, rlen + w)` are
    /// filled** - row `i`'s in-band columns are `[i - w, i + w]` intersected with `[0, rlen]`,
    /// which is empty past that - so every caller must bound its scan over `row_max` by that
    /// same `last_row`. Reading past it reads whatever the previous call left. This is the same
    /// rule `local` states at its own reduce, and it is the caller's to keep because `row_pass`
    /// does not know which rows its caller intends to ask about.
    #[allow(clippy::too_many_arguments)] // a kernel entry point: the scheme plus the bandwidth
    fn row_pass<const REVERSED: bool, const BANDED: bool>(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        match_score: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
        w: usize,
    ) {
        let (qlen, rlen) = (query.len(), refseq.len());
        if rlen == 0 {
            // No reference to align against: every row's best is its own column-0 cell, i.e.
            // spanning that much query as one pure insertion. Column 0 is load-bearing here -
            // for a query segment neither reference explains, it is the ONLY candidate.
            if self.row_max.len() < qlen + 1 + LANES {
                self.row_max.resize(qlen + 1 + LANES, 0);
            }
            for i in 0..=qlen {
                self.row_max[i] = if i == 0 {
                    0
                } else {
                    -(gap_open as i32) - (i as i32 - 1) * gap_extend as i32
                };
            }
            // Every row's best is its column-0 cell, so every argmax is 0 - including row 0's,
            // which is carried scalar and would otherwise be whatever the previous call left.
            self.row0_arg = 0;
            self.tb_base.clear();
            return;
        }
        self.prepare::<REVERSED>(query, refseq, match_score, gap_open, gap_extend);
        self.size_flags::<true, BANDED>(qlen, rlen, w);
        unsafe {
            self.run::<true, false, true, false, BANDED>(
                qlen,
                rlen,
                match_score,
                mismatch,
                gap_open,
                gap_extend,
                w,
            )
        };
    }

    /// The kernel behind every mode above.
    ///
    /// `TRACK_LAST_ROW` argmaxes `H[qlen][j]` instead of pinning the end cell at
    /// `(qlen, rlen)`; `REVERSED` runs the whole thing back-to-front. **Both are `const`
    /// generics, never runtime flags**, so that adding a mode costs the others nothing and each
    /// mode monomorphizes to exactly the code it needs.
    fn align<const TRACK_LAST_ROW: bool, const REVERSED: bool>(
        &mut self,
        query: &[u8],
        refseq: &[u8],
        match_score: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
    ) -> AlignmentResult {
        let (qlen, rlen) = (query.len(), refseq.len());

        // A zero-length side never reaches the kernel: there is no anti-diagonal to walk.
        if qlen == 0 || rlen == 0 {
            return self.degenerate::<TRACK_LAST_ROW, REVERSED>(qlen, rlen, gap_open, gap_extend);
        }

        self.prepare::<REVERSED>(query, refseq, match_score, gap_open, gap_extend);
        self.size_flags::<false, false>(qlen, rlen, 0);

        let (score, end_j) = unsafe {
            self.run::<true, TRACK_LAST_ROW, false, false, false>(
                qlen,
                rlen,
                match_score,
                mismatch,
                gap_open,
                gap_extend,
                0,
            )
        };
        let cigar = self.replay::<REVERSED, false>(query, refseq, qlen, end_j, 0);
        let (ref_start, ref_end) = ref_span::<REVERSED>(rlen, end_j);

        AlignmentResult {
            score,
            query_start: 0,
            query_end: qlen,
            ref_start,
            ref_end,
            cigar,
        }
    }

    /// The `qlen == 0` / `rlen == 0` cases, which skip the kernel entirely.
    ///
    /// Not merely a fast path - the modes genuinely differ here. With an empty query, `global`
    /// must still consume the reference (`rlen` deletions), while `local_reference_end` stops
    /// at `j = 0` and returns the empty alignment, because the reference tail is free.
    fn degenerate<const TRACK_LAST_ROW: bool, const REVERSED: bool>(
        &mut self,
        qlen: usize,
        rlen: usize,
        gap_open: u8,
        gap_extend: u8,
    ) -> AlignmentResult {
        let gap = |k: usize| -> i32 {
            if k == 0 {
                0
            } else {
                -(gap_open as i32) - (k as i32 - 1) * gap_extend as i32
            }
        };
        let mut cigar = Cigar::default();
        if qlen > 0 {
            // The query must always be spanned in full, in every mode here.
            cigar.push(CigarOperation::Insertion, qlen);
            let (ref_start, ref_end) = ref_span::<REVERSED>(rlen, 0);
            return AlignmentResult {
                score: gap(qlen),
                query_start: 0,
                query_end: qlen,
                ref_start,
                ref_end,
                cigar,
            };
        }
        // qlen == 0. Global must still consume the reference; a free reference end need not.
        let consumed = if TRACK_LAST_ROW { 0 } else { rlen };
        if consumed > 0 {
            cigar.push(CigarOperation::Deletion, consumed);
        }
        let (ref_start, ref_end) = ref_span::<REVERSED>(rlen, consumed);
        AlignmentResult {
            score: gap(consumed),
            query_start: 0,
            query_end: 0,
            ref_start,
            ref_end,
            cigar,
        }
    }

    /// Where cell `(i, j)`'s flags live, as `(word base, lane bit)`, with `p = i + j`.
    ///
    /// The address, not the flags: the replay's `H` state reads one word and sometimes two,
    /// and the `E`/`F` states read exactly one, so returning all four packed would build bits
    /// nobody looks at. The lane arithmetic is paid once per cell, the loads once per question,
    /// and the four words stay adjacent so they still share a cache line. Most of the replay's
    /// cost is the lane arithmetic and the branches rather than the loads, so this is a small
    /// win on a serial walk - see `replay`.
    #[inline(always)]
    fn flag_at<const BANDED: bool>(
        &self,
        p: usize,
        i: usize,
        qlen: usize,
        rlen: usize,
        w: usize,
    ) -> (*const u32, u32) {
        // `saturating_sub`, and it is load-bearing: `p - rlen` underflows for every
        // anti-diagonal before the band starts sliding right, which sends `lane` negative and
        // the unchecked read below off the end. Caught by a segfault, not by a wrong answer.
        //
        // The lane is `i - lo(p)`, so it moves with the band: the fill starts each
        // anti-diagonal's chunks at the band's top row, not the matrix's. Reading a banded
        // fill's flags with the unbanded `lo` would silently address the wrong cell.
        let (i_from, _) = band_rows::<BANDED>(p, qlen, rlen, w);
        let lane = i - i_from;
        debug_assert!(
            i >= i_from && (self.tb_base[p] as usize + lane / LANES) * TB_WORDS < self.tb.len(),
            "flag_at({p}, {i}) is outside the region the fill wrote"
        );
        // SAFETY: `p` is a live anti-diagonal and `i` one of its interior rows, so `base` is
        // within the region the fill wrote. The replay only ever asks about cells it walked to,
        // and the same argument covers `tb_base[p]`: `p <= qlen + rlen`, which is what
        // `size_flags` sized it to.
        unsafe {
            let base = (*self.tb_base.get_unchecked(p) as usize + lane / LANES) * TB_WORDS;
            (self.tb.as_ptr().add(base), 1u32 << (lane % LANES))
        }
    }

    /// Walk the flags back from `(qlen, rlen)` to `(0, 0)` and build the CIGAR.
    ///
    /// Scalar, and unavoidably so - a traceback is a serial walk down one path. Only the flag
    /// *generation* vectorises.
    ///
    /// # Diagonal-first priority minimises the edit distance
    ///
    /// Many alignments tie for optimal, and the score does not choose between them - this does.
    /// Among equally-optimal paths a caller wants the **fewest edits**: `10X` rather than
    /// `10I10D` at the same score, because ten substituted bases are what happened and a
    /// twenty-base structural event is not.
    ///
    /// Checking the diagonal first is what delivers that. It is **exact for strongly affine
    /// schemes** such as the default, where `gap_open` is large enough against `gap_extend` that
    /// the gappy alternative rarely ties the clean one at all. It is a heuristic outside them: a
    /// greedy backward walk commits to a tied step knowing nothing about what it costs
    /// downstream, so for schemes whose `gap_open` approaches `gap_extend` - where ties are
    /// dense - it can occasionally return an optimal alignment that is not the cleanest one.
    ///
    /// Exactness in general would need the lexicographic `(score, -edits)` DP, and the edit
    /// count is the one quantity here with **no scheme-only bound**: its differences grow with
    /// sequence length, where every value in the difference recurrence is bounded by the scheme
    /// alone. It would need an i32 lane and its own per-row state including anti-diagonal `p-2`,
    /// which this kernel never keeps.
    ///
    /// # The trap, and why `from_diag` is a stored bit
    ///
    /// Diagonal-first **cannot** be spelled "if neither gap achieved `A_G`'s max, go diagonal":
    /// that fallthrough is unreachable on precisely the ties this exists to resolve. It needs
    /// its own bit, and that bit needs the clamp correction - see [`TB_FROM_DIAG`].
    ///
    /// # Under a band the replay stays inside it, and does not check that it does
    ///
    /// Every step is decided by a flag, and the fill is what makes the flags safe to follow:
    /// an out-of-band predecessor is sealed to `0`, the smallest value a lane can hold, so it
    /// loses `A_G`'s max and its direction is never the one recorded. Two of the four
    /// directions need an argument beyond "it lost the max", because a sealed `0` can *tie* an
    /// `A_G` that is itself `0`:
    ///
    /// - **The band's lower edge is safe for free.** There `A_G = max(s_G, ΔE'_G, 0)`. If
    ///   `A_G == 0` then `ΔE'_G == 0` too, so `from_e` is set and the replay leaves via E
    ///   before it can fall through to F. If `A_G > 0` then either E achieved it (`from_e`
    ///   set) or `s_G` did, and an `s_G` that is strictly positive is one that did not clamp,
    ///   so `from_diag` is set and trustworthy. One of the two always fires, and both lead
    ///   inwards.
    /// - **The band's upper edge needs one bit cleared.** There the fallthrough F leads
    ///   *inwards* and E leads out, so the fill clears `from_e` on that lane - see the mask in
    ///   `run`. Discarding it costs nothing: it can only ever have been the seal tying `A_G`.
    ///
    /// So `(i, j)` is always in band, and `|i - j| <= w` with `1 <= i, j` is exactly the
    /// condition under which `flag_at` addresses a cell the fill wrote. The two closed-form
    /// exits are covered too: the replay can only reach `(i, 0)` or `(0, j)` in band, which
    /// bounds `i` and `j` by `w`, and the straight run to the origin it emits there is in band
    /// the whole way.
    fn replay<const REVERSED: bool, const BANDED: bool>(
        &self,
        query: &[u8],
        refseq: &[u8],
        start_i: usize,
        start_j: usize,
        w: usize,
    ) -> Cigar {
        let (qlen, rlen) = (query.len(), refseq.len());
        let mut cigar = Cigar::with_capacity(16);

        enum St {
            H,
            E,
            F,
        }
        let (mut i, mut j) = (start_i, start_j);
        let mut st = St::H;

        // The current run, held in registers and flushed only when the operation changes.
        //
        // `Cigar::push` already merges a repeated op into the run in progress - but it is the
        // **host's** type and reaches it through the `Vec`: two loads, a compare and a store,
        // per cell, to add one to a counter. Doing that here makes it a register compare and an
        // add, and `push` is then called once per *run* rather than once per *cell*. The
        // emitted CIGAR is identical (`push` merges either way; consecutive runs here differ by
        // construction, so it only ever takes its append path).
        //
        // `run_len == 0` means "nothing pending", which is why `run_op`'s initial value is
        // never read - guarded by that check, not by an `Option`, so the hot path has no
        // discriminant to test.
        let mut run_op = CigarOperation::Eq;
        let mut run_len = 0usize;
        macro_rules! emit {
            ($op:expr, $n:expr) => {{
                let op = $op;
                if run_len != 0 && op == run_op {
                    run_len += $n;
                } else {
                    if run_len != 0 {
                        cigar.push(run_op, run_len);
                    }
                    run_op = op;
                    run_len = $n;
                }
            }};
        }

        while i > 0 || j > 0 {
            // The boundaries carry no flags and need none: down column 0 the only path is all
            // insertions, along row 0 all deletions, whatever state we arrived in. And it is
            // the *whole* remaining path, so emit it in closed form and stop - the loop used to
            // walk it one base at a time to reach the same answer.
            if j == 0 {
                emit!(CigarOperation::Insertion, i);
                break;
            }
            if i == 0 {
                emit!(CigarOperation::Deletion, j);
                break;
            }
            let p = i + j;
            // SAFETY: `(w, bit)` addresses cell (i, j)'s flag words, which the fill wrote; see
            // `flag_at`. Each `read` below is one of the four adjacent words at that address.
            debug_assert!(
                !BANDED || i.abs_diff(j) <= w,
                "the replay stepped out of the band of {w}, to ({i}, {j})"
            );
            let (word, bit) = self.flag_at::<BANDED>(p, i, qlen, rlen, w);
            let flag = |which: usize| unsafe { (*word.add(which) & bit) != 0 };
            match st {
                St::H => {
                    // Diagonal FIRST: on a tie it spends one edit where a gap pair would spend
                    // two or more. See the tie-break note above.
                    //
                    // `from_f` is the fallthrough and is never stored: `A_G` is a max of three
                    // operands, so if neither the diagonal nor E achieved it, `ΔF'_G` did.
                    if flag(TB_FROM_DIAG) {
                        // `=` vs `X` is **literal byte identity**, re-derived from the original
                        // bytes rather than from how the kernel scored the cell. `U` scores as
                        // `T` but reports `X` against one; an unknown byte never matches, yet
                        // reports `=` against a byte-identical copy of itself. See `mod.rs`,
                        // "Alphabet" - this decoupling is deliberate.
                        // `=` vs `X` in the ORIGINAL frame. Under REVERSED the kernel ran on
                        // back-to-front codes, so cell (i, j) is the pair
                        // (query[qlen-i], refseq[rlen-j]) in the caller's bytes.
                        let (qb, rb) = if REVERSED {
                            (query[qlen - i], refseq[rlen - j])
                        } else {
                            (query[i - 1], refseq[j - 1])
                        };
                        let op = if qb == rb {
                            CigarOperation::Eq
                        } else {
                            CigarOperation::X
                        };
                        emit!(op, 1);
                        i -= 1;
                        j -= 1;
                    } else if flag(TB_FROM_E) {
                        // S[i,j] <- E[i-1,j]: the step consumes query base i.
                        emit!(CigarOperation::Insertion, 1);
                        i -= 1;
                        st = St::E;
                    } else {
                        emit!(CigarOperation::Deletion, 1);
                        j -= 1;
                        st = St::F;
                    }
                }
                St::E => {
                    if flag(TB_E_OPEN) {
                        // The run opened here: drop onto the score matrix without moving.
                        st = St::H;
                    } else {
                        emit!(CigarOperation::Insertion, 1);
                        i -= 1;
                    }
                }
                St::F => {
                    if flag(TB_F_OPEN) {
                        st = St::H;
                    } else {
                        emit!(CigarOperation::Deletion, 1);
                        j -= 1;
                    }
                }
            }
        }
        // The last run has no following op to flush it. An empty alignment (start_i == start_j
        // == 0, which `local_end` returns for a hopeless extension) never enters the loop and
        // leaves `run_len == 0`, so this correctly emits nothing.
        if run_len != 0 {
            cigar.push(run_op, run_len);
        }

        // The replay walks backwards and so builds the CIGAR back-to-front, needing one flip.
        // In the REVERSED frame, back-to-front IS forward order - so the flip is skipped and
        // the two reversals annihilate. This is what makes the `_start` modes free.
        if !REVERSED {
            cigar.reverse();
        }
        cigar
    }

    #[target_feature(enable = "avx2")]
    /// # `TRACK_ROW_MAX`: per-cell absolute scores, exactly, with nothing to bound
    ///
    /// `local_end`, `local_start` and `split_reference` argmax over **absolute** `H[i][j]`,
    /// which a difference kernel does not have. The paper's answer is a three-level
    /// decomposition whose 8-bit small delta **its own authors cannot bound** ("too complex for
    /// rigorous analysis"). Under a contract that reads "exact, always", that is unusable.
    ///
    /// We do not need it. `S[i,j] - S[i,j-1] = ΔV[i,j]` *is* the definition of ΔV, and on
    /// anti-diagonal `p` row `i` sits at column `j = p - i` while on `p-1` it sits at `j-1`. So
    /// `abs[i] += ΔV[i, p-i]` - one add per cell, from a value the fill already stores. Exact,
    /// per-lane (no scan, no cross-lane dependency), seeded from `S[i,0]` in closed form.
    ///
    /// The accumulator is `i32` and **that is not a fallback to a wider kernel**: the DP stays
    /// u8 at 32 lanes and `fits_u8_lanes` still has no length term. It is a reconstruction
    /// riding alongside. Its own bound is `(qlen + rlen) * max(match, gap_open) < i32::MAX`,
    /// i.e. `qlen + rlen < ~8.4M` at the worst possible u8 scores - and **memory binds first by
    /// an enormous margin**: an 8.4M-base alignment needs ~35 TB of traceback flags. There is
    /// no reachable input where this is the limit.
    ///
    /// The seed also disposes of a hazard: the partial chunk tail writes garbage to rows
    /// *above* the live range, and those rows **are** read later. But garbage only lands on row
    /// `i` at `p <= i`, and the boundary seed at `p = i` scrubs it - so still no mask, only the
    /// write ordering already required.
    ///
    /// `ARG_LARGEST` picks the tie-break on `j`: `local_end` wants the **largest** (longest
    /// extension), `split_reference`'s row maxima want the **smallest**.
    ///
    /// # `BANDED`: the strip, and the two rows either side of it that pay for it
    ///
    /// `BANDED` narrows each anti-diagonal's interior to `|i - j| <= w` (see [`band_rows`]).
    /// The chunk loop does not change at all - it is handed a shorter range - and that is the
    /// whole runtime win. What the band costs is at its **edges**, where the recurrence reads
    /// cells the fill never visited, and it costs three things:
    ///
    /// 1. **Two sealed bytes per anti-diagonal.** `A_G = max(s_G, ΔE'_G[i-1,j], ΔF'_G[i,j-1])`
    ///    reaches one row above the strip's top and one row below its bottom. Both are sealed
    ///    to `0`, the smallest value a `u8` lane can hold, so they lose every max they enter -
    ///    which is exactly what "unreachable" has to mean here. **This is not a sentinel**: it
    ///    is the bottom of the range every delta already lives in, so nothing about the
    ///    `[0, upper]` invariant, the plain `sub_epi8`s, or `fits_u8_lanes` moves an inch.
    ///
    /// 2. **Two garbage deltas per anti-diagonal, which are provably dead.** The *other* two
    ///    reads at each edge - `ΔV_G[i-1,j]` at the top, `ΔH_G[i,j-1]` at the bottom - are
    ///    subtrahends, and a seal cannot make a difference correct. So they are left as
    ///    garbage, and the band's geometry retires them: each bound moves at half of `p`'s
    ///    rate, so an edge that was flat at `p` **must** advance at `p+1`, and advancing is
    ///    precisely what makes the cell that consumed the garbage fall outside the next
    ///    anti-diagonal's reads. The values are computed and never read.
    ///
    /// 3. **One seed per row, replacing the column-0 one.** `abs[i]` accumulates `ΔV` along
    ///    row `i`, so it needs a starting absolute. Unbanded that is `H[i][0]`, seeded at the
    ///    matrix corner. Under a band, row `i > w` has no column-0 cell - it is `i` off the
    ///    diagonal - and enters at column `i - w` instead, where `ΔV`'s left neighbour is
    ///    outside the strip and its own `ΔV` is therefore the garbage of (2). But `ΔH` at that
    ///    same cell is **not** garbage: `ΔH` reads *upwards*, and the cell above is in band. So
    ///    `abs[i] = abs[i-1] + ΔH_G[i, i-w] - gap_open` seeds the row exactly, from a value the
    ///    fill already computed, one anti-diagonal's `abs` and one scalar add.
    ///
    /// That seed doubles as the scrub, and has to: the partial chunk's tail writes garbage to
    /// rows *below* the strip, and under a band those rows are reached far later than the
    /// unbanded corner seed would have scrubbed them. The row is seeded at the moment it enters
    /// the band, which is after the last garbage write to it and before its first real read.
    /// # Which mode compiles which flags
    ///
    /// Every public mode reaches this function through one of three wrappers. `REVERSED` is
    /// threaded separately (see `prepare`) and never reaches `run`.
    ///
    /// | mode | wrapper | TRACE | LAST_ROW | ROW_MAX | ARG_LARGEST | BANDED | score read from |
    /// |---|---|---|---|---|---|---|---|
    /// | `global`, exact | `align` | y | - | - | - | - | `anchor` |
    /// | `global`, banded | `row_pass` | y | - | y | - | y | `abs[qlen]` |
    /// | `local_reference_*`, exact | `align` | y | y | - | - | - | `best_score` |
    /// | `local_reference_*`, banded | `row_pass` | y | - | y | - | y | `row_max[qlen]` |
    /// | `local_end` / `local_start` | `local` | y | - | y | y | opt | reduce over `row_max` |
    /// | `split_reference` (two passes) | `row_pass` | y | - | y | - | opt | `row_max_b[k] + row_max[qlen-k]` |
    ///
    /// **`TRACK_LAST_ROW` and `BANDED` are never compiled together**: `align` is the only user
    /// of the former and is never banded, because its `anchor` walk cannot survive a band (see
    /// [`U8Probe::local_reference`]).
    ///
    #[allow(clippy::too_many_arguments)] // the scheme is five scalars, plus the bandwidth
    unsafe fn run<
        const TRACE: bool,
        const TRACK_LAST_ROW: bool,
        const TRACK_ROW_MAX: bool,
        const ARG_LARGEST: bool,
        const BANDED: bool,
    >(
        &mut self,
        qlen: usize,
        rlen: usize,
        match_score: u8,
        mismatch: u8,
        gap_open: u8,
        gap_extend: u8,
        w: usize,
    ) -> (i32, usize) {
        unsafe {
            let m = match_score as i64;
            let mm = mismatch as i64;
            let go = gap_open as i64 - gap_extend as i64;
            let ge = gap_extend as i64;
            let c = 2 * (go + ge); // == 2 * gap_open
            let upper = (m + c) as u8;

            // The two values of s_G, broadcast once. Only the mismatch side needs the clamp:
            // s_G(match) = M + C is exactly `upper`, which `fits_u8_lanes` already bounds and
            // which cannot be negative. s_G(mismatch) = C - mismatch **can** be negative (any
            // scheme with mismatch > 2*gap_open), and clamping it to zero is lossless because
            // A_G only ever maxes it against ΔE'_G/ΔF'_G, which are >= 0. Wrapping instead
            // would make the worst substitution score the largest value in the lane, win every
            // max, and return a score BETTER than optimal - the signature of a sentinel bug.
            let s_g_match = upper;
            let s_g_mismatch = (c - mm).max(0) as u8;
            let go_u8 = go as u8;

            // ΔV = ΔV_G - gap_open, un-offset, widened to the accumulator's lane type.
            let gap_open_vec32 = _mm256_set1_epi32(gap_open as i32);

            let match_vec = _mm256_set1_epi8(s_g_match as i8);
            let mismatch_vec = _mm256_set1_epi8(s_g_mismatch as i8);
            let go_vec = _mm256_set1_epi8(go_u8 as i8);
            // All-ones exactly when `s_G(mismatch)` did NOT clamp, i.e. when `cmpeq(A_G, s_G)`
            // means what it says at every cell. Loop-invariant, and a whole-vector constant, so
            // the `or` it feeds costs one uop on any port and sits off the critical path
            // (`a_g -> de_next/df_next -> store`). See TB_FROM_DIAG.
            let unclamped_vec = if c - mm >= 0 {
                _mm256_set1_epi8(-1)
            } else {
                _mm256_setzero_si256()
            };

            let q_ptr = self.q_codes.as_ptr();
            let r_ptr = self.r_rev.as_ptr();

            // Hoist every buffer to a raw pointer, ONCE, and ping-pong the pointers rather
            // than the `Vec`s.
            //
            // `mem::swap(&mut self.prev, &mut self.cur)` would swap two `Deltas`, i.e. eight
            // `Vec` headers per anti-diagonal, which dominates at small `n` where the ramp runs
            // more anti-diagonals than chunks of actual DP. Four pointer swaps are register
            // moves.
            let mut p_dh = self.prev.dh.as_mut_ptr();
            let mut p_dv = self.prev.dv.as_mut_ptr();
            let mut p_de = self.prev.de.as_mut_ptr();
            let mut p_df = self.prev.df.as_mut_ptr();
            let mut c_dh = self.cur.dh.as_mut_ptr();
            let mut c_dv = self.cur.dv.as_mut_ptr();
            let mut c_de = self.cur.de.as_mut_ptr();
            let mut c_df = self.cur.df.as_mut_ptr();
            let tb_ptr = self.tb.as_mut_ptr();
            let tb_base_ptr = self.tb_base.as_mut_ptr();
            let abs_ptr = self.abs.as_mut_ptr();
            let row_max_ptr = self.row_max.as_mut_ptr();
            let nm_ptr = self.nm.as_mut_ptr();

            // The running absolute value at the anti-diagonal's bottom row: the paper's 64-bit
            // large offset, reduced to what a global score actually needs - one i64, not a
            // per-lane middle delta.
            let mut anchor: i64 = 0;
            // `local_reference_end`'s objective, for free.
            //
            // `anchor` holds the absolute score at the anti-diagonal's *bottom* row,
            // hi(p) = min(qlen, p). So for every p >= qlen it is exactly S[qlen][p - qlen] -
            // **it already walks the last row**, one cell per anti-diagonal, which is precisely
            // the row this mode maximizes. No per-cell absolute values, no extra vector work:
            // one scalar compare per anti-diagonal, and only for p >= qlen.
            //
            // Ties resolve to the smallest `ref_end`, which falls out of updating on a strict
            // improvement only, since p (and therefore j) increases monotonically.
            let mut best_score = i64::MIN;
            let mut best_j = 0usize;

            // Row 0, in closed form.
            //
            // `H[0][j]` is "consume j reference bases, align no query": 0 at j=0, and
            // `-(gap_open + (j-1)*gap_extend) <= 0` beyond it. Both gap costs are unsigned, so
            // **row 0's max is always 0**, attained at j=0, and the only question is which `j`
            // the tie-break names:
            //
            //   * `!ARG_LARGEST` wants the smallest `j`, and j=0 attains the max, so: 0.
            //   * `gap_open > 0` makes every j>=1 strictly negative, so: 0.
            //   * `gap_open == 0 < gap_extend` leaves `H[0][1] == 0` tying with j=0 and
            //     everything past it negative, so: 1.
            //   * both zero ties the whole row at 0, so the largest j is `rlen` - and that is
            //     not a curiosity, it is the free-gap case `local_end` has a test for.
            //
            // `abs[0]` is deliberately not seeded: nothing reads it. The group loop starts at
            // `i_from`, which is `>= 1` always, so row 0 has no lane and never enters it.
            // The last case is the one the band touches: with both gap costs zero the whole of
            // row 0 ties at 0, and the largest tying `j` is the furthest cell the band still
            // reaches, `min(rlen, w)`, not the end of the reference.
            let row0_arg: i32 = if !ARG_LARGEST || gap_open > 0 {
                0
            } else if gap_extend > 0 {
                1.min(rlen) as i32
            } else if BANDED {
                rlen.min(w) as i32
            } else {
                rlen as i32
            };
            if TRACK_ROW_MAX {
                *row_max_ptr.add(0) = 0;
            }
            // The running chunk counter, so `tb_base` is filled by the fill itself rather
            // than by a separate O(qlen+rlen) pass.
            //
            // The flag address is recomputed per chunk rather than carried as a bumped
            // pointer: a live `*mut u32` across the chunk loop stops LLVM proving it does not
            // alias the delta buffers, so it stops keeping their pointers in registers.
            let mut chunk: usize = 0;

            // The interior row range, maintained INCREMENTALLY rather than recomputed.
            //
            // `i_from = max(1, p - rlen)` and `i_to = min(qlen, p - 1)` are three scalar ops
            // per anti-diagonal, charged on nearly every iteration of work at small `n`.
            //
            // Both are monotone step functions of `p`: `i_to` climbs to `qlen` and stops;
            // `i_from` sits at 1 until `p > rlen + 1` and then climbs. So carry them and
            // advance by one, with two perfectly-predicted branches (always-taken, then
            // never-taken) instead of three ALU ops on the dependency path.
            let mut i_from = 1usize;
            let mut i_to = 0usize;

            // The band's own two bounds, carried the same way and for the same reason.
            //
            // `b_lo = ceil((p - w) / 2)` and `b_hi = floor((p + w) / 2)`. Each advances by one
            // on alternating anti-diagonals - and, since their parities are opposite, **exactly
            // one of the two advances every step**. So the pair costs one branch on a parity
            // bit, not two divisions, on a path where the whole anti-diagonal may be only three
            // chunks of real work.
            //
            // These are `band_rows`' two band terms evaluated at `p = 1`, written the way they
            // are *defined* rather than the way they are shortest: clippy would have `(1 + w) /
            // 2` as `w.div_ceil(2)`, which is the same number and hides that this is
            // `floor((p + w) / 2)` at `p = 1`. The `debug_assert` against `band_rows` in the
            // loop is what checks the pair, and it can only do that if they stay recognisable.
            #[allow(clippy::manual_div_ceil)]
            let (mut b_lo, mut b_hi) = if BANDED {
                (1usize.saturating_sub(w).div_ceil(2), (1 + w) / 2)
            } else {
                (0, usize::MAX)
            };
            // `(lo, hi)` at `p - 1`, which is what says whether each edge just advanced - and
            // therefore whether the cell beyond it is one the last anti-diagonal wrote or one
            // the seal has to cover. Seeded at `p = 0`'s empty range.
            let (mut prev_lo, mut prev_hi) = (1usize, 0usize);

            for p in 1..=(qlen + rlen) {
                debug_assert_eq!(i_from, 1.max(p.saturating_sub(rlen)));
                debug_assert_eq!(i_to, qlen.min(p - 1));

                // The interior, band included. Under `!BANDED` these fold to `i_from`/`i_to`
                // and every branch guarded by `BANDED` below disappears with them.
                let lo = if BANDED { i_from.max(b_lo) } else { i_from };
                let hi = if BANDED { i_to.min(b_hi) } else { i_to };
                debug_assert_eq!((lo, hi), band_rows::<BANDED>(p, qlen, rlen, w));

                // Is this anti-diagonal's top (bottom) row sitting exactly ON the band's edge,
                // so that the cell beyond it is outside the strip?
                //
                // Row `i` sits at `i - j = 2i - p`, so the top row is at `+w` exactly when it is
                // the band's own bound *and* `p - w` is even - at odd `p - w` the bound rounds
                // inwards and the row above is still in the strip. `lo == b_lo` is the "band,
                // not the matrix corner, is what is limiting here" half; without it a top row
                // pinned by `max(1, p - rlen)` would be sealed against a neighbour that is
                // perfectly real.
                //
                // The two share a parity, and that is not a coincidence: the bounds move at
                // half of `p`'s rate on opposite phases, so the strip is `w + 1` rows wide with
                // both edges live, then `w` rows wide with neither. Every edge fix below fires
                // on the wide ones and none on the narrow ones.
                //
                // `live` is not pedantry. An anti-diagonal with no interior still *has* band
                // bounds, and they can sit anywhere - the strip runs off the end of the query
                // long before `p` does whenever `qlen` is short against the band. Without this
                // the edge fixes below would write a row nothing computed and, worse, `|=` a
                // flag word at `chunk_base` that belongs to the *next* anti-diagonal, because an
                // empty one advances `chunk` by zero.
                let live = lo <= hi;
                let edge_phase = BANDED && live && (p + w).is_multiple_of(2);
                let top_sealed = edge_phase && lo == b_lo;
                let bottom_sealed = edge_phase && hi == b_hi;
                debug_assert!(
                    !(top_sealed && bottom_sealed) || lo < hi,
                    "both edges live means a strip {} rows wide, so they cannot be the same row",
                    w + 1
                );
                // Once the strip leaves the matrix it never comes back, so a live anti-diagonal
                // always has a live predecessor to read from. `p <= 2` is the one exception, and
                // it is the ordinary one: `p = 1` has no interior in any mode, banded or not,
                // and `p = 2`'s reads land on the boundary cells rather than on `p = 1`'s
                // interior.
                debug_assert!(
                    !BANDED || !live || p <= 2 || prev_lo <= prev_hi,
                    "anti-diagonal {p} is live but {} was empty: its reads land on rows that \
                     nothing wrote",
                    p - 1
                );

                if TRACE {
                    *tb_base_ptr.add(p) = chunk as u32;
                }
                // `chunk` is advanced past this anti-diagonal before the TRACK_ROW_MAX pass
                // runs, so that pass cannot use it directly - it needs THIS anti-diagonal's
                // base, which is what the flags were written against.
                let chunk_base = chunk;

                let mut i0 = lo;
                while i0 <= hi {
                    // Two loads from anti-diagonal p-1, one byte apart. This is the trick that
                    // replaces the paper's shiftl/shiftr - the addressing does it.
                    let dv_up = _mm256_loadu_si256(p_dv.add(i0 - 1) as *const __m256i);
                    let de_up = _mm256_loadu_si256(p_de.add(i0 - 1) as *const __m256i);
                    let dh_left = _mm256_loadu_si256(p_dh.add(i0) as *const __m256i);
                    let df_left = _mm256_loadu_si256(p_df.add(i0) as *const __m256i);

                    // `rlen + i0 - p`, never `rlen - p + i0`: the index is provably
                    // non-negative on every live lane, but `rlen - p` alone underflows the
                    // moment the anti-diagonal stops growing and starts sliding right. The
                    // grouping is load-bearing.
                    let qv = _mm256_loadu_si256(q_ptr.add(i0 - 1) as *const __m256i);
                    let rv = _mm256_loadu_si256(r_ptr.add(rlen + i0 - p) as *const __m256i);
                    // Hoisted out of the blend because the traceback needs it too: on a clamped
                    // scheme `from_diag` is only trustworthy at match cells. See TB_FROM_DIAG.
                    let eq_qr = _mm256_cmpeq_epi8(qv, rv);
                    let s_g = _mm256_blendv_epi8(mismatch_vec, match_vec, eq_qr);

                    // Critical path 4: two maxes with nothing feeding them. ΔE'_G and ΔF'_G
                    // have already absorbed ΔV and ΔH, so A_G reads them directly rather than
                    // adding first. That is the entire reason Eq 3 exists.
                    let a_g = _mm256_max_epu8(_mm256_max_epu8(s_g, de_up), df_left);

                    let de_next = _mm256_max_epu8(a_g, _mm256_add_epi8(de_up, go_vec));
                    let df_next = _mm256_max_epu8(a_g, _mm256_add_epi8(df_left, go_vec));

                    // The traceback flags: equalities on values already in registers, extracted
                    // a bit per lane. `const TRACE` so the score-only path monomorphizes without
                    // them entirely.
                    if TRACE {
                        let tb = tb_ptr.add((chunk + (i0 - lo) / LANES) * TB_WORDS);
                        // `from_diag`/`from_e`: which operands achieved A_G's max. `from_f` is
                        // NOT stored - it is the fallthrough, recovered by the replay from "at
                        // least one of the three achieved the max". See TB_WORDS for why the
                        // diagonal and not a gap is the one that must be stored.
                        //
                        // `unclamped_vec` is all-ones for every scheme that does not clamp
                        // `s_G` (the default among them), which makes the `or` a no-op on the
                        // mask and the `and` an identity - so the correction is paid for only
                        // where it is needed, without a second monomorphization of this loop.
                        *tb.add(TB_FROM_DIAG) = _mm256_movemask_epi8(_mm256_and_si256(
                            _mm256_cmpeq_epi8(a_g, s_g),
                            _mm256_or_si256(eq_qr, unclamped_vec),
                        )) as u32;
                        *tb.add(TB_FROM_E) =
                            _mm256_movemask_epi8(_mm256_cmpeq_epi8(a_g, de_up)) as u32;
                        // `e_open`/`f_open`: whether the gap state opened from A_G rather than
                        // extending. Phrased against the max the fill already computed, so it
                        // is an equality and needs no unsigned compare (AVX2 has none for u8).
                        *tb.add(TB_E_OPEN) =
                            _mm256_movemask_epi8(_mm256_cmpeq_epi8(de_next, a_g)) as u32;
                        *tb.add(TB_F_OPEN) =
                            _mm256_movemask_epi8(_mm256_cmpeq_epi8(df_next, a_g)) as u32;
                    }

                    // Plain subs, no saturation: every result is provably in [0, upper], so
                    // the minuend is always >= the subtrahend. No sentinel exists to cross.
                    _mm256_storeu_si256(c_dh.add(i0) as *mut __m256i, _mm256_sub_epi8(a_g, dv_up));
                    _mm256_storeu_si256(
                        c_dv.add(i0) as *mut __m256i,
                        _mm256_sub_epi8(a_g, dh_left),
                    );
                    _mm256_storeu_si256(
                        c_de.add(i0) as *mut __m256i,
                        _mm256_sub_epi8(de_next, dh_left),
                    );
                    _mm256_storeu_si256(
                        c_df.add(i0) as *mut __m256i,
                        _mm256_sub_epi8(df_next, dv_up),
                    );

                    i0 += LANES;
                }
                if TRACE && lo <= hi {
                    chunk += (hi - lo) / LANES + 1;
                }

                // Does row `hi` enter the band here, at its own first in-band column rather
                // than at the matrix corner? Then its `abs` has to be seeded rather than
                // accumulated - see (3) in this function's docs.
                //
                // `hi > prev_hi` says the row is new; `hi > w` says it is too far off the
                // diagonal to have had a column-0 cell, which is what separates this from the
                // corner seed below. The two are exclusive by construction and cover every row
                // exactly once.
                // `bottom_sealed` IS "row `hi` enters the band here": the band's lower bound
                // reaches row `hi` for the first time exactly at `p = 2*hi - w`, which is the
                // anti-diagonal where it is the bound and the parity is even. Rows within `w` of
                // the origin enter at the matrix corner instead and are seeded there, and
                // `bottom_sealed` provably never fires for them - it needs `hi = (p + w) / 2`
                // and `hi <= p - 1`, which together force `hi > w`.
                let seeds_row = bottom_sealed && TRACK_ROW_MAX;
                debug_assert!(
                    !seeds_row || hi > w,
                    "a corner-seeded row cannot also band-enter"
                );
                // `abs[hi - 1]` **before** the group loop overwrites it: it must be row
                // `hi - 1`'s absolute at anti-diagonal `p - 1`, i.e. `S[hi-1, p-hi]`, which is
                // the cell diagonally up-left of the one being seeded. One anti-diagonal later
                // and it is the wrong cell.
                //
                // `hi > w >= 1` puts `hi - 1` at 1 or above, so this never reaches row 0 -
                // which has no `abs` and needs none.
                let seed_carry = if seeds_row { *abs_ptr.add(hi - 1) } else { 0 };

                // The per-row absolutes and their running maxima, for the interior of this
                // anti-diagonal. Deliberately a second pass over the same rows rather than
                // work folded into the chunk loop: it needs `cur.dv`, which the chunk loop is
                // still writing, and it is 8-lane i32 work against the loop's 32-lane u8 -
                // interleaving them would spill registers in the hot loop for the benefit of
                // the three modes that use this, and cost the three that do not.
                if TRACK_ROW_MAX {
                    // One 8-lane i32 group. A macro rather than a closure so both call sites
                    // inline identically - `&mut self` plus intrinsics makes a closure here
                    // fight the borrow checker for nothing.
                    macro_rules! group {
                        ($base:expr) => {{
                            let base = $base;
                            // 8 u8 deltas -> 8 i32, un-offset to true ΔV.
                            let dv = _mm256_cvtepu8_epi32(_mm_loadl_epi64(
                                c_dv.add(base) as *const __m128i
                            ));
                            let a = _mm256_add_epi32(
                                _mm256_loadu_si256(abs_ptr.add(base) as *const __m256i),
                                _mm256_sub_epi32(dv, gap_open_vec32),
                            );
                            _mm256_storeu_si256(abs_ptr.add(base) as *mut __m256i, a);

                            // No `j` or `i` is built here: the bit stream records *that* a cell
                            // took the max, not which, and `arg_for_row` reconstructs the
                            // coordinates from the bit's `(p, lane)` position afterwards.
                            let m = _mm256_loadu_si256(row_max_ptr.add(base) as *const __m256i);

                            // STAY OFF PORT 5. `blendv_epi8` and the `cvtepu8_epi32` above both
                            // issue only on the shuffle port, and this block saturates it long
                            // before it runs out of instructions. So `max_epi32` takes the max
                            // with no blend at all, and `cmpeq(new, a)` recovers the mask - both
                            // on p0/p1. A blend is fewer instructions here and slower.
                            let new_m = _mm256_max_epi32(m, a);
                            // `cmpeq(max(m,a), a)` is `a >= m` - which IS the largest-j rule.
                            // The smallest-j rule needs a strict `a > m`.
                            let take = if ARG_LARGEST {
                                _mm256_cmpeq_epi32(new_m, a)
                            } else {
                                _mm256_cmpgt_epi32(a, m)
                            };
                            _mm256_storeu_si256(
                                row_max_ptr.add(base) as *mut __m256i,
                                if ARG_LARGEST {
                                    new_m
                                } else {
                                    _mm256_max_epi32(m, a)
                                },
                            );
                            // One bit per lane instead of a `j` per lane, and the group
                            // *returns* its 8 bits rather than storing them - the caller packs
                            // all four into one word and stores once per 32 cells. See the
                            // store site for why the store count, not the store width, is what
                            // had to come down.
                            _mm256_movemask_ps(_mm256_castsi256_ps(take)) as u32
                        }};
                    }

                    // A CONSTANT trip count, so the four groups unroll. Bounding it to skip the
                    // dead lanes of a partial chunk takes the unroll with it and costs more than
                    // the lanes are worth.
                    let mut i0 = lo;
                    while i0 <= hi {
                        // Four groups, four 8-bit masks, ONE 32-bit store for the chunk.
                        //
                        // One full word per chunk, not four narrow stores: `nm` is written once
                        // and never re-read, so it streams, and four 4-byte streaming stores per
                        // chunk miss D1 badly where one 32-bit store does not.
                        let m0 = group!(i0);
                        let m1 = group!(i0 + 8);
                        let m2 = group!(i0 + 16);
                        let m3 = group!(i0 + 24);
                        *nm_ptr.add(chunk_base + (i0 - lo) / LANES) =
                            m0 | (m1 << 8) | (m2 << 16) | (m3 << 24);
                        i0 += LANES;
                    }
                }

                // The band's seals and its one seed, AFTER both loops and BEFORE the boundary
                // cells. Both halves of that ordering are load-bearing:
                //
                //  * after the chunk loop, for the reason the boundary cells are - a 32-lane
                //    store is unconditional, so the partial last chunk scribbles over the row
                //    just past the interior, which is exactly the row the bottom seal owns;
                //  * before the boundary cells, because at the matrix corners the seals aim at
                //    rows that ARE the boundary (`hi + 1 == p` while the anti-diagonal is still
                //    growing, `lo - 1 == 0` before it lifts off). The boundary is the truth
                //    there; letting it land second is cheaper than teaching the seals to duck.
                if BANDED && live {
                    // --- the seals: the two cells the recurrence reaches past the strip ---
                    //
                    // `live` is what keeps these in bounds, and it is not a formality: once the
                    // strip slides off the bottom of the query, `lo` keeps climbing with `p`
                    // while the buffers stop at `qlen + 1 + LANES`, so `c_de[lo - 1]` walks
                    // straight off the end of the allocation. Nothing reads those anti-diagonals
                    // - the strip never re-enters the matrix, which is what the `debug_assert`
                    // above pins - so there is nothing to seal. (Found as heap corruption, not
                    // as a wrong answer.)
                    //
                    // Given `live`: `1 <= lo <= hi <= qlen`, so `lo - 1` and `hi + 1` are both
                    // inside every delta buffer. Where they land on a boundary cell that
                    // legitimately exists - row 0 while the band still reaches it, row `p` while
                    // the anti-diagonal is still growing - the boundary block below runs second
                    // and takes it back.
                    //
                    // `p = 1` has no interior in any mode and needs no seals: everything `p = 2`
                    // reads is a boundary cell, which the block below writes.
                    *c_de.add(lo - 1) = 0;
                    *c_df.add(hi + 1) = 0;

                    // --- the edges: where a seal of `0` is not enough ---
                    //
                    // A seal makes the *three-way* max fall the right way: `0` is the smallest a
                    // lane can hold, so it loses `A_G = max(s_G, ΔE'_G, ΔF'_G)` and the
                    // out-of-band direction is never taken.
                    //
                    // **It does not survive the gap recurrence.** `ΔE'_G[i,j] = max(A_G,
                    // ΔE'_G[i-1,j] + go)` adds `go` to the seal *before* the max, so a sealed
                    // `0` arrives as `go` - and `go` beats `A_G` whenever `A_G < go`, which is
                    // most cells under a strongly affine scheme. The result: a `ΔE'_G` too large
                    // by `go - A_G`, and an `e_open` bit that says the gap extends from a cell
                    // outside the band. The score comes out wrong, and the replay walks up out
                    // of the strip and reads flags the fill never wrote.
                    //
                    // Both are fixed by saying what the band means: at the strip's top edge
                    // `E[i-1,j]` is unreachable, so the E-run **must open here**, and
                    // `ΔE'_G[i,j] = A_G[i,j] - ΔH_G[i,j-1]` - which is exactly what the chunk
                    // loop just stored in `ΔV_G[i,j]`. So the correction is a byte copy, not an
                    // arithmetic redo, and the F edge mirrors it: `ΔF'_G[i,j] = A_G[i,j] -
                    // ΔV_G[i-1,j]`, which is `ΔH_G[i,j]`.
                    //
                    // Each edge copies from the delta whose own inputs are in band - `ΔV_G` reads
                    // leftwards and the top edge's left neighbour is inside; `ΔH_G` reads upwards
                    // and the bottom edge's upper neighbour is inside. That is why this works at
                    // all, and why the copies do not cross: `top_sealed && bottom_sealed` implies
                    // a strip `w + 1 >= 2` rows wide, so `lo != hi`.
                    if top_sealed {
                        *c_de.add(lo) = *c_dv.add(lo);
                        if TRACE {
                            // In `St::H`: "H came from E[i-1,j]" - out of band, so clear it. It
                            // could only ever have been the seal tying an `A_G` of 0. The replay
                            // then falls through to F, which at this edge leads *inwards*.
                            *tb_ptr.add(chunk_base * TB_WORDS + TB_FROM_E) &= !1u32;
                            // In `St::E`: "the run opened here" - it must have.
                            *tb_ptr.add(chunk_base * TB_WORDS + TB_E_OPEN) |= 1u32;
                        }
                    }
                    if bottom_sealed {
                        *c_df.add(hi) = *c_dh.add(hi);
                        if TRACE {
                            // In `St::F`: "the run opened here". `from_f` itself needs no mask -
                            // it is the replay's fallthrough, and at this edge one of the two
                            // flags it falls through *from* is always set. See `replay`.
                            let lane = hi - lo;
                            *tb_ptr.add((chunk_base + lane / LANES) * TB_WORDS + TB_F_OPEN) |=
                                1u32 << (lane % LANES);
                        }
                    }

                    // Row `hi` enters the band at column `hi - w`: seed its running absolute
                    // from the cell above-left, which is in band, via a `ΔH` the chunk loop just
                    // computed from in-band inputs. See (3) in this function's docs.
                    //
                    // `row_max` is *assigned*, not maxed: this is the row's first candidate, and
                    // whatever is sitting there is the partial chunk's garbage.
                    if seeds_row {
                        let s = seed_carry + *c_dh.add(hi) as i32 - gap_open as i32;
                        *abs_ptr.add(hi) = s;
                        *row_max_ptr.add(hi) = s;
                        // No `nm` bit is written for it. `arg_for_row_banded` returns this
                        // column when its scan finds nothing above, which is precisely when this
                        // cell is still the row's argmax - so the bit's two possible values give
                        // the same answer. See that function.
                    }
                }

                // Boundary cells, AFTER the chunk loop - see the hazard note in the module
                // docs. The values are Eq 3's initial conditions.
                // Under a band the boundary cells are **conditional**, and this is where the
                // band's edge and the matrix's edge finally meet. `H[0][j]` sits `j` off the
                // diagonal and `H[i][0]` sits `i` off it, so the boundary leaves the strip after
                // `w` steps and past that these cells **do not exist**. Writing them anyway
                // would hand the recurrence a real, reachable-looking predecessor outside the
                // band - which is the one thing the seals are there to prevent, undone one line
                // later by the code that runs after them.
                //
                // Past `w` the seals above are the truth, and they already cover the only two
                // rows that get read: `c_de[lo - 1]` reaches row 0 while `lo` is still 1, and
                // `c_df[hi + 1]` reaches row `p` while the anti-diagonal is still growing.
                // Everything further out is never read at all - by the time `p > w + 1` the
                // band no longer reaches either boundary.
                if p <= rlen && (!BANDED || p <= w) {
                    let b = if p == 1 { 0 } else { go_u8 };
                    *c_dv.add(0) = b;
                    *c_de.add(0) = b;
                    // Row 0's max and argmax are hoisted out of the loop - see the closed form
                    // above it. `abs[0]` goes with them: it is never read, because `i_from >= 1`
                    // means the group loop never touches row 0.
                }
                if p <= qlen && (!BANDED || p <= w) {
                    let b = if p == 1 { 0 } else { go_u8 };
                    *c_dh.add(p) = b;
                    *c_df.add(p) = b;
                    // The corner seed, gated with the boundary it reads. Rows past `w` are
                    // seeded at their real first column instead, by `seeds_row` above; between
                    // the two, every row is seeded exactly once.
                    if TRACK_ROW_MAX {
                        // Row p enters the matrix here, at its column-0 cell. Seeding from
                        // S[p,0] is load-bearing twice over: H[i][0] - "span query[..i] as one
                        // insertion, consume no reference" - is a real argmax candidate and for
                        // a query segment nothing explains it is the only one; and this write
                        // is what scrubs the garbage the partial chunk tail left on this row at
                        // earlier anti-diagonals.
                        // No arg store: `j = 0` is what `arg_for_row` returns when no interior
                        // cell of the row ever took the max, so the column-0 seed is implicit.
                        let s_i0 = -(gap_open as i32) - (p as i32 - 1) * gap_extend as i32;
                        *abs_ptr.add(p) = s_i0;
                        *row_max_ptr.add(p) = s_i0;
                    }
                }

                // While the anti-diagonal still grows downwards the bottom cell is (p, 0), a
                // boundary known in closed form. Once it slides right the bottom row is pinned
                // at qlen and each step adds one ΔV, un-offset. At p = qlen + rlen that lands
                // exactly on S[qlen][rlen].
                //
                // `anchor` is dead in the `TRACK_ROW_MAX` modes - `local` and `row_pass`
                // discard `run`'s return value - but it needs no gate: `run` is private, every
                // call site is visible under LTO, and the unused return value is enough for LLVM
                // to kill the whole chain including the `c_dv[qlen]` load.
                anchor = if p <= qlen {
                    if p == 1 {
                        -(gap_open as i64)
                    } else {
                        anchor - ge
                    }
                } else {
                    anchor + *c_dv.add(qlen) as i64 - (go + ge)
                };

                if TRACK_LAST_ROW && p >= qlen && anchor > best_score {
                    best_score = anchor;
                    best_j = p - qlen;
                }

                // Four register moves, not eight Vec headers.
                std::mem::swap(&mut p_dh, &mut c_dh);
                std::mem::swap(&mut p_dv, &mut c_dv);
                std::mem::swap(&mut p_de, &mut c_de);
                std::mem::swap(&mut p_df, &mut c_df);

                // Advance the row range for p+1. AFTER the body, not before: the first cut
                // incremented on entry and so read min(qlen, p) where it wanted min(qlen, p-1)
                // - caught immediately by the debug_assert above, which is why it is there.
                if i_to < qlen {
                    i_to += 1;
                }
                if p >= rlen + 1 {
                    i_from += 1;
                }
                if BANDED {
                    (prev_lo, prev_hi) = (lo, hi);
                    // Exactly one of the two band edges advances per anti-diagonal, and which
                    // one is decided by a single parity bit: `b_hi = floor((p + w) / 2)` steps
                    // when `p + w` becomes even, `b_lo = ceil((p - w) / 2)` when it becomes odd.
                    // This is that test for `p + 1`.
                    //
                    // `p + 1 > w` is the saturation, carried rather than recomputed. Until the
                    // band's lower edge lifts off row 0 - `p <= w`, where `ceil((p - w) / 2)`
                    // is negative and `band_rows`' `saturating_sub` clamps it - it must not
                    // move at all. Without this the edge climbs from the very first
                    // anti-diagonal and the strip walks off the top of the matrix; the
                    // `debug_assert` against the closed form above is what says so.
                    if (p + 1 + w).is_multiple_of(2) {
                        b_hi += 1;
                    } else if p + 1 > w {
                        b_lo += 1;
                    }
                }
            }

            let _ = upper;
            // The bound `size_flags` computed, against the chunks this fill actually wrote.
            // `chunk` is one past the last one written, so `<=` is the exact condition. Every
            // store above is unchecked, so a bound that is too tight is heap corruption rather
            // than a wrong score - and the one time this crate got a band bound wrong, that is
            // exactly how it showed up.
            //
            // Read the lengths here rather than hoisting them beside `tb_ptr`: this chunk loop
            // is register-starved (see the note beside `chunk`), and down here the whole thing
            // is inside the `debug_assert!` so a release build reads nothing at all.
            debug_assert!(
                !TRACE || chunk * TB_WORDS <= self.tb.len(),
                "the fill wrote {chunk} chunks but tb holds {} ({TB_WORDS} words per chunk): \
                 the chunk bound is too tight, and every store past it was heap corruption",
                self.tb.len() / TB_WORDS,
            );
            debug_assert!(
                !(TRACE && TRACK_ROW_MAX) || chunk * NM_WORDS <= self.nm.len(),
                "the fill wrote {chunk} chunks but nm holds {} - see the tb assert above",
                self.nm.len() / NM_WORDS
            );
            if TRACK_ROW_MAX {
                self.row0_arg = row0_arg;
            }
            if TRACK_LAST_ROW {
                (best_score as i32, best_j)
            } else {
                (anchor as i32, rlen)
            }
        }
    }
}
