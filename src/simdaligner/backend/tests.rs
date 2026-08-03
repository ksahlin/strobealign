//! Conformance tests for the backends.
//!
//! Two layers, and the split is deliberate because of what can and cannot be run where.
//!
//! [`ops`] checks each backend's **21 operations** against a plain-Rust definition of what each
//! one means, on whatever machine the tests run on. This is the layer that validates a port: run
//! the suite on an aarch64 box and it checks the NEON transcription instruction by instruction,
//! including the synthesized movemask, which is the one operation with real logic in it.
//!
//! [`differential`] checks that two backends produce **byte-identical alignments** end to end.
//! On x86_64 that compares AVX2 against SSE4.1 and the portable kernel (32 lanes against 16),
//! which is what exercises the lane-width bookkeeping: chunk counts, the `nm` packing, and the
//! band's edge cases. NEON and simd128 are the same 16-lane shape, so a green run here says the
//! shared half is right for them too, and leaves only their intrinsics to [`ops`].
//!
//! Neither layer can substitute for the other, and neither runs NEON on x86. Together they mean
//! a NEON port is wrong only if a NEON intrinsic does not do what its name says.

use super::*;

/// A deterministic xorshift, so a failure reproduces exactly.
///
/// Deliberately not `fastrand`: a test that fails once every few hundred seeds and cannot be
/// replayed is worse than no test.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn byte(&mut self) -> u8 {
        (self.next() >> 24) as u8
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() >> 16) as usize % n
    }

    fn dna(&mut self, len: usize) -> Vec<u8> {
        // Includes bytes outside ACGT on purpose: the unknown-base path is part of the kernel
        // and has its own scoring rule.
        const ALPHABET: &[u8] = b"ACGTACGTACGTNU";
        (0..len)
            .map(|_| ALPHABET[self.below(ALPHABET.len())])
            .collect()
    }
}

mod ops {
    use super::*;

    /// Read a `V8` back out as bytes.
    unsafe fn bytes<B: Backend>(v: B::V8) -> Vec<u8> {
        let mut out = vec![0u8; B::LANES];
        unsafe { B::store8(out.as_mut_ptr(), v) };
        out
    }

    /// Read a `V32` back out as `i32`s.
    unsafe fn words<B: Backend>(v: B::V32) -> Vec<i32> {
        let mut out = vec![0i32; B::LANES32];
        unsafe { B::store32(out.as_mut_ptr(), v) };
        out
    }

    /// Every `u8` operation, against its plain-Rust definition.
    fn check_u8<B: Backend>() {
        let mut rng = Rng(0x2545F491_4F6CDD1D);
        for round in 0..64 {
            // Round 0 pins the boundary values the random rounds would only hit by luck: 0x00
            // and 0xFF matter because they are what comparisons produce, and 0x80 matters
            // because it is the MSB that `movemask8` and `blend8` read.
            let (a, b): (Vec<u8>, Vec<u8>) = if round == 0 {
                let pattern = [0u8, 0xFF, 0x80, 0x7F, 1, 254, 128, 127];
                (
                    (0..B::LANES).map(|i| pattern[i % pattern.len()]).collect(),
                    (0..B::LANES)
                        .map(|i| pattern[(i + 3) % pattern.len()])
                        .collect(),
                )
            } else {
                // Full-range bytes, with `b` copying every third lane from `a` so that `eq8`
                // still sees both outcomes in every vector.
                //
                // The collisions used to come from masking both to `& 0x0F` instead. That also
                // forced every MSB to zero, which made `movemask8` return 0 on every random
                // round and left the bit-order contract (the one thing NEON has to synthesize)
                // tested only by round 0, whose pattern happens to be a palindrome and so
                // survives a reversed movemask. Caught by mutating `<< i` to `<< (15 - i)`.
                let a: Vec<u8> = (0..B::LANES).map(|_| rng.byte()).collect();
                let b: Vec<u8> = (0..B::LANES)
                    .enumerate()
                    .map(|(i, _)| if i % 3 == 0 { a[i] } else { rng.byte() })
                    .collect();
                (a, b)
            };

            unsafe {
                let va = B::load8(a.as_ptr());
                let vb = B::load8(b.as_ptr());

                assert_eq!(bytes::<B>(va), a, "{}: load8/store8 round trip", B::NAME);

                let splat = bytes::<B>(B::splat8(0xAB));
                assert_eq!(splat, vec![0xAB; B::LANES], "{}: splat8", B::NAME);

                let expect: Vec<u8> = (0..B::LANES).map(|i| a[i].max(b[i])).collect();
                assert_eq!(
                    bytes::<B>(B::max8(va, vb)),
                    expect,
                    "{}: max8 (unsigned)",
                    B::NAME
                );

                let expect: Vec<u8> = (0..B::LANES).map(|i| a[i].wrapping_add(b[i])).collect();
                assert_eq!(
                    bytes::<B>(B::add8(va, vb)),
                    expect,
                    "{}: add8 wraps",
                    B::NAME
                );

                let expect: Vec<u8> = (0..B::LANES).map(|i| a[i].wrapping_sub(b[i])).collect();
                assert_eq!(
                    bytes::<B>(B::sub8(va, vb)),
                    expect,
                    "{}: sub8 wraps",
                    B::NAME
                );

                let expect: Vec<u8> = (0..B::LANES)
                    .map(|i| if a[i] == b[i] { 0xFF } else { 0x00 })
                    .collect();
                let eq = B::eq8(va, vb);
                assert_eq!(
                    bytes::<B>(eq),
                    expect,
                    "{}: eq8 saturates the mask",
                    B::NAME
                );

                let expect: Vec<u8> = (0..B::LANES).map(|i| a[i] & b[i]).collect();
                assert_eq!(bytes::<B>(B::and8(va, vb)), expect, "{}: and8", B::NAME);

                let expect: Vec<u8> = (0..B::LANES).map(|i| a[i] | b[i]).collect();
                assert_eq!(bytes::<B>(B::or8(va, vb)), expect, "{}: or8", B::NAME);

                // Selection under a saturated mask, which is the only kind the kernel passes.
                let expect: Vec<u8> = (0..B::LANES)
                    .map(|i| if a[i] == b[i] { a[i] } else { b[i] })
                    .collect();
                assert_eq!(
                    bytes::<B>(B::blend8(vb, va, eq)),
                    expect,
                    "{}: blend8 picks if_true where the mask is set",
                    B::NAME
                );

                // The operation with actual logic on NEON. Lane i must land in bit i.
                let mut expect_mask = 0u32;
                for (i, lane) in a.iter().enumerate() {
                    expect_mask |= ((lane >> 7) as u32) << i;
                }
                assert_eq!(
                    B::movemask8(va),
                    expect_mask,
                    "{}: movemask8 must take each lane's MSB into bit i (round {round})",
                    B::NAME
                );
            }
        }
    }

    /// Every `i32` operation, against its plain-Rust definition.
    fn check_i32<B: Backend>() {
        let mut rng = Rng(0x9E3779B9_7F4A7C15);
        for round in 0..64 {
            let (a, b): (Vec<i32>, Vec<i32>) = if round == 0 {
                // Negatives matter: these are absolute scores, and `max32`/`gt32` are signed.
                let pattern = [0i32, -1, i32::MIN, i32::MAX, -5, 5, -100000, 100000];
                (
                    (0..B::LANES32)
                        .map(|i| pattern[i % pattern.len()])
                        .collect(),
                    (0..B::LANES32)
                        .map(|i| pattern[(i + 2) % pattern.len()])
                        .collect(),
                )
            } else {
                let a: Vec<i32> = (0..B::LANES32)
                    .map(|_| (rng.next() as i32) % 1000 - 500)
                    .collect();
                let b: Vec<i32> = (0..B::LANES32)
                    .enumerate()
                    .map(|(i, _)| {
                        if i % 3 == 0 {
                            a[i]
                        } else {
                            (rng.next() as i32) % 1000 - 500
                        }
                    })
                    .collect();
                (a, b)
            };

            unsafe {
                let va = B::load32(a.as_ptr());
                let vb = B::load32(b.as_ptr());

                assert_eq!(words::<B>(va), a, "{}: load32/store32 round trip", B::NAME);
                assert_eq!(
                    words::<B>(B::splat32(-77)),
                    vec![-77; B::LANES32],
                    "{}: splat32",
                    B::NAME
                );

                let expect: Vec<i32> = (0..B::LANES32).map(|i| a[i].wrapping_add(b[i])).collect();
                assert_eq!(words::<B>(B::add32(va, vb)), expect, "{}: add32", B::NAME);

                let expect: Vec<i32> = (0..B::LANES32).map(|i| a[i].wrapping_sub(b[i])).collect();
                assert_eq!(words::<B>(B::sub32(va, vb)), expect, "{}: sub32", B::NAME);

                let expect: Vec<i32> = (0..B::LANES32).map(|i| a[i].max(b[i])).collect();
                assert_eq!(
                    words::<B>(B::max32(va, vb)),
                    expect,
                    "{}: max32 must be SIGNED",
                    B::NAME
                );

                let expect: Vec<i32> = (0..B::LANES32)
                    .map(|i| if a[i] == b[i] { -1 } else { 0 })
                    .collect();
                assert_eq!(words::<B>(B::eq32(va, vb)), expect, "{}: eq32", B::NAME);

                let expect: Vec<i32> = (0..B::LANES32)
                    .map(|i| if a[i] > b[i] { -1 } else { 0 })
                    .collect();
                assert_eq!(
                    words::<B>(B::gt32(va, vb)),
                    expect,
                    "{}: gt32 must be SIGNED",
                    B::NAME
                );

                let mut expect_mask = 0u32;
                for (i, lane) in a.iter().enumerate() {
                    expect_mask |= (((*lane as u32) >> 31) & 1) << i;
                }
                assert_eq!(
                    B::movemask32(va),
                    expect_mask,
                    "{}: movemask32 must take each lane's sign bit into bit i (round {round})",
                    B::NAME
                );

                // Widening load: `LANES32` bytes, zero-extended.
                let src: Vec<u8> = (0..B::LANES32 + 8).map(|_| rng.byte()).collect();
                let expect: Vec<i32> = (0..B::LANES32).map(|i| src[i] as i32).collect();
                assert_eq!(
                    words::<B>(B::widen8to32(src.as_ptr())),
                    expect,
                    "{}: widen8to32 must ZERO-extend",
                    B::NAME
                );
            }
        }
    }

    fn check<B: Backend>() {
        assert_eq!(
            B::LANES,
            4 * B::LANES32,
            "{}: the nm packing needs LANES == 4 * LANES32",
            B::NAME
        );
        assert!(B::LANES <= 32, "{}: a lane per bit of a u32", B::NAME);
        check_u8::<B>();
        check_i32::<B>();
    }

    #[test]
    fn scalar_backend_matches_its_definition() {
        check::<Scalar>();
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_backend_matches_its_definition() {
        if is_x86_feature_detected!("avx2") {
            check::<Avx2>();
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn sse41_backend_matches_its_definition() {
        if is_x86_feature_detected!("sse4.1") {
            check::<Sse41>();
        }
    }

    #[cfg(all(
        target_arch = "aarch64",
        target_feature = "neon",
        target_endian = "little"
    ))]
    #[test]
    fn neon_backend_matches_its_definition() {
        check::<Neon>();
    }

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    #[test]
    fn simd128_backend_matches_its_definition() {
        check::<Simd128>();
    }
}

mod differential {
    use super::*;
    use crate::simdaligner::kernel::U8Probe;

    /// Run every alignment mode on both backends and require identical results.
    ///
    /// Identical, not merely equally-scoring: the tie-break rules are part of the contract, so
    /// two backends that pick different equally-optimal CIGARs are a bug, not a wash.
    fn compare<A: Backend, B: Backend>() {
        let mut rng = Rng(0xDEADBEEF_CAFEF00D);
        let mut pa = U8Probe::<A>::new();
        let mut pb = U8Probe::<B>::new();

        // Lengths that straddle a chunk boundary at both widths (16 and 32), so the partial
        // last chunk, which is where the write-ordering hazards live, is exercised at every
        // offset rather than only at the ones a round number happens to hit.
        let lengths = [
            1usize, 2, 3, 15, 16, 17, 31, 32, 33, 47, 63, 64, 65, 100, 129, 200,
        ];

        for &qlen in &lengths {
            for &rlen in &lengths {
                let q = rng.dna(qlen);
                let r = rng.dna(rlen);
                for bandwidth in [None, Some(0), Some(1), Some(8), Some(33)] {
                    let (m, mm, go, ge, eb) = (2u8, 8u8, 12u8, 1u8, 10u32);

                    macro_rules! same {
                        ($mode:ident $(, $extra:expr)*) => {
                            assert_eq!(
                                pa.$mode(&q, &r, m, mm, go, ge $(, $extra)*, bandwidth),
                                pb.$mode(&q, &r, m, mm, go, ge $(, $extra)*, bandwidth),
                                "{} vs {} disagree on {} (qlen={qlen}, rlen={rlen}, band={bandwidth:?})\n\
                                 query={:?}\nref={:?}",
                                A::NAME, B::NAME, stringify!($mode),
                                String::from_utf8_lossy(&q), String::from_utf8_lossy(&r),
                            );
                        };
                    }

                    same!(global_alignment);
                    same!(local_reference_end_alignment);
                    same!(local_reference_start_alignment);
                    same!(local_end_alignment, eb);
                    same!(local_start_alignment, eb);
                }

                // Split-reference takes two references, so it gets its own call.
                let r2 = rng.dna(rlen);
                for bandwidth in [None, Some(8)] {
                    assert_eq!(
                        pa.split_reference_alignment(&q, &r, &r2, 2, 8, 12, 1, bandwidth),
                        pb.split_reference_alignment(&q, &r, &r2, 2, 8, 12, 1, bandwidth),
                        "{} vs {} disagree on split_reference_alignment \
                         (qlen={qlen}, rlen={rlen}, band={bandwidth:?})",
                        A::NAME,
                        B::NAME
                    );
                }
            }
        }
    }

    /// The load-bearing one on x86: 32 lanes against 16, which is what says the lane-width
    /// bookkeeping is right before NEON or simd128 depend on it.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_and_sse41_agree() {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("sse4.1") {
            compare::<Avx2, Sse41>();
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_and_the_portable_kernel_agree() {
        if is_x86_feature_detected!("avx2") {
            compare::<Avx2, Scalar>();
        }
    }

    #[cfg(all(
        target_arch = "aarch64",
        target_feature = "neon",
        target_endian = "little"
    ))]
    #[test]
    fn neon_and_the_portable_kernel_agree() {
        compare::<Neon, Scalar>();
    }

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    #[test]
    fn simd128_and_the_portable_kernel_agree() {
        compare::<Simd128, Scalar>();
    }
}

/// The dispatch itself: that a backend is picked, reported honestly, and aligns correctly
/// whichever one it turned out to be.
mod selection {
    use crate::simdaligner::{Scores, SimdAligner};

    #[test]
    fn the_chosen_backend_is_reported_and_works() {
        let mut aligner = SimdAligner::new(Scores::default());
        let name = aligner.backend_name();

        // `is_supported` means "the pick was vectorized", so it must agree with the name.
        assert_eq!(
            SimdAligner::is_supported(),
            name != "scalar",
            "is_supported disagrees with the backend actually selected ({name})"
        );

        // Whichever kernel was picked, it still has to align. This is the fallback's contract:
        // identical answers, not merely a graceful degradation.
        let r = aligner.global_alignment(b"ACGTACGT", b"ACGTTACGT", None);
        assert_eq!(
            r.cigar.to_string(),
            "3=1D5=",
            "backend {name} produced a wrong CIGAR"
        );
    }

    /// The portable kernel is reachable by name on every machine, whatever was auto-selected,
    /// and agrees with the selected one through the public API rather than at the `U8Probe`
    /// level the `differential` suite works at.
    #[test]
    fn the_portable_kernel_is_reachable_and_agrees() {
        let mut scalar = SimdAligner::with_scalar_kernel(Scores::default());
        assert_eq!(scalar.backend_name(), "scalar");

        let mut selected = SimdAligner::new(Scores::default());
        for (query, reference) in [
            (&b"ACGTACGT"[..], &b"ACGTTACGT"[..]),
            (&b"ACGTACGT"[..], &b"ACGTACGT"[..]),
            (&b"AACCGGTT"[..], &b"AACCTTGGTT"[..]),
            (&b""[..], &b"ACGT"[..]),
        ] {
            let a = scalar.global_alignment(query, reference, None);
            let b = selected.global_alignment(query, reference, None);
            assert_eq!(
                a,
                b,
                "scalar and {} disagree on {:?} vs {:?}",
                selected.backend_name(),
                String::from_utf8_lossy(query),
                String::from_utf8_lossy(reference)
            );
        }
    }
}
