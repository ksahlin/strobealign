// This was extremely simplified from xxhash-rust

const PRIME_1: u64 = 0x9E3779B185EBCA87;
const PRIME_2: u64 = 0xC2B2AE3D27D4EB4F;
const PRIME_3: u64 = 0x165667B19E3779F9;
const PRIME_4: u64 = 0x85EBCA77C2B2AE63;
const PRIME_5: u64 = 0x27D4EB2F165667C5;

/// xxh64, but it can only be used for a single u64
#[inline]
pub fn xxh64(input: u64) -> u64 {
    let mut result= PRIME_5.wrapping_add(8);
    result ^= input.wrapping_mul(PRIME_2).rotate_left(31).wrapping_mul(PRIME_1);
    result = result.rotate_left(27).wrapping_mul(PRIME_1).wrapping_add(PRIME_4);
    result ^= result >> 33;
    result = result.wrapping_mul(PRIME_2);
    result ^= result >> 29;
    result = result.wrapping_mul(PRIME_3);
    result ^= result >> 32;
    result
}

#[test]
fn test() {
    // "hello123"
    assert_eq!(xxh64(0x3332316f6c6c6568u64), 0x2d119fced5ee39ce);
}
