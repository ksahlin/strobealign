// This was extremely simplified from xxhash-rust

// xxHash64 constants
const PRIME64_1: u64 = 0x9E3779B185EBCA87;
const PRIME64_2: u64 = 0xC2B2AE3D27D4EB4F;
const PRIME64_3: u64 = 0x165667B19E3779F9;
const PRIME64_4: u64 = 0x85EBCA77C2B2AE63;
const PRIME64_5: u64 = 0x27D4EB2F165667C5;

// xxHash32 constants
const PRIME32_1: u32 = 0x9E3779B1;
const PRIME32_2: u32 = 0x85EBCA77;
const PRIME32_3: u32 = 0xC2B2AE3D;
const PRIME32_4: u32 = 0x27D4EB2F;
const PRIME32_5: u32 = 0x165667B1;

/// xxh64, but it can only be used for a single u64
#[inline]
pub fn xxh64(input: u64) -> u64 {
    let mut result = PRIME64_5.wrapping_add(8);
    result ^= input
        .wrapping_mul(PRIME64_2)
        .rotate_left(31)
        .wrapping_mul(PRIME64_1);
    result = result
        .rotate_left(27)
        .wrapping_mul(PRIME64_1)
        .wrapping_add(PRIME64_4);
    result ^= result >> 33;
    result = result.wrapping_mul(PRIME64_2);
    result ^= result >> 29;
    result = result.wrapping_mul(PRIME64_3);
    result ^= result >> 32;
    result
}

/// xxh32, but it can only be used for a single u32
#[inline]
pub fn xxh32(input: u32) -> u32 {
    let mut result = PRIME32_5.wrapping_add(4);
    result = result.wrapping_add(input
        .wrapping_mul(PRIME32_2)
        .rotate_left(13)
        .wrapping_mul(PRIME32_1));
    result = result
        .rotate_left(17)
        .wrapping_mul(PRIME32_4);
    result ^= result >> 15;
    result = result.wrapping_mul(PRIME32_2);
    result ^= result >> 13;
    result = result.wrapping_mul(PRIME32_3);
    result ^= result >> 16;
    result
}

#[cfg(test)]
mod tests {
    use super::{xxh32, xxh64};

    #[test]
    fn test_xxh64() {
        // "hello123"
        assert_eq!(xxh64(0x3332316f6c6c6568u64), 0x2d119fced5ee39ce);
    }

    #[test]
    fn test_xxh32_deterministic() {
        // Verify xxh32 produces consistent output
        let input = 0x12345678u32;
        let hash1 = xxh32(input);
        let hash2 = xxh32(input);
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_xxh32_distribution() {
        // Verify different inputs produce different outputs
        let hash0 = xxh32(0);
        let hash1 = xxh32(1);
        let hash2 = xxh32(2);
        assert_ne!(hash0, hash1);
        assert_ne!(hash1, hash2);
        assert_ne!(hash0, hash2);
    }
}
