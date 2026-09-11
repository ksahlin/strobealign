//! The "modes" that strobealign knows about:
//!
//! - Mapping-only mode (`-x`)
//! - Abundance estimation mode (`--aemb`)
//! - Extension alignment mode (default)

use std::cmp::min;

use crate::chain::Chain;

pub mod abundance;
pub mod align;
pub mod map;

/// Return mapping quality for the top chain
pub fn mapping_quality(chains: &[Chain]) -> u8 {
    if chains.len() <= 1 {
        return 60;
    }
    let s1 = chains[0].score;
    let s2 = chains[1].score;
    // from minimap2: MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    let min_matches = min(chains[0].anchors.len(), 10) as f32 / 10.0;
    let uncapped_mapq = 40.0 * (1.0 - s2 / s1) * min_matches * s1.ln();

    uncapped_mapq.min(60.0) as u8
}
