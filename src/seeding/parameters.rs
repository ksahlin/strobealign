use thiserror::Error;

use super::strobes::{DEFAULT_AUX_LEN, RandstrobeParameters};
use super::syncmers::SyncmerParameters;

/// Pre-defined seeding parameters that work well for a certain
/// "canonical" read length (and similar read lengths)
struct Profile {
    canonical_read_length: usize,
    r_threshold: usize,
    k: usize,
    s_offset: isize,
    w_min: usize,
    w_max: usize,
}

static PROFILES: [Profile; 7] = [
    Profile {
        canonical_read_length: 50,
        r_threshold: 70,
        k: 18,
        s_offset: -4,
        w_min: 1,
        w_max: 4,
    },
    Profile {
        canonical_read_length: 75,
        r_threshold: 90,
        k: 20,
        s_offset: -4,
        w_min: 1,
        w_max: 6,
    },
    Profile {
        canonical_read_length: 100,
        r_threshold: 110,
        k: 20,
        s_offset: -4,
        w_min: 2,
        w_max: 6,
    },
    Profile {
        canonical_read_length: 125,
        r_threshold: 135,
        k: 20,
        s_offset: -4,
        w_min: 3,
        w_max: 8,
    },
    Profile {
        canonical_read_length: 150,
        r_threshold: 175,
        k: 20,
        s_offset: -4,
        w_min: 5,
        w_max: 11,
    },
    Profile {
        canonical_read_length: 250,
        r_threshold: 375,
        k: 22,
        s_offset: -4,
        w_min: 6,
        w_max: 16,
    },
    Profile {
        canonical_read_length: 400,
        r_threshold: usize::MAX,
        k: 23,
        s_offset: -6,
        w_min: 5,
        w_max: 15,
    },
];

/* Settings that influence seeding (creation */
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedingParameters {
    pub canonical_read_length: usize,
    pub syncmer: SyncmerParameters,
    pub randstrobe: RandstrobeParameters,
}

impl SeedingParameters {
    /// Return a parameter-specific filename extension such as ".r100.sti"
    ///
    /// If any of the parameters deviate from the defaults for the current
    /// canonical read length, the returned extension is just ".sti".
    pub fn filename_extension(&self) -> String {
        if *self != SeedingParameters::default_from_read_length(self.canonical_read_length) {
            ".sti".to_string()
        } else {
            format!(".r{}.sti", self.canonical_read_length)
        }
    }
}

#[derive(Error, Debug)]
pub enum InvalidSeedingParameter {
    #[error("Invalid seeding parameter: {0}")]
    InvalidParameter(&'static str),
}

impl SeedingParameters {
    pub fn try_new(
        canonical_read_length: usize,
        k: usize,
        s: usize,
        w_min: usize,
        w_max: usize,
        q: u64,
        max_dist: u8,
        aux_len: u8,
    ) -> Result<Self, InvalidSeedingParameter> {
        if aux_len > 63 {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "aux length must be less than 64",
            ));
        }

        let main_hash_mask = !0u64 << (10 + aux_len);
        Ok(SeedingParameters {
            canonical_read_length,
            syncmer: SyncmerParameters::try_new(k, s)?,
            randstrobe: RandstrobeParameters::try_new(w_min, w_max, q, max_dist, main_hash_mask)?,
        })
    }

    /// Create an IndexParameters instance based on a given read length.
    /// k, s, l, u, c and max_seed_len can be used to override determined parameters
    pub fn from_read_length(
        read_length: usize,
        mut k: Option<usize>,
        mut s: Option<usize>,
        mut w_min: Option<usize>,
        mut w_max: Option<usize>,
        c: Option<u32>,
        max_seed_len: Option<usize>,
        aux_len: u8,
    ) -> Result<SeedingParameters, InvalidSeedingParameter> {
        let default_c = 8;
        let mut canonical_read_length = 50;
        for profile in &PROFILES {
            if read_length <= profile.r_threshold {
                if k.is_none() {
                    k = Some(profile.k);
                }
                if s.is_none() {
                    s = Some((k.unwrap() as isize + profile.s_offset) as usize);
                }
                if w_min.is_none() {
                    w_min = Some(profile.w_min);
                }
                if w_max.is_none() {
                    w_max = Some(profile.w_max);
                }
                canonical_read_length = profile.canonical_read_length;
                break;
            }
        }

        let k = k.unwrap();
        let s = s.unwrap();
        let w_min = w_min.unwrap();
        let w_max = w_max.unwrap();

        let max_dist = match max_seed_len {
            Some(max_seed_len) => {
                if max_seed_len < k || max_seed_len - k > 255 {
                    dbg!(max_seed_len);
                    return Err(InvalidSeedingParameter::InvalidParameter(
                        "max seed length must be between k and k + 255",
                    ));
                }
                (max_seed_len - k) as u8
            }
            None => usize::clamp(canonical_read_length.saturating_sub(70), k, 255) as u8,
        };
        let q = 2u64.pow(c.unwrap_or(default_c)) - 1;

        SeedingParameters::try_new(
            canonical_read_length,
            k,
            s,
            w_min,
            w_max,
            q,
            max_dist,
            aux_len,
        )
    }

    pub fn default_from_read_length(read_length: usize) -> SeedingParameters {
        Self::from_read_length(
            read_length,
            None,
            None,
            None,
            None,
            None,
            None,
            DEFAULT_AUX_LEN,
        )
        .unwrap()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_seeding_parameters() {
        let canonical_read_length = 250;
        let k = 22;
        let s = 18;
        let w_min = 6;
        let w_max = 16;
        let max_dist = 180;
        let q = 255;
        let main_hash_mask = 0xfffffffff8000000;
        let aux_len = 17;
        let sp = SyncmerParameters::try_new(k, s).unwrap();
        let rp = RandstrobeParameters {
            w_min,
            w_max,
            q,
            max_dist,
            main_hash_mask,
        };
        let seeding_parameters = SeedingParameters::try_new(
            canonical_read_length,
            k,
            s,
            w_min,
            w_max,
            q,
            max_dist,
            aux_len,
        )
        .unwrap();
        assert_eq!(
            seeding_parameters.canonical_read_length,
            canonical_read_length
        );
        assert_eq!(seeding_parameters.randstrobe, rp);
        assert_eq!(seeding_parameters.syncmer, sp);

        let ip = SeedingParameters::default_from_read_length(canonical_read_length + 1);
        assert_eq!(ip.canonical_read_length, canonical_read_length);
        assert_eq!(ip.randstrobe, rp);
        assert_eq!(ip.syncmer, sp);
    }

    #[test]
    fn test_seeding_parameters_same_read_length() {
        let sp150a = SeedingParameters::default_from_read_length(150);
        let sp150b = SeedingParameters::default_from_read_length(150);

        assert_eq!(sp150a, sp150b);
    }

    #[test]
    fn test_seeding_parameters_similar_read_length() {
        let sp150 = SeedingParameters::default_from_read_length(150);
        let sp149 = SeedingParameters::default_from_read_length(149);
        let sp151 = SeedingParameters::default_from_read_length(151);

        assert_eq!(sp150, sp149);
        assert_eq!(sp150, sp151);
    }
}
