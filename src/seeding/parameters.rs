use std::fmt::Display;

use thiserror::Error;

use super::strobes::{DEFAULT_AUX_LEN, RandstrobeParameters};
use super::syncmers::SyncmerParameters;

// A preset for seeding parameters (k, s, l, u, etc.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Profile {
    Noisy,
    Canonical50,
    Canonical75,
    Canonical100,
    Canonical125,
    Canonical150,
    Canonical250,
    Canonical400,
}

impl From<usize> for Profile {
    /// Given a read length, returns the corresponding closest canonical `Profile`
    ///
    /// This never returns the `Noisy` variant (construct that explicitly).
    fn from(read_length: usize) -> Self {
        // unwrap is fine because the last r_threshold is usize::MAX
        READ_LENGTH_SETTINGS
            .iter()
            .find(|&profile| read_length <= profile.r_threshold)
            .unwrap()
            .profile
            .clone()
    }
}

impl Profile {
    /// Returns the identifier for the profile (used in the index file name).
    fn identifier(&self) -> String {
        match self {
            Profile::Noisy => "noisy",
            Profile::Canonical50 => "r50",
            Profile::Canonical75 => "r75",
            Profile::Canonical100 => "r100",
            Profile::Canonical125 => "r125",
            Profile::Canonical150 => "r150",
            Profile::Canonical250 => "r250",
            Profile::Canonical400 => "r400",
        }
        .to_string()
    }
}

/// A u32 value is stored in `.sti` files to encode the profile.
impl TryFrom<u32> for Profile {
    type Error = ();

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        match value as u32 {
            x if x == Profile::Noisy as u32 => Ok(Profile::Noisy),
            x if x == Profile::Canonical50 as u32 => Ok(Profile::Canonical50),
            x if x == Profile::Canonical75 as u32 => Ok(Profile::Canonical75),
            x if x == Profile::Canonical100 as u32 => Ok(Profile::Canonical100),
            x if x == Profile::Canonical125 as u32 => Ok(Profile::Canonical125),
            x if x == Profile::Canonical150 as u32 => Ok(Profile::Canonical150),
            x if x == Profile::Canonical250 as u32 => Ok(Profile::Canonical250),
            x if x == Profile::Canonical400 as u32 => Ok(Profile::Canonical400),
            _ => Err(()),
        }
    }
}

impl Display for Profile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Profile::Noisy => write!(f, "noisy"),
            Profile::Canonical50 => write!(f, "canonical 50 nt"),
            Profile::Canonical75 => write!(f, "canonical 75 nt"),
            Profile::Canonical100 => write!(f, "canonical 100 nt"),
            Profile::Canonical125 => write!(f, "canonical 125 nt"),
            Profile::Canonical150 => write!(f, "canonical 150 nt"),
            Profile::Canonical250 => write!(f, "canonical 250 nt"),
            Profile::Canonical400 => write!(f, "canonical 400 nt"),
        }
    }
}

/// Preset seeding parameters based on read length (the "canonical read length").
struct ReadLengthSettings {
    profile: Profile,
    canonical_read_length: usize,
    r_threshold: usize,
    k: usize,
    s_offset: isize,
    w_min: usize,
    w_max: usize,
}

static READ_LENGTH_SETTINGS: [ReadLengthSettings; 7] = [
    ReadLengthSettings {
        profile: Profile::Canonical50,
        canonical_read_length: 50,
        r_threshold: 70,
        k: 18,
        s_offset: -4,
        w_min: 1,
        w_max: 4,
    },
    ReadLengthSettings {
        profile: Profile::Canonical75,
        canonical_read_length: 75,
        r_threshold: 90,
        k: 20,
        s_offset: -4,
        w_min: 1,
        w_max: 6,
    },
    ReadLengthSettings {
        profile: Profile::Canonical100,
        canonical_read_length: 100,
        r_threshold: 110,
        k: 20,
        s_offset: -4,
        w_min: 2,
        w_max: 6,
    },
    ReadLengthSettings {
        profile: Profile::Canonical125,
        canonical_read_length: 125,
        r_threshold: 135,
        k: 20,
        s_offset: -4,
        w_min: 3,
        w_max: 8,
    },
    ReadLengthSettings {
        profile: Profile::Canonical150,
        canonical_read_length: 150,
        r_threshold: 175,
        k: 20,
        s_offset: -4,
        w_min: 5,
        w_max: 11,
    },
    ReadLengthSettings {
        profile: Profile::Canonical250,
        canonical_read_length: 250,
        r_threshold: 375,
        k: 22,
        s_offset: -4,
        w_min: 6,
        w_max: 16,
    },
    ReadLengthSettings {
        profile: Profile::Canonical400,
        canonical_read_length: 400,
        r_threshold: usize::MAX,
        k: 23,
        s_offset: -6,
        w_min: 5,
        w_max: 15,
    },
];

/// Settings that influence seed creation
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedingParameters {
    pub profile: Profile,
    pub syncmer: SyncmerParameters,
    pub randstrobe: RandstrobeParameters,
}

impl From<&Profile> for SeedingParameters {
    fn from(profile: &Profile) -> Self {
        let rp = READ_LENGTH_SETTINGS
            .iter()
            .find(|&rlp| rlp.profile == *profile)
            .expect("missing profile");

        // TODO duplicated code
        let q = 2u64.pow(8) - 1;

        // TODO duplicated code
        let main_hash_mask = !0u64 << (9 + DEFAULT_AUX_LEN);

        // TODO duplicated code
        let max_dist = usize::clamp(rp.canonical_read_length.saturating_sub(70), rp.k, 255) as u8;

        SeedingParameters {
            profile: profile.clone(),
            syncmer: SyncmerParameters::try_new(rp.k, (rp.k as isize + rp.s_offset) as usize)
                .unwrap(),
            randstrobe: RandstrobeParameters::try_new(
                rp.w_min,
                rp.w_max,
                q,
                max_dist,
                main_hash_mask,
            )
            .unwrap(),
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
        profile: Profile,
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

        let main_hash_mask = !0u64 << (9 + aux_len);
        Ok(SeedingParameters {
            profile,
            syncmer: SyncmerParameters::try_new(k, s)?,
            randstrobe: RandstrobeParameters::try_new(w_min, w_max, q, max_dist, main_hash_mask)?,
        })
    }

    /// Create an IndexParameters instance based on a given read length.
    /// k, s, l, u, c and max_seed_len can be used to override determined parameters
    pub fn noisy(
        k: Option<usize>,
        s: Option<usize>,
        w_min: Option<usize>,
        w_max: Option<usize>,
        c: Option<u32>,
        max_seed_len: Option<usize>,
        aux_len: u8,
    ) -> Result<Self, InvalidSeedingParameter> {
        struct P {
            k: usize,
            s_offset: isize,
            w_min: usize,
            w_max: usize,
            max_dist: u8,
        }

        let p = P {
            k: 16,
            s_offset: -4,
            w_min: 2,
            w_max: 2,
            max_dist: 84,
        };

        // TODO this is duplicated in from_read_length
        let k = k.unwrap_or(p.k);
        let s = s.unwrap_or((k as isize + p.s_offset) as usize);
        let w_min = w_min.unwrap_or(p.w_min);
        let w_max = w_max.unwrap_or(p.w_max);
        let q = 2u64.pow(c.unwrap_or(8)) - 1;
        let max_dist = match max_seed_len {
            Some(max_seed_len) => {
                if max_seed_len < k || max_seed_len - k > 255 {
                    return Err(InvalidSeedingParameter::InvalidParameter(
                        "max seed length must be between k and k + 255",
                    ));
                }
                (max_seed_len - k) as u8
            }
            None => p.max_dist,
        };

        SeedingParameters::try_new(Profile::Noisy, k, s, w_min, w_max, q, max_dist, aux_len)
    }

    /// Create an IndexParameters instance based on a given read length.
    /// k, s, l, u, c and max_seed_len can be used to override determined parameters
    pub fn from_read_length(
        read_length: usize,
        k: Option<usize>,
        s: Option<usize>,
        w_min: Option<usize>,
        w_max: Option<usize>,
        c: Option<u32>,
        max_seed_len: Option<usize>,
        aux_len: u8,
    ) -> Result<Self, InvalidSeedingParameter> {
        let read_length_profile = READ_LENGTH_SETTINGS
            .iter()
            .find(|&rp| read_length <= rp.r_threshold)
            .expect("missing profile");

        let k = k.unwrap_or(read_length_profile.k);
        let s = s.unwrap_or((k as isize + read_length_profile.s_offset) as usize);
        let w_min = w_min.unwrap_or(read_length_profile.w_min);
        let w_max = w_max.unwrap_or(read_length_profile.w_max);
        let q = 2u64.pow(c.unwrap_or(8)) - 1;
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
            None => usize::clamp(
                read_length_profile.canonical_read_length.saturating_sub(70),
                k,
                255,
            ) as u8,
        };

        SeedingParameters::try_new(
            read_length_profile.profile.clone(),
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

    pub fn default_from_profile(profile: &Profile) -> SeedingParameters {
        match profile {
            Profile::Noisy => {
                Self::noisy(None, None, None, None, None, None, DEFAULT_AUX_LEN).unwrap()
            }
            Profile::Canonical50 => Self::default_from_read_length(50),
            Profile::Canonical75 => Self::default_from_read_length(75),
            Profile::Canonical100 => Self::default_from_read_length(100),
            Profile::Canonical125 => Self::default_from_read_length(125),
            Profile::Canonical150 => Self::default_from_read_length(150),
            Profile::Canonical250 => Self::default_from_read_length(250),
            Profile::Canonical400 => Self::default_from_read_length(400),
        }
    }

    pub fn is_custom(&self) -> bool {
        Self::default_from_profile(&self.profile) != *self
    }

    /// Returns a filename extension such as ".r100.sti",
    /// where `r100` is the profile name.
    ///
    /// If the profile was modified from the default, returns just ".sti".
    pub fn filename_extension(&self) -> String {
        if *self != SeedingParameters::from(&self.profile) {
            ".sti".to_string()
        } else {
            format!(".{}.sti", self.profile.identifier())
        }
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
        let main_hash_mask = 0xfffffffffc000000;
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
            Profile::Canonical250,
            k,
            s,
            w_min,
            w_max,
            q,
            max_dist,
            aux_len,
        )
        .unwrap();
        assert_eq!(seeding_parameters.profile, Profile::Canonical250);
        assert_eq!(seeding_parameters.randstrobe, rp);
        assert_eq!(seeding_parameters.syncmer, sp);

        let ip = SeedingParameters::default_from_read_length(canonical_read_length + 1);
        assert_eq!(ip.profile, Profile::Canonical250);
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
