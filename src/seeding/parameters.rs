use std::fmt::Display;

use thiserror::Error;

use super::strobes::{DEFAULT_AUX_LEN, RandstrobeParameters};
use super::syncmers::SyncmerParameters;

/// A preset for seeding parameters
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Profile {
    Noisy,
    ReadLength50,
    ReadLength75,
    ReadLength100,
    ReadLength125,
    ReadLength150,
    ReadLength250,
    ReadLength400,
}

/*impl From<usize> for Profile {
    /// Given a read length, returns the corresponding closest canonical `Profile`
    ///
    /// This never returns the `Noisy` variant (to get it, construct it explicitly
    /// instead).
    fn from(read_length: usize) -> Self {
        // unwrap is fine because the last r_threshold is usize::MAX
        READ_LENGTH_SETTINGS
            .iter()
            .find(|&profile| read_length <= profile.r_threshold)
            .unwrap()
            .profile
            .clone()
    }
}*/

impl Profile {
    /// If this is a read-length based profile, returns the canonical length.
    /// Otherwise, returns None.
    fn length(&self) -> Option<usize> {
        match *self {
            Profile::Noisy => None,
            Profile::ReadLength50 => Some(50),
            Profile::ReadLength75 => Some(75),
            Profile::ReadLength100 => Some(100),
            Profile::ReadLength125 => Some(125),
            Profile::ReadLength150 => Some(150),
            Profile::ReadLength250 => Some(250),
            Profile::ReadLength400 => Some(400),
        }
    }

    /// Returns the identifier for the profile (used in the index file name).
    fn identifier(&self) -> String {
        if let Some(length) = self.length() {
            format!("r{}", length)
        } else {
            "noisy".to_string()
        }
    }
}

/// A u32 value is stored in `.sti` files to encode the profile. This trait
/// implementation is used to convert it back to a `Profile`.
impl TryFrom<u32> for Profile {
    type Error = ();

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        match value as u32 {
            x if x == Profile::Noisy as u32 => Ok(Profile::Noisy),
            x if x == Profile::ReadLength50 as u32 => Ok(Profile::ReadLength50),
            x if x == Profile::ReadLength75 as u32 => Ok(Profile::ReadLength75),
            x if x == Profile::ReadLength100 as u32 => Ok(Profile::ReadLength100),
            x if x == Profile::ReadLength125 as u32 => Ok(Profile::ReadLength125),
            x if x == Profile::ReadLength150 as u32 => Ok(Profile::ReadLength150),
            x if x == Profile::ReadLength250 as u32 => Ok(Profile::ReadLength250),
            x if x == Profile::ReadLength400 as u32 => Ok(Profile::ReadLength400),
            _ => Err(()),
        }
    }
}

impl Display for Profile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(length) = self.length() {
            write!(f, "canonical {} nt", length)
        } else {
            write!(f, "noisy")
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
        profile: Profile::ReadLength50,
        canonical_read_length: 50,
        r_threshold: 70,
        k: 18,
        s_offset: -4,
        w_min: 1,
        w_max: 4,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength75,
        canonical_read_length: 75,
        r_threshold: 90,
        k: 20,
        s_offset: -4,
        w_min: 1,
        w_max: 6,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength100,
        canonical_read_length: 100,
        r_threshold: 110,
        k: 20,
        s_offset: -4,
        w_min: 2,
        w_max: 6,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength125,
        canonical_read_length: 125,
        r_threshold: 135,
        k: 20,
        s_offset: -4,
        w_min: 3,
        w_max: 8,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength150,
        canonical_read_length: 150,
        r_threshold: 175,
        k: 20,
        s_offset: -4,
        w_min: 5,
        w_max: 11,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength250,
        canonical_read_length: 250,
        r_threshold: 375,
        k: 22,
        s_offset: -4,
        w_min: 6,
        w_max: 16,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength400,
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

#[derive(Error, Debug)]
pub enum InvalidSeedingParameter {
    #[error("Invalid seeding parameter: {0}")]
    InvalidParameter(&'static str),
}

impl SeedingParameters {
    pub fn new(profile: &Profile) -> Self {
        let c = 8;
        let q = 2u64.pow(c) - 1;
        let main_hash_mask = !0u64 << (9 + DEFAULT_AUX_LEN);

        match profile {
            Profile::Noisy => SeedingParameters {
                profile: profile.clone(),
                syncmer: SyncmerParameters::try_new(16, 12).unwrap(),
                randstrobe: RandstrobeParameters::try_new(2, 2, q, 84, main_hash_mask).unwrap(),
            },
            profile => {
                let read_length = profile.length().unwrap();
                let read_length_profile = READ_LENGTH_SETTINGS
                    .iter()
                    .find(|&rp| read_length <= rp.r_threshold)
                    .expect("missing profile");

                let s = (read_length_profile.k as isize + read_length_profile.s_offset) as usize;
                let max_dist =
                    usize::clamp(read_length.saturating_sub(70), read_length_profile.k, 255) as u8;

                SeedingParameters {
                    profile: profile.clone(),
                    syncmer: SyncmerParameters::try_new(read_length_profile.k, s).unwrap(),
                    randstrobe: RandstrobeParameters::try_new(
                        read_length_profile.w_min,
                        read_length_profile.w_max,
                        q,
                        max_dist,
                        main_hash_mask,
                    )
                    .unwrap(),
                }
            }
        }
    }

    pub fn with_aux_len(mut self, aux_len: u8) -> Result<Self, InvalidSeedingParameter> {
        if aux_len > 63 {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "aux length must be less than 64",
            ));
        }
        self.randstrobe.main_hash_mask = !0u64 << (9 + aux_len);

        Ok(self)
    }

    pub fn try_new(
        profile: Profile,
        k: usize,
        s: usize,
        w_min: usize,
        w_max: usize,
        q: u64,
        max_dist: u8,
    ) -> Result<Self, InvalidSeedingParameter> {
        let main_hash_mask = Self::new(&profile).randstrobe.main_hash_mask;
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

        SeedingParameters::try_new(Profile::Noisy, k, s, w_min, w_max, q, max_dist)
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
        )
    }

    pub fn default_from_read_length(read_length: usize) -> SeedingParameters {
        Self::from_read_length(read_length, None, None, None, None, None, None).unwrap()
    }

    pub fn default_from_profile(profile: &Profile) -> SeedingParameters {
        if let Some(length) = profile.length() {
            Self::default_from_read_length(length)
        } else {
            Self::noisy(None, None, None, None, None, None).unwrap()
        }
    }

    /// Returns whether the settings differ from the profile they are based on.
    pub fn is_custom(&self) -> bool {
        Self::default_from_profile(&self.profile) != *self
    }

    /// Returns a filename extension such as ".r100.sti",
    /// where `r100` is the profile name.
    ///
    /// If the profile was modified from the default, returns just ".sti".
    pub fn filename_extension(&self) -> String {
        if self.is_custom() {
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
        let seeding_parameters =
            SeedingParameters::try_new(Profile::ReadLength250, k, s, w_min, w_max, q, max_dist)
                .unwrap()
                .with_aux_len(aux_len)
                .unwrap();
        assert_eq!(seeding_parameters.profile, Profile::ReadLength250);
        assert_eq!(seeding_parameters.randstrobe, rp);
        assert_eq!(seeding_parameters.syncmer, sp);

        let ip = SeedingParameters::default_from_read_length(canonical_read_length + 1);
        assert_eq!(ip.profile, Profile::ReadLength250);
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
