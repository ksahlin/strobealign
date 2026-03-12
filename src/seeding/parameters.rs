use std::fmt::Display;

use thiserror::Error;

use super::strobes::{DEFAULT_AUX_LEN, RandstrobeParameters};
use super::syncmers::SyncmerParameters;

/// A preset for seeding parameters
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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

impl From<usize> for Profile {
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
}

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

/// Seeding parameters for a read-length-based profile
struct ReadLengthSettings {
    profile: Profile,
    r_threshold: usize,
    k: usize,
    s: usize,
    w_min: usize,
    w_max: usize,
}

static READ_LENGTH_SETTINGS: [ReadLengthSettings; 7] = [
    ReadLengthSettings {
        profile: Profile::ReadLength50,
        r_threshold: 70,
        k: 18,
        s: 14,
        w_min: 1,
        w_max: 4,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength75,
        r_threshold: 90,
        k: 20,
        s: 16,
        w_min: 1,
        w_max: 6,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength100,
        r_threshold: 110,
        k: 20,
        s: 16,
        w_min: 2,
        w_max: 6,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength125,
        r_threshold: 135,
        k: 20,
        s: 16,
        w_min: 3,
        w_max: 8,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength150,
        r_threshold: 175,
        k: 20,
        s: 16,
        w_min: 5,
        w_max: 11,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength250,
        r_threshold: 375,
        k: 22,
        s: 18,
        w_min: 6,
        w_max: 16,
    },
    ReadLengthSettings {
        profile: Profile::ReadLength400,
        r_threshold: usize::MAX,
        k: 23,
        s: 17,
        w_min: 5,
        w_max: 15,
    },
];

/// Settings that influence seed creation
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedingParameters {
    /// The base profile provided at construction time from which the settings are derived.
    /// If parameters were changed after construction, `is_custom()` returns false.
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
    /// Constructs a new SeedingParameters object with default parameters
    /// specific to the given profile. If a `usize` is provided, it is treated
    /// as a read length and the appropriate `Profile::ReadLength...` value is used.
    pub fn new<T: Into<Profile>>(t: T) -> Self {
        let profile = t.into();
        let bitcount = 8;
        let q = 2u64.pow(bitcount) - 1;
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
                    .find(|&rlp| read_length <= rlp.r_threshold)
                    .expect("missing profile"); // should not happen because the last r_threshold is usize::MAX

                let max_dist =
                    usize::clamp(read_length.saturating_sub(70), read_length_profile.k, 255) as u8;

                SeedingParameters {
                    profile: profile.clone(),
                    syncmer: SyncmerParameters::try_new(
                        read_length_profile.k,
                        read_length_profile.s,
                    )
                    .unwrap(),
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

    /// A shortcut for `SeedingParameters::new(Profile::Noisy)`
    pub fn noisy() -> Self {
        Self::new(Profile::Noisy)
    }

    /// Returns parameters with updated k, s and also max_dist. The latter is
    /// updated because it depends on k.
    pub fn with_k_s(
        mut self,
        k: Option<usize>,
        s: Option<usize>,
    ) -> Result<Self, InvalidSeedingParameter> {
        let new_k = k.unwrap_or(self.syncmer.k);
        let new_s = if let Some(s) = s {
            s
        } else {
            // When only k, but no s is not given, adjust s such that k - s stays the same.
            let diff = self.syncmer.k - self.syncmer.s;
            assert!(new_k >= diff);

            new_k - diff
        };
        self.syncmer = SyncmerParameters::try_new(new_k, new_s)?;

        self.randstrobe.max_dist = if let Some(length) = self.profile.length() {
            usize::clamp(length.saturating_sub(70), new_k, 255) as u8
        } else {
            // noisy
            100 - new_k as u8
        };

        Ok(self)
    }

    /// Returns parameters with updated w_min and w_max.
    pub fn with_window(
        mut self,
        w_min: Option<usize>,
        w_max: Option<usize>,
    ) -> Result<Self, InvalidSeedingParameter> {
        self.randstrobe = self.randstrobe.with_window(w_min, w_max)?;

        Ok(self)
    }

    /// Returns parameters with updated max_dist.
    ///
    /// Change k before this as k is used to compute max_dist from max_seed_length.
    pub fn with_max_seed_length(
        mut self,
        max_seed_length: usize,
    ) -> Result<Self, InvalidSeedingParameter> {
        if max_seed_length < self.syncmer.k || max_seed_length - self.syncmer.k > 255 {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "max seed length must be between k and k + 255",
            ));
        }

        self.randstrobe.max_dist = (max_seed_length - self.syncmer.k) as u8;

        Ok(self)
    }

    /// Returns parameters with updated main_hash_mask, computed from aux_len.
    pub fn with_aux_len(mut self, aux_len: u8) -> Result<Self, InvalidSeedingParameter> {
        if aux_len > 63 {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "aux length must be less than 64",
            ));
        }
        self.randstrobe.main_hash_mask = !0u64 << (10 + aux_len);

        Ok(self)
    }

    /// Returns parameters with updated q mask, computed from bitcount.
    pub fn with_bitcount(mut self, bitcount: u32) -> Result<Self, InvalidSeedingParameter> {
        if bitcount < 2 || bitcount > 63 {
            return Err(InvalidSeedingParameter::InvalidParameter(
                "bitcount must be between 2 and 63",
            ));
        }

        self.randstrobe.q = 2u64.pow(bitcount) - 1;

        Ok(self)
    }

    /// Returns whether the settings differ from the base profile (that was given at
    /// construction time).
    pub fn is_custom(&self) -> bool {
        Self::new(self.profile) != *self
    }

    /// Returns an index filename extension such as ".r100.sti",
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
        let seeding_parameters = SeedingParameters::new(250)
            .with_k_s(Some(k), Some(s))
            .unwrap()
            .with_window(Some(w_min), Some(w_max))
            .unwrap()
            .with_max_seed_length(max_dist as usize + k)
            .unwrap()
            .with_aux_len(aux_len)
            .unwrap()
            .with_bitcount(8)
            .unwrap();
        assert_eq!(seeding_parameters.profile, Profile::ReadLength250);
        assert_eq!(seeding_parameters.randstrobe, rp);
        assert_eq!(seeding_parameters.syncmer, sp);

        let ip = SeedingParameters::new(canonical_read_length + 1);
        assert_eq!(ip.profile, Profile::ReadLength250);
        assert_eq!(ip.randstrobe, rp);
        assert_eq!(ip.syncmer, sp);
    }

    #[test]
    fn test_seeding_parameters_similar_read_length() {
        let sp150 = SeedingParameters::new(150);
        let sp149 = SeedingParameters::new(149);
        let sp151 = SeedingParameters::new(151);

        assert_eq!(sp150, sp149);
        assert_eq!(sp150, sp151);
    }

    #[test]
    fn test_seeding_parameters_is_custom() {
        assert!(!SeedingParameters::new(100).is_custom());
        assert!(
            !SeedingParameters::new(100)
                .with_window(None, None)
                .unwrap()
                .is_custom()
        );

        assert!(
            SeedingParameters::new(100)
                .with_k_s(Some(17), None)
                .unwrap()
                .is_custom()
        );
        assert!(
            SeedingParameters::new(100)
                .with_k_s(None, Some(8))
                .unwrap()
                .is_custom()
        );
        assert!(
            SeedingParameters::new(100)
                .with_window(Some(3), Some(19))
                .unwrap()
                .is_custom()
        );
        assert!(
            SeedingParameters::new(100)
                .with_max_seed_length(123)
                .unwrap()
                .is_custom()
        );
    }
}
