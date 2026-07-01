use std::fmt::{Display, Formatter};

#[derive(Default, Debug, PartialEq, Eq, Copy, Clone, clap::ValueEnum)]
pub enum McsStrategy {
    // For each strobemer, do a full lookup. If that did not generate a hit,
    // try a partial lookup.
    #[default]
    Always,

    // For each strobemer, do a full lookup. If after processing the entire
    // query, no hits were generated, do partial lookups of each strobemer.
    Rescue,

    // Do full lookups only.
    Off,

    // Do partial lookups only, that is, use only the first strobe.
    FirstStrobe,
}

impl Display for McsStrategy {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            McsStrategy::Always => write!(f, "Always"),
            McsStrategy::Rescue => write!(f, "Rescue"),
            McsStrategy::Off => write!(f, "Off"),
            McsStrategy::FirstStrobe => write!(f, "FirstStrobe"),
        }
    }
}
