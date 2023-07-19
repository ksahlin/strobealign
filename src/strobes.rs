use std::collections::VecDeque;
use crate::syncmers::{Syncmer, SyncmerIterator};

pub struct RandstrobeParameters {
    pub w_min: usize,
    pub w_max: usize,
    pub q: u64,
    pub max_dist: u8,

// TODO
// void verify() const {
//     if (max_dist > 255) {
//         throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
//     }
// }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Randstrobe {
    pub hash: u64,
    pub strobe1_pos: usize,
    pub strobe2_pos: usize,
}

pub struct RandstrobeIterator<'a> {
    parameters: &'a RandstrobeParameters,
    syncmers: VecDeque<Syncmer>,
    syncmer_iterator: &'a mut SyncmerIterator<'a>,
}

impl<'a> RandstrobeIterator<'a> {
    pub fn new(syncmer_iterator: &'a mut SyncmerIterator<'a>, parameters: &'a RandstrobeParameters) -> RandstrobeIterator<'a> {
        RandstrobeIterator {
            parameters,
            syncmers: VecDeque::<Syncmer>::new(),
            syncmer_iterator,
        }
    }
}

impl<'a> Iterator for RandstrobeIterator<'a> {
    type Item = Randstrobe;
    fn next(&mut self) -> Option<Self::Item> {
        while self.syncmers.len() <= self.parameters.w_max {
            if let Some(syncmer) = self.syncmer_iterator.next() {
                self.syncmers.push_back(syncmer);
            } else {
                break;
            }
        }
        if self.syncmers.len() <= self.parameters.w_min {
            return None;
        }
        let strobe1 = self.syncmers[0];
        let max_position = strobe1.position + self.parameters.max_dist as usize;
        let mut min_val = u64::MAX;
        let mut strobe2 = self.syncmers[0]; // Defaults if no nearby syncmer

        for i in self.parameters.w_min .. self.syncmers.len() {
            debug_assert!(i <= self.parameters.w_max);
            if self.syncmers[i].position > max_position {
                break;
            }
            let b = (strobe1.hash ^ self.syncmers[i].hash) & self.parameters.q;
            let ones = b.count_ones() as u64;
            if ones < min_val {
                min_val = ones;
                strobe2 = self.syncmers[i];
            }
        }
        self.syncmers.pop_front();
        return Some(Randstrobe {
            hash: strobe1.hash.wrapping_add(strobe2.hash),
            strobe1_pos: strobe1.position,
            strobe2_pos: strobe2.position,
        })
    }
}
