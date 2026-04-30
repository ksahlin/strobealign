use fastrand::Rng;

/// Shuffles the best scoring items in a vector, assuming the data is sorted
/// This helps to ensure we pick a random best item in case there are multiple
/// equally good ones.
/// Returns the number of best items
pub fn shuffle_best<T, F, S>(items: &mut [T], score: F, rng: &mut Rng) -> usize
where
    F: Fn(&T) -> S,
    S: PartialEq + Copy,
{
    let Some(best) = items.first() else {
        return 0;
    };

    let best_score = score(best);

    let n_best = items.iter().take_while(|x| score(x) == best_score).count();

    if n_best > 1 {
        let idx = rng.usize(..n_best);
        items.swap(0, idx);
    }

    n_best
}
