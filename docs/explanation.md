# Explanation

## Multi-context seeds

Strobealign uses randstrobes as seeds, which in our case consist of two k-mers
("strobes") that are somewhat close to each other. When a seed is looked up
in the index, it is only found if both strobes match. By changing the way in
which the index is stored in v0.15.0, it became  possible to support
*multi-context seeds*. With those changes, strobealign falls back to looking
up only one of the strobes (a "partial seed") if the full seed cannot be found.
This results in better mapping rate and accuracy for read lengths of up to
about 200 nt.

Usage of multi-context seeds is enabled by default in strobealign since v0.16.0.
The strategy is to first search for all full seeds of the query and fall back to
partial seeds if *no* seeds could be found.

A slightly more accurate, but slower mode of using multi-context seeds is
available by using option `--mcs`: With it, the strategy is changed to a
fallback *per seed*: If an individual full seed cannot be found, its partial
version is looked up in the index.

## Collinear chaining

Strobealign uses collinear chaining as its default mapping and alignment method,
replacing the previous NAM approach. The collinear chaining algorithm reproduces
the method used in Minimap2.

Collinear chaining works by splitting strobemer hits into anchors (two anchors for full
hits, one for partial hits) and then constructing chains using these anchors with
a scoring function. Chains are created in O(N×h) time complexity, where N is the
number of anchors and h is a constant set at 50 by default.

Several command line options allow fine tuning of the chaining behavior: `-H`
controls the chaining look-back window heuristic, `--gd` sets the diagonal gap cost,
`--gl` sets the gap length cost, `--vp` determines the best-chain score threshold,
and `--sg` controls the maximum skip distance on the reference.
