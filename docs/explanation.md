# How it works

Strobealign achieves its speedup by using a dynamic seed size obtained from [syncmer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7869670/)-thinned [strobemers](https://github.com/ksahlin/strobemers#what-is-a-strobemer).

For details, refer to [Strobealign: flexible seed size enables ultra-fast and accurate read alignment](https://doi.org/10.1186/s13059-022-02831-7). The paper describes v0.7.1 of the program.

For an introduction, see also the 📺 [RECOMB-Seq video from 2022: “Flexible seed size enables ultra-fast and accurate read alignment”](https://www.youtube.com/watch?v=cn32telW63w) (12 minutes). For a more detailed presentation of the underlying seeding mechanism in strobealign (strobemers) see 📺 [“Efficient sequence similarity searches with strobemers”](https://www.youtube.com/watch?v=DS4tURz1Wio) (73 minutes).


(read-profiles)=
## Read profiles and canonical read lengths

Strobealign needs to build an index of the reference before it can map reads to it.

The optimal indexing parameters depend on the type of reads that are to be mapped.

Strobealign by default uses pre-defined sets of parameters that are optimized
for different read lengths. These *canonical read lengths* are
50, 75, 100, 125, 150, 250 and 400. When deciding which of the pre-defined
indexing parameter sets to use, strobealign chooses one whose canonical
read length is close to the average read length of the input.

The average read length of the input is normally estimated from the first
500 reads, but can also be explicitly set with `-r`.

In addition, it is possible to choose a profile optimized for “noisy”, that is,
error-prone reads. Use `-P noisy` on the command line to select it.
This is equivalent to `-k 16 -s 12 -l 2 -u 2 -m 100`.
We found these settings to increase accuracy on error-prone reads,
but note this comes at the cost of runtime.


## Multi-context seeds

Strobealign uses randstrobes as seeds, which in our case consist of two $k$-mers
("strobes") that are somewhat close to each other. When a seed is looked up
in the index, it is only found if both strobes match. By changing the way in
which the index is stored in v0.15.0, it became  possible to support
*multi-context seeds*. With those changes, strobealign falls back to looking
up only one of the strobes (a "partial seed") if the full seed cannot be found.
This results in better mapping rate and accuracy for read lengths of up to
about 200 nt.

Usage of multi-context seeds is enabled by default in strobealign since v0.16.0.

The default strategy (`--mcs=always`) is to attempt to look up each full seed and
then, if that lookup failed, to look up the partial seed.

A slightly faster, but less accurate mode of using multi-context seeds is
available by using option `--mcs=rescue`: With it, the strategy is changed to a
fallback *per read*: Only if *all* full seed lookups fail for a query will the
partial seeds be looked up.


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

Chains scoring at least `--vp` (0.7 by default) times the best chain are used
for the primary alignment. For single-end reads, lower-scoring chains are also
kept for the supplementary search if they score at least `--vs` (40 by
default). Chains are grouped by the read interval they cover, and such a chain
is only kept if it covers a new part of the read or scores at least `--vp`
times the best chain of its group.


### Supplementary alignments

Strobealign searches for a single supplementary alignment when the primary
alignment does not explain the whole read, which happens when a read crosses a
large variation. Once the primary is extended, every chain on the same contig
and the same strand as the primary is a candidate, unless it overlaps the
primary on the read by more than a fraction `--so` (0.5 by default) of what the
primary covers. Translocations and inversions are not detected.

Candidates are extended like the primary, but only over the read interval the
primary does not cover, so that the end bonus `-L` pulls the extension towards
the primary. Strobealign reports the best extension as a supplementary record
if it scores at least `--sm` (80 by default) and meets the primary within
`--ss` bases (50 by default), either on the read or on the reference. The
supplementary gets its own mapping quality, computed from the second best
candidate.

Duplication detection is run on the supplementary alignment the same way it is
run on the primary.

The supplementary alignment search runs for single-end reads only.

Use `--supplementary` to turn it on.


(aemb)=
## Abundance estimation mode

The command-line option `--aemb` switches strobealign into abundance estimation
mode, intended for metagenomic binning.
In this mode, strobealign outputs a single table with abundance values in
tab-separated value format instead of SAM or PAF.

Paired-end example:
```
strobealign -t 8 --aemb ref.fa reads.1.fastq.gz reads.2.fastq.gz > abundances.tsv
```
The output table contains one row for each contig of the reference.
The first column is the reference/contig id and the second its abundance.

The abundance is the number of bases mapped to a contig, divided by the length
of the contig. When a read maps to multiple locations, each of the *n* locations
is weighted 1/*n*. Note that this is in contrast to alignment mode, where by
default strobealign would output only a single alignment for each read
(this can be changed with `-N`).

Example output:
```tsv
chr1    0.030627
chr2    0.030627
chr3    0.030314
chr4    0.030489
chr5    0.030550
```

Further columns may be added to this table in future versions of strobealign.
