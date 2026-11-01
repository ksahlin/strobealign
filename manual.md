# strobealign manual

This document contains the detailed documentation for strobealign.
See the [README](README.md) for installation and a quick start.

- [Background](#background)
- [Mapping-only mode](#mapping-only-mode)
- [Abundance estimation mode](#abundance-estimation-mode-for-metagenomic-binning)
- [Command-line options](#command-line-options)
- [Read profiles and canonical read lengths](#read-profiles-and-canonical-read-lengths)
- [Index files](#index-files)
- [Explanation](#explanation)
- [Python bindings](#python-bindings)
- [Evaluation](#evaluation)

## Background

Strobealign achieves its speedup by using a dynamic seed size obtained from [syncmer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7869670/)-thinned [strobemers](https://github.com/ksahlin/strobemers#what-is-a-strobemer).

For details, refer to [Strobealign: flexible seed size enables ultra-fast and accurate read alignment](https://doi.org/10.1186/s13059-022-02831-7). The paper describes v0.7.1 of the program.

For an introduction, see also the 📺 [RECOMB-Seq video from 2022: “Flexible seed size enables ultra-fast and accurate read alignment”](https://www.youtube.com/watch?v=cn32telW63w) (12 minutes). For a more detailed presentation of the underlying seeding mechanism in strobealign (strobemers) see 📺 [“Efficient sequence similarity searches with strobemers”](https://www.youtube.com/watch?v=DS4tURz1Wio) (73 minutes).

## Mapping-only mode

The command-line option `-x` switches strobealign into mapping-only mode,
in which it will output [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md)
files instead of SAM. For example:
```
strobealign -x -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | igzip > output.paf.gz
```
`igzip` is a faster version of gzip that is part of
[ISA-L](https://github.com/intel/isa-l/).
If it is not available, replace it with `pigz` or regular `gzip` in the
command.

PAF output includes only mapped reads. Unmapped reads are omitted. This is also
true for paired-end reads: If one of the reads is unmapped, only the mapped one
is output.

## Abundance estimation mode (for metagenomic binning)

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

Further columns may be added to this table in future versions of strobealign.


## Command-line options

Please run `strobealign --help` to see the most up-to-date list of command-line
options. Some important ones are:

* `-r`: Mean read length. If given, this overrides the read length estimated
  from the input file(s). This is usually only required in combination with
  `--create-index`, see [index files](#index-files).
* `-t N`, `--threads=N`: Use N threads (both for mapping and indexing).
* `--eqx`: Emit `=` and `X` CIGAR operations instead of `M`.
* `-x`: Only map reads, do not do no base-level alignment. This switches the
  output format from SAM to [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md).
* `--aemb`: Output estimated abundance value of each contig, see section above.
* `--rg-id=ID`: Add RG tag to each SAM record.
* `--rg=TAG:VALUE`: Add read group metadata to the SAM header. This can be
  specified multiple times. Example: `--rg-id=1 --rg=SM:mysamle --rg=LB:mylibrary`.
* `-N INT`: Output up to INT secondary alignments. By default, no secondary
  alignments are output.
* `-U`: Suppress output of unmapped single-end reads and pairs in which both
  reads are unmapped.
* `--use-index`: Use a pre-generated index instead of generating a new one.
* `--create-index`, `-i`: Generate a strobemer index file (`.sti`) and write it
  to disk next to the input reference FASTA. Do not map reads. If read files are
  provided, they are used to estimate read length. See [index files](#index-files).

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

## Index files

### Pre-computing an index (`.sti`)

By default, strobealign creates a new index every time the program is run.
On current CPUs (and using multiple cores), indexing a human-sized genome takes
less than 1 minute, which is not long compared to mapping many millions of reads.
However, for repeatedly mapping small libraries, it may be faster to pre-generate
an index on disk and use that.

To create an index, use the `--create-index` option.
Since strobealign needs to know the read profile, either provide some reads on
the command line as if you wanted to map them:

    strobealign --create-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz

This will use a read-length based profile based on the estimated read length.
You can also set the read length explicitly with `-r`:

    strobealign --create-index -t 8  -r 150 ref.fa

Or use the “noisy” profile with `-P noisy`:

    strobealign --create-index -t 8 -P noisy ref.fa

This creates a file named `ref.fa.X.sti` containing the strobemer index,
where `X` is an identifier for the read profile (such as `r50`, `r100`, `noisy`).
To use the index when mapping, provide option `--use-index` when doing the
actual mapping:

    strobealign --use-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools ...

If you want to use the “noisy” profile, you also need to specify `-P noisy`
during mapping.

- Note that the `.sti` files are usually tied to a specific strobealign version.
  That is, when upgrading strobealign, the `.sti` files need to be regenerated.
  Strobealign detects whether the index was created with an incompatible
  version and refuses to load it.
- Index files are about four times as large as the reference.


## Explanation

### Multi-context seeds

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

### Collinear chaining

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

## Python bindings

Experimental Python bindings were available in the C++ version of strobealign
(until version 0.17.0) and have not been ported to Rust.
We may add them back, but until then,
you will need to get strobealign 0.17.0,
and can then install them with `pip install .`.
The only documentation for the moment are the tests in `tests/*.py`.

## Evaluation

See [Performance evaluation](evaluation.md) for some measurements of mapping
accuracy and runtime using strobealign 0.7.
