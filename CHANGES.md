# Strobealign Changelog

## development version

* #386: Parallelize indexing even more by using @alugowski’s
  [poolSTL](https://github.com/alugowski/) `pluggable_sort`.
  Indexing a human reference (measured on CHM13) now takes only ~45 s on a
  recent machine (using 8 threads).
* #376: Improve accuracy for read length 50 by optimizing the default
  indexing parameters. Paired-end accuracy increases by 0.3 percentage
  points on average. Single-end accuracy increases by 1 percentage point.
* #395: Previously, read length 75 used the same indexing parameters as length
  50, but the improved settings for length 50 are not the best for length 75.
  To avoid a decrease in accuracy, we introduced a new set of pre-defined
  indexing parameters for read length 75 (a new canonical read length).
* If `--details` is used, output `X0:i` SAM tag with the number of
  identically-scored best alignments
* #378: Added `-C` option for appending the FASTA or FASTQ comment to SAM
  output. (Idea and name of the option taken from BWA-MEM.)
* #371: Added `--no-PG` option for not outputting the PG SAM header
* Include [ZStr](https://github.com/mateidavid/zstr/) in our own repository
  instead of downloading it at build time. This should make it possible to
  build strobealign without internet access.

## v0.12.0 (2023-11-23)

* #293: Fix: When mapping single-end reads, many multimappers were previously
  assigned a high mapping quality. They now get assigned mapping quality zero
  as intended.
* #321: Fix: For paired-end reads that cannot be placed as proper pairs, we now
  prefer placing them onto the same chrosome instead of on different ones if
  there is a choice.
* #328: Adjust MAPQ computation for single-end reads.
* #318: Added a `--details` option mainly intended for debugging. When used,
  some strobealign-specific tags are added to the SAM output that inform about
  things like no. of seeds found, whether mate rescue was performed etc.
* #333: Fix matches ending too early in PAF output.
* #359, #367: Assign (single-end and paired-end) multimappers randomly to one
  of the candidate mapping locations to reduce biases.
* #347: Reduce memory usage by avoiding an unnecessary copy of reference
  contigs.

## v0.11.0 (2023-06-22)

* #278: Memory usage was reduced drastically due to a redesigned strobemer
  index memory layout. For the human genome, for example, RAM usage was reduced
  from 23 to 13 GiB. (Other changes increased RAM usage again slightly, see
  below.)

  Idea and implementation for this substantial improvement were contributed by
  Shaojun Pan (@psj1997) (supervised by Luis Pedro @luispedro) and originate in
  his work on a "strobealign-lm" (low memory) branch of strobealign. Thanks!
* #277, #285, PR #306: Support for very large references (exceeding ~20 Gbp) was
  added by switching from 32 bit to 64 bit strobemer indices.
  This was also enabled and made simpler by the memory layout changes.
  This increases RAM usage by 1 GiB for human-sized genomes.
* #313: Increased accuracy (especially on short single-end reads) due to
  "more random" syncmers. This increases memory usage again slightly so that we
  are at 14.7 GiB RAM usage for the human genome for this version of
  strobealign.
* #307: Indexing was further parallelized, cutting the time for index generation
  in about half for many cases.
* #289: Fixed missing CIGAR for secondary alignments.
* #212: SEQ and QUAL are set to `*` for secondary alignments as recommended
  by the SAM specification.
* #294: Updated the alignment library (SSW), which fixes some incorrect
  alignments.

## v0.10.0 (2023-06-05)

* #258: Fixed compilation on MinGW. Thanks @teepean.
* #260: Include full command line in the SAM PG header. Thanks @telmin.
* #20: By default, emit `M` CIGAR operations instead of `=` and `X`.
  Added option `--eqx` to use `=` and `X` as before.
* #265: Fixed overflowing read count statistics when processing $2^{31}$ reads
  or more. Thanks @telmin.
* #273: Fix handling of interleaved files using `/1` or `/2` suffixes

## v0.9.0 (2023-03-16)

* Added progress report (only shown if output is not a terminal; can be
  disabled with `--no-progress`)
* PR #250: Avoid overeager soft clipping by adding an “end bonus” to the
  alignment score if the alignment reaches the 5' or 3' end of the read.
  This is equivalent to penalizing soft-clipping and improves mapping
  accuracy, in particular for short reads, as candidate mapping sites with and
  without soft clipping are compared more fairly. Use `-L` to change the end
  bonus. (This emulates a feature found in BWA-MEM.)
* Issue #238: Fixed occasionally incorrect soft clipping.
* PR #239: Fixed an uninitialized variable that could lead to nondeterministic
  results.
* Issue #137: Compute TLEN (in SAM output) correctly
* PR #255: Added support for reading gzip-compressed reference FASTA files.
* Issue #222: Made it possible again to build strobealign from the release
  tarball (not only from the Git repository).

## v0.8.0 (2023-02-01)

### Summary

This is a large release with over 600+ commits since the previous one (0.7.1).
Much of the work that went into it was enabled by a
[Bioinformatics Long-Term Support](https://nbis.se/support/) grant through [National Bioinformatics Structure Sweden (NBIS)](https://nbis.se/), which is the [SciLifeLab Bioinformatics platform](https://www.scilifelab.se/units/nbis/).

A majority of the commits was focused on reorganizing the code to make it
easier to maintain, to read, to test and to change. This has already paid off
in the form of external contributions that would have been more difficult
without those changes.

Another focus was on usability and standards compliance: SAM output follows
the SAM specification more closely, error messages are better, there is
now a `--help` command-line option, some irrelevant logging output was
hidden, and the documentation was updated.

Mapping speed and accuracy remain mostly unaffected in this release, except for
one bugfix that increases mapping rate and accuracy at short read lengths
(<100 bp). Also, some unintended coverage spikes no longer occur.

Memory usage was decreased due to switching to a modified in-memory
representation of the index. For the human genome, for example, RAM usage went
from 28 GiB to 21 GiB. The reduction is smaller for more repetitive genomes.

Strobealign also gained the ability to pre-generate an index and save it to disk.
As indexing is quite fast, this is not as relevant as it is for other read
mappers, but important when processing many small libraries.

### Features

* PR #191 and PR #192: RAM usage was reduced significantly. For the
  human genome, for example, it went from from 28 GiB to 21 GiB. (Mapping
  runtime is unaffected.)
* Issue #121: Deterministic output: Mapped reads are now output in the order
  that they had in the input (this was previously not the case when using
  more than one thread).
* PR #48, PR #195: The index can now be pre-generated and saved to disk
  (`--create-index` and `--use-index`).
* Added `-h`/`--help` options.
* Added `--version` option.
* Issue #31: Added `--rg-id` option for adding read group IDs to each SAM
  record, and added `--rg` option for adding read group metadata
  (sample, library etc.) to the SAM header (example: `--rg 1 --rg SM:mysample`
  sets read group id to "1" and the sample name to "mysample").
* Suppress some logging output by default, and added option `-v` for showing
  full logging output.
* Issue #206: Added option `-U` for suppressing output of unaligned reads. Thanks @sjaenick.
* PR #213: Add support for interleaved reads. Thanks @luispedro.
* PR #213: Add support for reading from a pipe (including stdin). Thanks @luispedro.

### Bug fixes

* Issue #121: When more than one thread is used, changed behavior so that
  alignments are output in the order that the reads had in the input file.
* PR #114: Fix read-length estimate on input files with few reads.
* PR #78: Fixed incorrect template length (TLEN) sign.
* Issue #56: Unmapped reads did not have quality values in the SAM output.
* PR #65: Fixed a logic error. This changes some mapping qualites.
* Issue #141: Remove `/1` and `/2` suffixes from read names if they exist.
* Issue #142: Template length (TLEN) should be zero when the two reads in a pair
  map to different contigs.
* PR #149: Set mapping quality of unmapped reads to 0 for compatibility
  with Picard.
* Issue #148: Fixed SAM output when one read in a pair is mapped and the other
  not.
* Issue #141: Remove any `/1` and `/2` suffixes from read names.
* Invalid or missing input files no longer lead to a crash.
* Issue #199: No longer (inadvertently) ignore the `-l` and `-u` command-line
  options.

### Developer-oriented changes

* Switched to using a third-party library for parsing command-line arguments
* Reorganized the source code (a lot).
* In particular, moved unused code into `unused.cpp`.
* Started using a unit testing library
* Added Continuous Integration tests (GitHub Actions) for Linux and macOS.
* Switched to CMake as main build system.
* Maintain the version number in only one place (`CMakeLists.txt`).
* Add a simple `dumpstrobes` command-line utility that produces a BED file
  with the locations of all seeds (randstrobes).

### Other changes

* There is now a Conda package for strobealign on Bioconda.
* Issue #34: Disabled AVX2 CPU instructions by default as they are sometimes not
  available. Re-enable by running `cmake` with `-DENABLE_AVX=ON`.
* Added a changelog.

## v0.7.1 (2022-04-17)

* Changed to MIT license.
