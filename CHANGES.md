# Strobealign Changelog

## v0.9.0 (2023-03-16)

* Add progress report (only shown if output is not a terminal; can be
  disabled with `--no-progress`)
* PR #250: Avoid overeager soft clipping by adding an “end bonus” to the
  alignment score if the alignment reaches the 5' or 3' end of the read.
  This is equivalent to penalizing soft-clipping and improves mapping
  accuracy, in particular for short reads, as candidate mapping sites with and
  without soft clipping are compared more fairly. Use `-L` to change the end
  bonus. (This emulates a feature found in BWA-MEM.)
* Issue #238: Fix occasionally incorrect soft clipping.
* PR #239: Fix an uninitialized variable that could lead to nondeterministic
  results.
* Issue #137: Compute TLEN (in SAM output) correctly
* PR #255: Add support for reading gzip-compressed reference FASTA files.
* Issue #222: Make it possible again to build strobealign from the release
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
