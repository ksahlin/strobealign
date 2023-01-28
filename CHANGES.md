# Strobealign Changelog

## development version

### Bug fixes

* Issue #121: When more than one thread is used, changed behavior so that
  alignments are output in the order that the reads had in the input file.
* PR #114: Fix read-length estimate on input files with few reads.
* PR #78: Fixed incorrect template length (TLEN) sign.
* Issue #56: Unmapped reads did not have quality values in the SAM output.
* PR #65: Fixed a logic error. This changes some mapping qualites.
* Issue #142: Template length (TLEN) should be zero when the two reads in a pair
  map to different contigs.
* Issue #149: Changed mapping quality for unmapped reads to 0 for compatibility
  with Picard.
* Issue #148: Fixed SAM output when one read in a pair is mapped and the other
  not.

### Features

* PR #191 and PR #192: RAM usage has been decreased significantly. For the
  human genome, for example, it went from from 28 GiB to 21 GiB. (Mapping
  runtime is unaffected.)
* There is now a Conda package for StrobeAlign on Bioconda.
* PR #48, PR #195: The index can now be pre-generated and saved to disk.
* Invalid or missing input files no longer lead to a crash.
* Added `-h`/`--help` options.
* Added `--version` option.
* Added `--rg-id` option for adding read group IDs to each SAM record.
* Added `--rg` option for adding read group metadata (sample, library etc.)
  to the SAM header (example: `--rg 1 --rg SM:mysample` sets read group id to
  "1" and the sample name to "mysample").
* Suppress some logging output by default.
* Added option `-v` for showing full logging output.
* Issue #206: Added option `-U` for suppressing output of unaligned reads. Thanks @sjaenick.

### Developer-oriented changes

* Use a third-party library for parsing command-line arguments.
* Reorganized the source code (a lot).
* In particular, moved unused code into `unused.cpp`.
* Added Continuous Integration tests (GitHub Actions) for Linux and macOS.
* Switched to CMake as main build system.
* Maintain the version number in only one place (`CMakeLists.txt`).

### Other changes

* Issue #34: Disabled AVX2 CPU instructions by default as they are sometimes not
  available. Re-enable by running `cmake` with `-DENABLE_AVX=ON`.


## v0.7.1 (2022-04-17)

* Changed to MIT license.
