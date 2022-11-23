# StrobeAlign Changelog

## development version

### Bug fixes

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

* There is now a Conda package for StrobeAlign on Bioconda.
* The index can now be pre-generated and safed to disk.
* Invalid or missing input files no longer lead to a crash.
* Added `-h`/`--help` options.
* Added `--version` option.
* Suppress some logging output by default.
* Added option `-v` for showing full logging output.

### Developer-oriented changes

* Use a third-party library for parsing command-line arguments.
* Reorganized the source code (a lot).
* In particular, moved unused code into `unused.cpp`.
* Added Continuous Integration tests (GitHub Actions) for Linux and macOS.
* Switched to CMake as main build system.
* Maintain the version number in only one place (`CMakeLists.txt`).

### Other changes

* Issue #34: Disabled AVX2 CPU instructions by default.


## v0.7.1 (2022-04-17)

* Changed to MIT license.
