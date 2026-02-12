# Contributing

Contributions to strobealign are very welcome!

- To report a bug,
  [open an issue](https://github.com/ksahlin/strobealign/issues/new).
- To suggest a small change, just submit a pull request.
- To suggest a larger change, it may be a good idea to first open an issue or
  to send an e-mail so we can discuss it.

Instead of one large PR, consider submitting multiple small, logically
self-contained PRs if it makes sense. This facilitates review and allows for a
more focused discussion.

## Building strobealign

Use `cargo build --release`. The compiled binary is at
`target/release/strobealign`.

Without `--release`, the compiler uses the default `debug` profile, which does
not enable as many optimizations. This makes compilation faster,
but results in an unnecessarily slow binary.

To create a slightly faster binary with optimizations specific to your CPU,
run `cargo` with the `RUSTFLAGS` environment variable set like this:
```
RUSTFLAGS='-C target-cpu=native' cargo build --release
```

You can also build and run the program in one step using `cargo run`, which may
look like this:
```
cargo run -r -- -t 8 tests/phix.fasta tests/phix.1.fastq
```
Option `-r` is short for `--release`. Options after `--` are passed to
strobealign.

Run the tests with `cargo test`.

There are additional tests that should be run with `tests/run.sh`.
(These have have not yet been ported to native Rust tests.)

## Building the legacy C++ version of strobealign

Until version 0.17.0, strobealign was written in C++. The code currently still
exists in the repository. To compile that code, run `cmake -B build`.
You can add specify a build type to set some useful compiler options by adding
one of these options to the command:
- Use `-DCMAKE_BUILD_TYPE=Release` when building a release (`-O3 -DNDEBUG`) (this is the default)
- Use `-DCMAKE_BUILD_TYPE=Debug` for development (`-O2 -g`)
- Use `-DCMAKE_BUILD_TYPE=RelWithDebInfo` for profiling (`-O3 -g -DNDEBUG`)

`-g` gives you debug symbols and `-DNDEBUG` disables assertions.

## Testing

Whenever you make changes that could potentially affect mapping results, you can
run a more elaborate test that compares strobealign against a “baseline”
(know good) commit. Just run this script:
```
tests/compare-baseline.sh
```
The first time, it will download a *D. melanogaster* genome and some reads from
the Sequence Read Archive (SRA). Since the dataset is truncated to the first
100'000 reads, mapping it should take less than 30 seconds.

The baseline commit is the most recent commit that contains the trailer
`Is-new-baseline: yes`.
The script builds strobealign from that commit
and runs it against the downloaded test data,
then builds strobealign as it is in your working copy and compares the
two produced BAM files. The baseline BAM is cached and re-used as long as the
baseline commit does not change.

## Making a release

* Update changelog (adjust section header)
* Bump version in `Cargo.toml`
* Bump version in `setup.py`
* Commit
* Push and wait for CI to pass
* Do `git tag ...` and `git push --tags`
* Make a release via the GitHub releases page
* Wait for the Bioconda bot to pick up the new release and then approve its
  version bump PR

## Detailed SAM output (`--details`)

When `--details` is provided, the following additional SAM tags are output for
mapped reads.

`na`: Number of NAMs (seeds) found
`nr`: The number of NAMs found during NAM rescue or -1 if no rescue was attempted
`mr`: Number of times mate rescue was attempted (local alignment of a mate to
  the expected region given its mate)
`al`: Number of times an attempt was made to extend a seed (by gapped or ungapped alignment)
`ga`: Number of times an attempt was made to extend a seed by gapped alignment
`X0`: Number of equally scored best alignments (greater than 1 for multimappers).
 For paired-end reads, the tag is output for both reads, but the value is
 identical and is the number of equally scored best alignment *pairs*.
