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

## Debugging

When compiling strobealign, you can add `-DCMAKE_BUILD_TYPE=RelWithDebInfo` to
the `cmake` options to get debug symbols.

If needed, run `make` with `VERBOSE=1` to get more logging output.

## Testing

After CMake has been run, you can use this one-liner to compile strobealign and
run the tests:
```
make -j -C build && tests/run.sh
```

Whenever you make changes that could potentially affect mapping results, you can
run a more elaborate test that compares strobealign against a “baseline”
(know good) commit. Just run this script:
```
tests/compare-baseline.sh
```
The first time, it will download a *D. melanogaster* genome and some reads from
the Sequence Read Archive (SRA). Since the dataset is truncated to the first
100'000 reads, mapping it should take less than 30 seconds.

The baseline commit is configured in `tests/baseline-commit.txt`. The script
builds strobealign from that commit and runs it against the downloaded test
data, then builds strobealign as it is in your working copy and compares the
two produced BAM files. The baseline BAM is cached and re-used as long as the
baseline commit does not change.

## Making a release

* Update changelog (adjust section header)
* Bump version in `CMakeLists.txt`
* Bump version in `setup.py`
* Commit
* Push and wait for CI to pass
* Do `git tag ...` and `git push --tags`
* Make a release via the GitHub releases page
* Wait for the Bioconda bot to pick up the new release and then approve its
  version bump PR

## Style guide

Existing code in strobealign is currently not necessarily consistent with this
style guide, but new code should follow it.

* Use the `clang-format` code formatter with the `-style=file` option.
* Use Python-like naming of functions, methods, variables, etc. That is, it
  should be `ClassName`, `variable_name`, `method_name`, `CONSTANT`.
* The header guard of a file named `xyz.hpp` should be named
  `STROBEALIGN_XYZ_HPP`.

## Detailed SAM output (`--details`)

When `--details` is provided, the following additional SAM tags are output for
mapped reads.

`na`: Number of NAMs (seeds) found
`nr`: Whether NAM rescue was needed (1 if yes, 0 if not)
`mr`: Number of times mate rescue was attempted (local alignment of a mate to
  the expected region given its mate)
`al`: Number of times an attempt was made to extend a seed (by gapped or ungapped alignment)
`ga`: Number of times an attempt was made to extend a seed by gapped alignment
`X0`: Number of equally scored best alignments (greater than 1 for multimappers).
 For paired-end reads, the tag is output for both reads, but the value is
 identical and is the number of equally scored best alignment *pairs*.
