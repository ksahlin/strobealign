![CI](https://github.com/ksahlin/strobealign/actions/workflows/ci.yml/badge.svg?branch=main)
[![Crates.io](https://img.shields.io/crates/v/strobealign.svg)](https://crates.io/crates/strobealign)
[![install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/strobealign/README.html)
[![Bioconda downloads](https://anaconda.org/bioconda/strobealign/badges/downloads.svg)](https://anaconda.org/bioconda/strobealign)

# Introduction

Strobealign is a read mapper that is typically significantly faster than other read mappers while achieving comparable or better accuracy, see the [performance evaluation](evaluation.md).

## Features

- Map single-end and paired-end reads
- Multithreading support
- Fast indexing (<1 minute for a human-sized reference genome using four cores)
- On-the-fly indexing by default. Optionally create an on-disk index.
- Output in standard SAM format or produce even faster results by writing PAF (without alignments)
- Strobealign is most suited for read lengths between 100 and 500 bp

## Background

Strobealign achieves its speedup by using a dynamic seed size obtained from [syncmer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7869670/)-thinned [strobemers](https://github.com/ksahlin/strobemers#what-is-a-strobemer).

For details, refer to [Strobealign: flexible seed size enables ultra-fast and accurate read alignment](https://doi.org/10.1186/s13059-022-02831-7). The paper describes v0.7.1 of the program.

### Videos

* For an introduction, see the 📺 [RECOMB-Seq video from 2022: “Flexible seed size enables ultra-fast and accurate read alignment”](https://www.youtube.com/watch?v=cn32telW63w) (12 minutes).
* For a more detailed presentation of the underlying seeding mechanism in strobealign (strobemers) see 📺 [“Efficient sequence similarity searches with strobemers”](https://www.youtube.com/watch?v=DS4tURz1Wio) (73 minutes).



## Evaluation

See [Performance evaluation](evaluation.md) for some measurements of mapping
accuracy and runtime using strobealign 0.7.

