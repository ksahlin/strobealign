# Installation

![CI](https://github.com/ksahlin/strobealign/actions/workflows/ci.yml/badge.svg?branch=main)
[![Crates.io](https://img.shields.io/crates/v/strobealign.svg)](https://crates.io/crates/strobealign)
[![install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/strobealign/README.html)
[![Bioconda downloads](https://anaconda.org/bioconda/strobealign/badges/downloads.svg)](https://anaconda.org/bioconda/strobealign)

## Conda/Bioconda

Strobealign is [available from Bioconda](https://anaconda.org/bioconda/strobealign/).
Follow the [Bioconda setup instructions](https://bioconda.github.io/#usage). Then
```sh
conda create -n strobealign strobealign samtools  # install strobealign
conda activate strobealign                        # activate environment
strobealign --version                             # test the installation
```
You can also use other Conda-compatible package managers such as
[pixi](https://pixi.prefix.dev/) or [Mamba](https://mamba.readthedocs.io/).

## Homebrew

Strobealign is available in Homebrew core:
```sh
brew install strobealign
```

## Cargo

Strobealign [is available as a crate on crates.io](https://crates.io/crates/strobealign/).

To install it, you need a recent Rust version, which you can get via
[Rustup](https://rust-lang.org/tools/install/) (the Rust packaged with
Debian and Ubuntu is usually too old). This will give you the `cargo` command.
```sh
RUSTFLAGS='-C target-cpu=native' cargo install --locked strobealign
```

## From source

This will allow you to install the in-development version of strobealign.

You need to have Rust installed; see the previous section.
```sh
git clone https://github.com/ksahlin/strobealign
cd strobealign
RUSTFLAGS='-C target-cpu=native' cargo build --release
```
The binary is then available at `target/release/strobealign`.

## Python bindings

Experimental Python bindings were available in the C++ version of strobealign
(until version 0.17.0) and have not been ported to Rust.
We may add them back, but until then,
you will need to get strobealign 0.17.0,
and can then install them with `pip install .`.
The only documentation for the moment are the tests in `tests/*.py`.
