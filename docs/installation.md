# Installation

### Conda

Strobealign is available from [Bioconda](https://bioconda.github.io/).

1. Follow the [Bioconda setup instructions](http://bioconda.github.io/#usage)
2. Install strobealign into a new Conda environment:
   ```
   conda create -n strobealign strobealign samtools
   ```
3. Activate the environment that was just created:
   ```
   conda activate strobealign
   ```
4. Run strobealign:
   ```
   strobealign --version
   ```

### From source

To compile strobealign from sources, you need a somewhat recent Rust
version, which you can obtain via
[Rustup](https://rust-lang.org/tools/install/).
At the time of writing,
the Rust versions included in current Debian and Ubuntu releases are too old.

When Rust is available, you can compile strobealign with
```
RUSTFLAGS='-C target-cpu=native' cargo build --release
```
The resulting binary will then be available at `target/release/strobealign`.


See the [contributing instructions](#contributing) for how to compile
strobealign as a developer.

### Python bindings

Experimental Python bindings were available in the C++ version of strobealign
(until version 0.17.0) and have not been ported to Rust.
We may add them back, but until then,
you will need to get strobealign 0.17.0,
and can then install them with `pip install .`.
The only documentation for the moment are the tests in `tests/*.py`.
