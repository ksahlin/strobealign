![CI](https://github.com/ksahlin/strobealign/actions/workflows/ci.yml/badge.svg?branch=main)
[![install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/strobealign/README.html)

# strobealign: A fast short-read aligner

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

For an introduction, see also the 📺 [RECOMB-Seq video from 2022: “Flexible seed size enables ultra-fast and accurate read alignment”](https://www.youtube.com/watch?v=cn32telW63w) (12 minutes). For a more detailed presentation of the underlying seeding mechanism in strobealign (strobemers) see 📺 [“Efficient sequence similarity searches with strobemers”](https://www.youtube.com/watch?v=DS4tURz1Wio) (73 minutes).

## Table of contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Command-line options](#command-line-options)
4. [Index files](#index-files)
5. [Changelog](#changelog)
6. [Contributing](#contributing)
7. [Evaluation](#evaluation)
8. [Citation](#citation)
9. [Version info](#version-info)
10. [License](#license)

## Installation

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

To compile from the source, you need CMake,
a recent `g++` (tested with version 8), [zlib](https://zlib.net/), `pkg-config` and
[ISA-L](https://github.com/intel/isa-l/) (Intel Intelligent Storage Acceleration Library).
On Debian/Ubuntu,
running `sudo apt-get install build-essential libisal-dev cmake` should take care of this.
On macOS, use `brew install pkg-config isa-l`.

Then do the following:
```
git clone https://github.com/ksahlin/strobealign
cd strobealign
cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native"
cmake --build build -j $(nproc)
```
The resulting binary is `build/strobealign`.

The binary is tailored to the CPU the compiler runs on.
If it needs to run on other machines, use this `cmake` command instead for compatibility with most x86-64 CPUs in use today:
```
cmake -B build -DCMAKE_C_FLAGS="-msse4.2" -DCMAKE_CXX_FLAGS="-msse4.2"
```

See the [contributing instructions](#contributing) for how to compile strobealign
as a developer.

### Python bindings

Experimental Python bindings can be installed with
`pip install .`. The only documentation for the moment are the tests in
`tests/*.py`.

## Usage

To align paired-reads against a reference FASTA and produce a sorted BAM file,
using eight threads:
```
strobealign -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools sort -o sorted.bam
```

For single-end reads:
```
strobealign -t 8 ref.fa reads.fastq.gz | samtools sort -o sorted.bam
```

For mixed reads (the input file can contain both single and paired-end reads):
```
strobealign -t 8 ref.fa --interleaved reads.fastq.gz | samtools sort -o sorted.bam
```

In alignment mode, strobealign produces SAM output. By piping the output
directly into `samtools`, the above commands avoid creating potentially large
intermediate SAM files and also reduce disk I/O.

To produce unsorted BAM, use `samtools view` instead of `samtools sort`.


### Mapping-only mode

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

## Index files

### Background

Strobealign needs to build an index (strobemer index) of the reference before
it can map reads to it.
The optimal indexing parameters depend on the length of the input reads.
There are pre-defined sets of parameters that are optimized for different read
lengths. These *canonical read lengths* are
50, 75, 100, 125, 150, 250 and 400. When deciding which of the pre-defined
indexing parameter sets to use, strobealign chooses one whose canonical
read length is close to the average read length of the input.

The average read length of the input is normally estimated from the first
500 reads, but can also be explicitly set with `-r`.

### Pre-computing an index (`.sti`)

By default, strobealign creates a new index every time the program is run.
Depending on CPU, indexing a human-sized genome takes
1 to 2 minutes, which is not long compared to mapping many millions of reads.
However, for repeatedly mapping small libraries, it is faster to pre-generate
an index on disk and use that.

To create an index, use the `--create-index` option.
Since strobealign needs to know the read length, either provide it with
read file(s) as if you wanted to map them:

    strobealign --create-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz

Or set the read length explicitly with `-r`:

    strobealign --create-index -t 8 ref.fa -r 150

This creates a file named `ref.fa.rX.sti` containing the strobemer index,
where `X` is the canonical read length that the index is optimized for (see
above).
To use the index when mapping, provide option `--use-index` when doing the
actual mapping:

    strobealign --use-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools ...

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


## Changelog

See [Changelog](CHANGES.md).


## Contributing

See [Contributing](CONTRIBUTING.md).


## Evaluation

See [Performance evaluation](evaluation.md) for some measurements of mapping
accuracy and runtime using strobealign 0.7.


## Citation

Sahlin, K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. Genome Biol 23, 260 (2022). [https://doi.org/10.1186/s13059-022-02831-7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02831-7)


## License

Strobealign is available under the MIT license,
see [LICENSE](https://github.com/ksahlin/strobealign/blob/main/LICENSE).
