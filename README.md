![CI](https://github.com/ksahlin/strobealign/actions/workflows/ci.yml/badge.svg?branch=main)
[![Crates.io](https://img.shields.io/crates/v/strobealign.svg)](https://crates.io/crates/strobealign)
[![install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/strobealign/README.html)
[![Bioconda downloads](https://anaconda.org/bioconda/strobealign/badges/downloads.svg)](https://anaconda.org/bioconda/strobealign)

# strobealign

Strobealign is a fast and accurate read mapper. It aligns single-end and
paired-end reads to a reference genome and outputs SAM, PAF, or contig
abundances. Indexing is fast enough (under a minute for a human genome) that it is
done on the fly by default.

## Getting started

```sh
# Install (see below for building from source)
conda create -n strobealign strobealign samtools
conda activate strobealign

# Paired-end reads -> sorted BAM
strobealign -t 8 ref.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o aln.bam

# Single-end reads
strobealign -t 8 ref.fa reads.fq.gz | samtools sort -o aln.bam

# Mapping only, no base-level alignment (PAF output, faster)
strobealign -x -t 8 ref.fa reads_1.fq.gz reads_2.fq.gz > aln.paf

# Short reads for variant calling: add read group, exclude GRCh38 ALT contigs
strobealign -t 8 --rg-id=S1 --rg=SM:sample1 --exclude-alt GRCh38.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o aln.bam

# Long reads
# TODO

# Metagenomic abundance estimation (for binning)
strobealign -t 8 --aemb contigs.fa reads_1.fq.gz reads_2.fq.gz > abundances.tsv
```

## Installation

Strobealign is available from [Bioconda](https://bioconda.github.io/).
Follow the [Bioconda setup instructions](https://bioconda.github.io/#usage). Then
```sh
conda create -n strobealign strobealign samtools  # install strobealign
conda activate strobealign                        # activate environment
strobealign --version                             # test the installation
```

To compile from source, you need a recent Rust version, which you can get via
[Rustup](https://rust-lang.org/tools/install/) (the Rust packaged with
Debian and Ubuntu is usually too old). Then
```sh
git clone https://github.com/ksahlin/strobealign
cd strobealign
RUSTFLAGS='-C target-cpu=native' cargo build --release
```
The binary is then available at `target/release/strobealign`.

## Usage

### Short-read mapping

```sh
strobealign -t 8 ref.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o aln.bam  # paired-end
strobealign -t 8 ref.fa reads.fq.gz | samtools sort -o aln.bam                  # single-end
strobealign -t 8 --interleaved ref.fa reads.fq.gz | samtools sort -o aln.bam    # interleaved or mixed
```
Piping directly into `samtools` avoids writing large intermediate SAM files.
Strobealign works best for read lengths between 100 and 500 bp. It estimates the
read length from the input and picks indexing parameters to match (see
[read profiles](manual.md#read-profiles-and-canonical-read-lengths)).

Useful options: `-N INT` outputs up to INT secondary alignments, `-U` suppresses
unmapped reads, `--eqx` writes `=`/`X` instead of `M` in CIGAR strings, and `-x`
skips base-level alignment and outputs [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md)
(see [mapping-only mode](manual.md#mapping-only-mode)).

### Short-read mapping for SNP and variant calling

Add read group information so that downstream tools (GATK, bcftools, DeepVariant)
can identify the sample:
```sh
strobealign -t 8 --rg-id=S1 --rg=SM:sample1 --rg=LB:lib1 --rg=PL:ILLUMINA \
    ref.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o aln.bam
```

When mapping to GRCh38 with ALT contigs (such as the reference bundled with
BWA or GATK), use `--exclude-alt`. ALT contigs duplicate regions of the primary
assembly, so without this option, reads in these regions get MAPQ 0 and may be
assigned to the ALT contig. Strobealign recognizes ALT contigs by name
(`*_alt`, `*_decoy`, `HLA-*`) and prints a warning if it detects them while
`--exclude-alt` is not given.

### Long-read mapping

TODO

### Metagenomic abundance estimation

With `--aemb`, strobealign outputs a table of contig abundances instead of
alignments, for use with metagenomic binners:
```sh
strobealign -t 8 --aemb contigs.fa reads_1.fq.gz reads_2.fq.gz > abundances.tsv
```
Each row holds a contig name and its abundance: the number of bases mapped to the
contig divided by its length, where a read mapping to *n* locations counts
1/*n* toward each.

### Pre-computing an index

The index is built on the fly by default. When mapping many small libraries
against the same reference, you can save time by writing it to disk once:
```sh
strobealign --create-index -t 8 -r 150 ref.fa       # creates ref.fa.r150.sti
strobealign --use-index -t 8 ref.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o aln.bam
```
See [index files](manual.md#index-files) for details.

## Further documentation

- [Manual](manual.md): all command-line options, read profiles, index files,
  and how strobealign works
- [Changelog](CHANGES.md)
- [Contributing](CONTRIBUTING.md)
- Run `strobealign --help` for the full list of options

## Citation

If using v0.17.0 or later:

*Tolstoganov, I., Martin, M., Buchin, N., and Sahlin, K. Multi-context seeds enable fast and high-accuracy read mapping. Genome Biol (2026). [https://doi.org/10.1186/s13059-026-04017-x](https://link.springer.com/article/10.1186/s13059-026-04017-x)*

If using v0.16.1 or earlier:

*Sahlin, K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. Genome Biol 23, 260 (2022). [https://doi.org/10.1186/s13059-022-02831-7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02831-7)*

## License

Strobealign is available under the MIT license, see [LICENSE](LICENSE).
