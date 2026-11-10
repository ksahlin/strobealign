# Usage


## Short-read mapping for variant calling

Map paired-end reads using 8 threads while creating a sorted BAM file:
```sh
strobealign -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools sort -o aln.bam
```
Piping directly into `samtools` avoids writing large intermediate SAM files.

Strobealign works best for read lengths between 100 and 500 bp. It estimates the
read length from the input and picks indexing parameters to match, see
[read profiles](#read-profiles).

In many cases,
you also need to add read group information so that downstream tools (GATK, bcftools, DeepVariant)
can identify the sample:
```sh
strobealign -t 8 --exclude-alt --rg-id=S1 --rg=SM:sample1 --rg=LB:lib1 --rg=PL:ILLUMINA \
    GRCh38.fasta reads.1.fastq.gz reads.2.fastq.gz | samtools sort -o aln.bam
```

```{admonition} GRCh38 with ALT contigs

When mapping to GRCh38 with ALT contigs such as the reference bundled with
BWA or GATK, you should use option `--exclude-alt`.

ALT contigs duplicate regions of the primary assembly, so without this option,
reads in these regions get MAPQ 0 and may be assigned to the ALT contig.

Strobealign recognizes ALT contigs by name
(`*_alt`, `*_decoy`, `HLA-*`) and prints a warning if it detects them while
`--exclude-alt` is not given.
```


## Mapping single-end reads

To map single-end reads, provide only the name of one FASTQ file:
```sh
strobealign -t 8 ref.fa reads.fq.gz | samtools sort -o aln.bam
```


## Mapping interleaved or mixed reads

To map reads that are interleaved (R1, R2, R1, R2, etc.),
use option `--interleaved`:
```sh
strobealign -t 8 --interleaved ref.fa reads.fq.gz | samtools sort -o aln.bam
```
Even mixed reads are accepted (where some reads are single-end and some are paired).


## Useful command-line options

* `-N INT` outputs up to INT secondary alignments
* `-U` suppresses unmapped reads
* `--eqx` writes `=`/`X` instead of `M` in CIGAR strings

See [Command-line options](#command-line-options) for a more complete list.


## Mapping reads without aligning them

The command-line option `-x` switches strobealign into mapping-only mode,
which is much faster because no base-level alignments are done.
The output will be [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md)
files instead of SAM. For example:
```sh
strobealign -x -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | igzip > output.paf.gz
```
Because compression with regular `gzip` can become a bottleneck when mapping with many
threads, we suggest using `igzip` to compresse the output, which is much faster.
It is part of [ISA-L](https://github.com/intel/isa-l/).
If it is not available, replace it with `pigz` or regular `gzip` in the
command.

PAF output includes only mapped reads. Unmapped reads are omitted. This is also
true for paired-end reads: If one of the reads is unmapped, only the mapped one
is output.


(abundances)=
## Estimating metagenomic abundances

With `--aemb`, strobealign outputs a table of contig abundances instead of
alignments, for use with metagenomic binners:
```sh
strobealign -t 8 --aemb contigs.fa reads.1.fastq.gz reads.2.fastq.gz > abundances.tsv
```
Each row holds a contig name and its abundance: the number of bases mapped to the
contig divided by its length, where a read mapping to *n* locations counts
1/*n* toward each. See [Abundance estimation mode](#aemb).

(index-files)=
## Pre-computing an index

An index of the reference FASTA file is built on the fly by default, that is,
every time the program is run.

This is often fine since indexing a human-sized genome takes less than 1 minute using multiple cores.
However, for repeatedly mapping small libraries, it may be faster to write the index to
disk first:
```sh
strobealign --create-index -t 8 -r 150 ref.fasta
strobealign --use-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools sort -o aln.bam
```

The read length needs to be provided because the index differs depending on the [read profile](#read-profiles).

Instead of using option `-r` to set the read length,
you can also provide the provide the read files on the command line as if you wanted to map them:

```sh
strobealign --create-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz
```
This will use a profile based on the estimated read length.

To create an index for the [noisy profile](#read-profiles), use `-P noisy`:
```sh
strobealign --create-index -t 8 -P noisy ref.fa
```

In all cases, this creates a file named `ref.fa.X.sti` containing the strobemer index,
where `X` is an identifier for the read profile (such as `r50`, `r100`, `noisy`).

To use the index when mapping, provide the FASTA file as usual and also
option `--use-index` when doing the actual mapping:
```sh
strobealign --use-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools ...
```
The FASTA file is needs to be available for mapping
(the reference is not stored in the `.sti` file).

If you want to use the “noisy” profile, you also need to specify `-P noisy`
during mapping.

- Note that the `.sti` files are usually tied to a specific strobealign version.
  That is, when upgrading strobealign, the `.sti` files need to be regenerated.
  Strobealign detects whether the index was created with an incompatible
  version and refuses to load it.
- Index files are about four times as large as the reference.
