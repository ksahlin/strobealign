# Usage

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


## Mapping-only mode

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

## Abundance estimation mode

Abundance estimation is used for metagenomic binning.

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


## Read profiles and canonical read lengths

Strobealign needs to build an index of the reference before it can map reads to it.

The optimal indexing parameters depend on the type of reads that are to be mapped.

Strobealign by default uses pre-defined sets of parameters that are optimized
for different read lengths. These *canonical read lengths* are
50, 75, 100, 125, 150, 250 and 400. When deciding which of the pre-defined
indexing parameter sets to use, strobealign chooses one whose canonical
read length is close to the average read length of the input.

The average read length of the input is normally estimated from the first
500 reads, but can also be explicitly set with `-r`.

In addition, it is possible to choose a profile optimized for “noisy”, that is,
error-prone reads. Use `-P noisy` on the command line to select it.
This is equivalent to `-k 16 -s 12 -l 2 -u 2 -m 100`.
We found these settings to increase accuracy on error-prone reads,
but note this comes at the cost of runtime.

## Pre-computing index files

By default, strobealign creates a new index every time the program is run.
On current CPUs (and using multiple cores), indexing a human-sized genome takes
less than 1 minute, which is not long compared to mapping many millions of reads.
However, for repeatedly mapping small libraries, it may be faster to pre-generate
an index on disk and use that.

To create an index, use the `--create-index` option.
Since strobealign needs to know the read profile, either provide some reads on
the command line as if you wanted to map them:

    strobealign --create-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz

This will use a read-length based profile based on the estimated read length.
You can also set the read length explicitly with `-r`:

    strobealign --create-index -t 8  -r 150 ref.fa

Or use the “noisy” profile with `-P noisy`:

    strobealign --create-index -t 8 -P noisy ref.fa

This creates a file named `ref.fa.X.sti` containing the strobemer index,
where `X` is an identifier for the read profile (such as `r50`, `r100`, `noisy`).
To use the index when mapping, provide option `--use-index` when doing the
actual mapping:

    strobealign --use-index -t 8 ref.fa reads.1.fastq.gz reads.2.fastq.gz | samtools ...

If you want to use the “noisy” profile, you also need to specify `-P noisy`
during mapping.

- Note that the `.sti` files are usually tied to a specific strobealign version.
  That is, when upgrading strobealign, the `.sti` files need to be regenerated.
  Strobealign detects whether the index was created with an incompatible
  version and refuses to load it.
- Index files are about four times as large as the reference.
