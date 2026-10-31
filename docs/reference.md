# Reference

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
