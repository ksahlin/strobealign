#!/bin/bash
set -euo pipefail

# should fail when unknown command-line option used
#if strobealign -G; then false; fi

# should succeed when only printing help
strobealign -h

d=tests
strobealign $d/phix.fasta $d/phix.1.fastq $d/phix.2.fastq | \
  samtools sort --no-PG -o phix.bam -
test $(samtools view -c -F 4 phix.bam) = 200
