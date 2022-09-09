#!/bin/bash
set -xeuo pipefail

# should fail when unknown command-line option used
if strobealign -G > /dev/null; then false; fi

# should succeed when only printing help
strobealign -h > /dev/null

d=tests

# Paired-end test
strobealign $d/phix.fasta $d/phix.1.fastq $d/phix.2.fastq > phix.pe.sam
diff -u tests/phix.pe.sam phix.pe.sam
rm phix.pe.sam

# Single-end test
strobealign $d/phix.fasta $d/phix.1.fastq > phix.se.sam
diff -u tests/phix.se.sam phix.se.sam
rm phix.se.sam
