#!/bin/bash
set -xeuo pipefail

# should fail when unknown command-line option used
if strobealign -G > /dev/null 2> /dev/null; then false; fi

# should succeed when only printing help
strobealign -h > /dev/null

d=tests

# Single-end SAM
strobealign -v $d/phix.fasta $d/phix.1.fastq > phix.se.sam
diff -u tests/phix.se.sam phix.se.sam
rm phix.se.sam

# Paired-end SAM
strobealign $d/phix.fasta $d/phix.1.fastq $d/phix.2.fastq > phix.pe.sam
diff -u tests/phix.pe.sam phix.pe.sam
rm phix.pe.sam

# Single-end PAF
strobealign -x $d/phix.fasta $d/phix.1.fastq | tail > phix.se.paf
diff -u tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Paired-end PAF
strobealign -x $d/phix.fasta $d/phix.1.fastq $d/phix.2.fastq | tail > phix.pe.paf
diff -u tests/phix.pe.paf phix.pe.paf
rm phix.pe.paf

# Unit tests
build/test-strobealign
