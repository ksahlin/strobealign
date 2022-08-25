#!/bin/bash
set -xeuo pipefail

# should fail when unknown command-line option used
#if strobealign -G; then false; fi

# should succeed when only printing help
strobealign -h

d=tests
strobealign $d/phix.fasta $d/phix.1.fastq $d/phix.2.fastq > phix.sam
diff -u tests/phix.sam phix.sam
rm phix.sam
