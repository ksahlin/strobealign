#!/bin/bash
set -xeuo pipefail

if [[ $OSTYPE = linux-gnu ]]; then
    color="--color=always"
else
    color=""
fi

function diff() {
    env diff -u ${color} "$1" "$2"
}

# Unit tests
build/test-strobealign

# should fail when unknown command-line option used
if strobealign -G > /dev/null 2> /dev/null; then false; fi

# should succeed when only printing help
strobealign -h > /dev/null

# Single-end SAM
strobealign --chunk-size 3 --rg-id 1 --rg SM:sample --rg LB:library -v tests/phix.fasta tests/phix.1.fastq | grep -v '^@PG' > phix.se.sam
diff tests/phix.se.sam phix.se.sam
rm phix.se.sam

# Paired-end SAM
strobealign --chunk-size 3 --rg-id 1 --rg SM:sample --rg LB:library tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq | grep -v '^@PG' > phix.pe.sam
diff tests/phix.pe.sam phix.pe.sam
rm phix.pe.sam

# Single-end PAF
strobealign -x tests/phix.fasta tests/phix.1.fastq | tail > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Single-end PAF (stdin input)
cat tests/phix.1.fastq | strobealign -x tests/phix.fasta - | tail > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Paired-end PAF
strobealign -x tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq | tail > phix.pe.paf
diff tests/phix.pe.paf phix.pe.paf
rm phix.pe.paf

# Build a separate index
strobealign -r 150 tests/phix.fasta tests/phix.1.fastq | grep -v '^@PG' > without-sti.sam
strobealign -r 150 -i tests/phix.fasta
strobealign -r 150 --use-index tests/phix.fasta tests/phix.1.fastq | grep -v '^@PG' > with-sti.sam
diff without-sti.sam with-sti.sam
rm without-sti.sam with-sti.sam

# Create index requires -r or reads file
if strobealign --create-index tests/phix.fasta > /dev/null 2> /dev/null; then false; fi
