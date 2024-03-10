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

function rstrobes() {
    target/debug/rstrobes "${@}"
}

cargo build

# Unit tests
cargo test

# should fail when unknown command-line option used
if rstrobes -G > /dev/null 2> /dev/null; then false; fi

# should succeed when only printing help
rstrobes -h > /dev/null

# Ensure the binary is available
samtools --version > /dev/null

# Single-end SAM
# TODO --chunk-size 3
# TODO -v
rstrobes --no-PG --eqx --rg-id 1 --rg SM:sample --rg LB:library tests/phix.fasta tests/phix.1.fastq > phix.se.sam
diff tests/phix.se.sam phix.se.sam
rm phix.se.sam

# Single-end SAM, M CIGAR operators
rstrobes --no-PG tests/phix.fasta tests/phix.1.fastq > phix.se.m.sam
if samtools view phix.se.m.sam | cut -f6 | grep -q '[X=]'; then false; fi

rm phix.se.m.sam

# Paired-end SAM
# TODO --chunk-size 3
rstrobes --no-PG --eqx --rg-id 1 --rg SM:sample --rg LB:library tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq > phix.pe.sam
diff tests/phix.pe.sam phix.pe.sam
rm phix.pe.sam

# Single-end PAF
rstrobes -x tests/phix.fasta tests/phix.1.fastq | tail -n 11 > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Single-end PAF (stdin input)
cat tests/phix.1.fastq | rstrobes -x tests/phix.fasta - | tail -n 11 > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Paired-end PAF
rstrobes -x tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq | tail -n 11 > phix.pe.paf
diff tests/phix.pe.paf phix.pe.paf
rm phix.pe.paf

# Build a separate index
rstrobes --no-PG -r 150 tests/phix.fasta tests/phix.1.fastq > without-sti.sam
rstrobes -r 150 -i tests/phix.fasta
rstrobes --no-PG -r 150 --use-index tests/phix.fasta tests/phix.1.fastq > with-sti.sam
diff without-sti.sam with-sti.sam
rm without-sti.sam with-sti.sam

# Create index requires -r or reads file
if rstrobes --create-index tests/phix.fasta > /dev/null 2> /dev/null; then false; fi

# --details output is proper SAM
rstrobes --details tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq 2> /dev/null | samtools view -o /dev/null
rstrobes --details tests/phix.fasta tests/phix.1.fastq 2> /dev/null | samtools view -o /dev/null

# Secondary alignments

# No secondary alignments on phix
rstrobes --no-PG tests/phix.fasta tests/phix.1.fastq > no-secondary.sam
rstrobes --no-PG -N 5 tests/phix.fasta tests/phix.1.fastq > with-secondary.sam
test $(samtools view -f 0x100 -c with-secondary.sam) -eq 0
rm no-secondary.sam with-secondary.sam

# Secondary alignments for repeated phiX
cp tests/phix.fasta repeated-phix.fasta
echo ">repeated_NC_001422" >> repeated-phix.fasta
sed 1d tests/phix.fasta >> repeated-phix.fasta
rstrobes --no-PG repeated-phix.fasta tests/phix.1.fastq > no-secondary.sam
rstrobes --no-PG -N 5 repeated-phix.fasta tests/phix.1.fastq > with-secondary.sam
test $(samtools view -f 0x100 -c with-secondary.sam) -gt 0

# Removing secondary alignments gives same result as not producing them in the first place
samtools view -h --no-PG -F 0x100 with-secondary.sam > with-secondary-only-primary.sam
diff no-secondary.sam with-secondary-only-primary.sam
rm no-secondary.sam with-secondary.sam with-secondary-only-primary.sam repeated-phix.fasta
