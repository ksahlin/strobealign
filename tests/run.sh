#!/bin/bash
set -euo pipefail

if [[ $OSTYPE = linux-gnu ]]; then
    color="--color=always"
else
    color=""
fi

function strobealign() {
    echo "Testing '${@}'" >&2
    build/strobealign "${@}" 2> testlog.txt
}

function diff() {
    if ! env diff -u ${color} "$1" "$2"; then
      echo "Failure running 'diff $1 $2'"
      return 1
    fi
}


trap 'echo -e "\e[1;31mFailure\e[0m"' ERR

build/test-strobealign --no-intro

# should fail when unknown command-line option used
if strobealign -G > /dev/null 2> /dev/null; then false; fi

# should succeed when only printing help
strobealign -h > /dev/null

# Ensure the binary is available
samtools --version > /dev/null

# Single-end SAM
strobealign --no-PG --eqx --chunk-size 3 --rg-id 1 --rg SM:sample --rg LB:library -v tests/phix.fasta tests/phix.1.fastq > phix.se.sam
diff tests/phix.se.sam phix.se.sam
rm phix.se.sam

# Single-end SAM, M CIGAR operators
strobealign --no-PG tests/phix.fasta tests/phix.1.fastq > phix.se.m.sam
if samtools view phix.se.m.sam | cut -f6 | grep -q '[X=]'; then false; fi

rm phix.se.m.sam

# Paired-end SAM
strobealign --no-PG --eqx --chunk-size 3 --rg-id 1 --rg SM:sample --rg LB:library tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq > phix.pe.sam
diff tests/phix.pe.sam phix.pe.sam
rm phix.pe.sam

# Single-end PAF
strobealign -x tests/phix.fasta tests/phix.1.fastq | tail -n 11 > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Single-end PAF (stdin input)
cat tests/phix.1.fastq | strobealign -x tests/phix.fasta - | tail -n 11 > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Paired-end PAF
strobealign -x tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq | tail -n 11 > phix.pe.paf
diff tests/phix.pe.paf phix.pe.paf
rm phix.pe.paf

# Single-end abundance estimation 
strobealign --aemb tests/phix.fasta tests/phix.1.fastq > phix.abun.se.txt
diff tests/phix.abun.se.txt phix.abun.se.txt
rm phix.abun.se.txt

# Paired-end abundance estimation
strobealign --aemb tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq > phix.abun.pe.txt
diff tests/phix.abun.pe.txt phix.abun.pe.txt
rm phix.abun.pe.txt

# Build a separate index
strobealign --no-PG -r 150 tests/phix.fasta tests/phix.1.fastq > without-sti.sam
strobealign -r 150 -i tests/phix.fasta
strobealign --no-PG -r 150 --use-index tests/phix.fasta tests/phix.1.fastq > with-sti.sam
diff without-sti.sam with-sti.sam
rm without-sti.sam with-sti.sam

# Create index requires -r or reads file
if strobealign --create-index tests/phix.fasta > /dev/null 2> /dev/null; then false; fi

# --details output is proper SAM
strobealign --details tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq 2> /dev/null | samtools view -o /dev/null
strobealign --details tests/phix.fasta tests/phix.1.fastq 2> /dev/null | samtools view -o /dev/null

# Secondary alignments

# No secondary alignments on phix
strobealign --no-PG tests/phix.fasta tests/phix.1.fastq > no-secondary.sam
strobealign --no-PG -N 5 tests/phix.fasta tests/phix.1.fastq > with-secondary.sam
test $(samtools view -f 0x100 -c with-secondary.sam) -eq 0
rm no-secondary.sam with-secondary.sam

# Secondary alignments for repeated phiX
cp tests/phix.fasta repeated-phix.fasta
echo ">repeated_NC_001422" >> repeated-phix.fasta
sed 1d tests/phix.fasta >> repeated-phix.fasta
strobealign --no-PG repeated-phix.fasta tests/phix.1.fastq > no-secondary.sam
strobealign --no-PG -N 5 repeated-phix.fasta tests/phix.1.fastq > with-secondary.sam
test $(samtools view -f 0x100 -c with-secondary.sam) -gt 0

# Removing secondary alignments gives same result as not producing them in the first place
samtools view -h --no-PG -F 0x100 with-secondary.sam > with-secondary-only-primary.sam
diff no-secondary.sam with-secondary-only-primary.sam
rm no-secondary.sam with-secondary.sam with-secondary-only-primary.sam repeated-phix.fasta

strobealign -C --no-PG --rg-id=1 --rg=SM:sample --rg=LB:library tests/phix.fasta tests/phix.tags.fastq > with-tags.sam
diff tests/phix.tags.sam with-tags.sam
rm with-tags.sam

# Multi-block gzip
( head -n 4 tests/phix.1.fastq | gzip ; tail -n +5 tests/phix.1.fastq | gzip ) > multiblock.fastq.gz
strobealign --no-PG --eqx --rg-id 1 --rg SM:sample --rg LB:library tests/phix.fasta multiblock.fastq.gz > multiblock.sam
diff tests/phix.se.sam multiblock.sam
rm multiblock.fastq.gz multiblock.sam

# Single-end PAF with multi-context seeds
strobealign --mcs=always -x tests/phix.fasta tests/phix.1.fastq | tail -n 11 > phix.mcs.se.paf
diff tests/phix.mcs.se.paf phix.mcs.se.paf
rm phix.mcs.se.paf

echo -e "\e[32mSuccess\e[0m"
