#!/bin/bash
#
# Run strobealign on some test data and compare the output against what the
# baseline version produces.
#
# - Test data is automatically downloaded if needed
#   (and put into tests/drosophila)
# - The baseline BAM file is generated if necessary

set -euo pipefail

# Fail early if pysam is missing
python3 -c 'import pysam'

ends="pe"
while getopts "s" opt; do
  case "${opt}" in
    s)
      ends=se  # single-end reads
      ;;
    \?)
      exit 1
      ;;
  esac
done

ref=tests/drosophila/ref.fasta
reads=(tests/drosophila/reads.1.fastq.gz)
if [[ ${ends} = "pe" ]]; then
  reads+=(tests/drosophila/reads.2.fastq.gz)
fi

# Ensure test data is available
tests/download.sh

source tests/baseline-commit.txt

baseline_bam=baseline/baseline-${baseline_commit}.${ends}.bam
baseline_binary=baseline/strobealign-${baseline_commit}
cmake_options=-DCMAKE_BUILD_TYPE=RelWithDebInfo

# Generate the baseline BAM if necessary
mkdir -p baseline
if ! test -f ${baseline_bam}; then
  if ! test -f ${baseline_binary}; then
    srcdir=$(mktemp -p . -d compile.XXXXXXX)
    git clone . ${srcdir}
    ( cd ${srcdir} && git checkout ${baseline_commit} )
    cmake ${srcdir} -B ${srcdir}/build ${cmake_options}
    if ! make -j 4 -C ${srcdir}/build strobealign; then
      exit 1
    fi
    mv ${srcdir}/build/strobealign ${baseline_binary}
    rm -rf "${srcdir}"
  fi
  ${baseline_binary} -t 4 ${ref} ${reads[@]} | samtools view -o ${baseline_bam}
fi

# Run strobealign. This recompiles from scratch to ensure we use consistent
# compiler options.
builddir=$(mktemp -p . -d build.XXXXXXX)
cmake . -B ${builddir} ${cmake_options}
make -j 4 -C ${builddir} strobealign
set -x
${builddir}/strobealign -t 4 ${ref} ${reads[@]} | samtools view -o head.bam
rm -rf ${builddir}

# Do the actual comparison
tests/samdiff.py ${baseline_bam} head.bam
