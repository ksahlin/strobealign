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
threads=4
mcs=0
while getopts "st:m" opt; do
  case "${opt}" in
    t)
      threads=$OPTARG
      ;;
    s)
      ends=se  # single-end reads
      ;;
    m)
      mcs=1
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

baseline_commit=$(< tests/baseline-commit.txt)

baseline_binary=baseline/strobealign-${baseline_commit}
cmake_options=-DCMAKE_BUILD_TYPE=RelWithDebInfo
extra_ext=""
strobealign_options="-t ${threads}"
if [[ ${mcs} == 1 ]]; then
  extra_ext=".mcs"
  strobealign_options="${strobealign_options} --mcs"
fi
baseline_bam=baseline/bam/${baseline_commit}.${ends}${extra_ext}.bam

# Generate the baseline BAM if necessary
mkdir -p baseline/bam
if ! test -f ${baseline_bam}; then
  if ! test -f ${baseline_binary}; then
    srcdir=$(mktemp -p . -d compile.XXXXXXX)
    git clone . ${srcdir}
    ( cd ${srcdir} && git checkout -d ${baseline_commit} )
    cmake ${srcdir} -B ${srcdir}/build ${cmake_options}
    if ! make -s -j 4 -C ${srcdir}/build strobealign; then
      exit 1
    fi
    mv ${srcdir}/build/strobealign ${baseline_binary}
    rm -rf "${srcdir}"
  fi
  ${baseline_binary} ${strobealign_options} ${ref} ${reads[@]} | samtools view -o ${baseline_bam}
fi

# Run strobealign. This recompiles from scratch to ensure we use consistent
# compiler options.
builddir=$(mktemp -p . -d build.XXXXXXX)
cmake . -B ${builddir} ${cmake_options}
make -s -j 4 -C ${builddir} strobealign
set -x
${builddir}/strobealign ${strobealign_options} ${ref} ${reads[@]} | samtools view -o head.bam
rm -rf ${builddir}

# Do the actual comparison
tests/samdiff.py ${baseline_bam} head.bam
