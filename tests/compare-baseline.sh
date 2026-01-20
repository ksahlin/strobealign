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
baseline_commit=$(git --no-pager log -n1 --pretty=format:%H --grep='^Is-new-baseline: yes')

while getopts "b:st:m" opt; do
  case "${opt}" in
    b)
      baseline_commit=$(git rev-parse "${OPTARG}")
      ;;
    t)
      threads=$OPTARG
      ;;
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

baseline_binary=baseline/strobealign-${baseline_commit}
extra_ext=""
baseline_bam=baseline/bam/${baseline_commit}.${ends}${extra_ext}.bam

# Generate the baseline BAM if necessary
mkdir -p baseline/bam
if ! test -f ${baseline_bam}; then
  if ! test -f ${baseline_binary}; then
    srcdir=$(mktemp -p . -d compile.XXXXXXX)
    git clone . ${srcdir}
    pushd ${srcdir}
    git checkout -d ${baseline_commit}
    cargo build --release
    popd
    mv ${srcdir}/target/release/strobealign ${baseline_binary}
    rm -rf "${srcdir}"
  fi
  ${baseline_binary} -v -t ${threads} ${ref} ${reads[@]} | samtools view -o ${baseline_bam}.tmp.bam
  mv ${baseline_bam}.tmp.bam ${baseline_bam}
fi

# Build and run strobealign
cargo build --release
set -x
target/release/strobealign -v -t ${threads} ${ref} ${reads[@]} | samtools view -o head.bam

# Do the actual comparison
tests/samdiff.py ${baseline_bam} head.bam
