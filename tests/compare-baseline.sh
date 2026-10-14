#!/bin/bash
#
# Run strobealign on some test data and compare the output against what the
# baseline version produces.
#
# - Test data is automatically downloaded if needed
#   (and put into tests/drosophila)
# - The baseline BAM or PAF file is generated if necessary

set -euo pipefail

# Fail early if pysam is missing
python3 -c 'import pysam'

ends="pe"
threads=4
baseline_commit=$(git --no-pager log -n1 --pretty=format:%H --grep='^Is-new-baseline: yes')

mode=align
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
    m)
      mode=map
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

if [[ ${mode} = align ]]; then ext=bam; else ext=paf.gz; fi

# Ensure test data is available
tests/download.sh

baseline_binary=baseline/strobealign-${baseline_commit}
baseline_file=baseline/bampaf/${baseline_commit}.${ends}.${ext}

# Generate the baseline BAM if necessary
mkdir -p baseline/bampaf
if ! test -f ${baseline_file}; then
  if ! test -f ${baseline_binary}; then
    srcdir=$(mktemp -p . -d compile.XXXXXXX)
    git clone . ${srcdir}
    pushd ${srcdir}
    git checkout -d ${baseline_commit}
    cargo build
    popd
    mv ${srcdir}/target/debug/strobealign ${baseline_binary}
    rm -rf "${srcdir}"
  fi
  if [[ ${mode} = align ]]; then
    ${baseline_binary} -N 2 -v -t ${threads} ${ref} ${reads[@]} | samtools view -o ${baseline_file}.tmp.${ext}
  else
    ${baseline_binary} -N 2 -v -x -t ${threads} ${ref} ${reads[@]} | gzip > ${baseline_file}.tmp.${ext}
  fi
  mv ${baseline_file}.tmp.${ext} ${baseline_file}
fi

# Build and run strobealign
cargo build
set -x
if [[ ${mode} = align ]]; then
  RUST_BACKTRACE=1 target/debug/strobealign -N 2 -v -t ${threads} ${ref} ${reads[@]} | samtools view -o head.bam
  tests/samdiff.py ${baseline_file} head.bam
else
  RUST_BACKTRACE=1 target/debug/strobealign -N 2 -v -x -t ${threads} ${ref} ${reads[@]} | gzip > head.paf.gz
  diff -u <(zcat ${baseline_file}) <(zcat head.paf.gz) | head -n 100
fi
