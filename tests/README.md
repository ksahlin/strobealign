# Tiny test dataset

- `phix.fasta` is the complete sequence of NC_001422.1 (phiX-174)
- `phix.{1,2}.fastq` are the first 100 reads of SRR1377138:

       fastq-dump --split-3 -X 100 --defline-seq '@$ac.$si' --defline-qual '+' SRR1377138
