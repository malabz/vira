#!/usr/bin/env bash
set -euo pipefail

vira --help >/dev/null
command -v mafft >/dev/null

printf '>reference\nACGTACGTACGT\n' > reference.fasta
printf '>query1\nACGTACGTACGT\n' > queries.fasta

vira \
  --ref reference.fasta \
  --query queries.fasta \
  --out result.fasta \
  --max-instances 1 \
  --threads-per-mafft 1 \
  --chunk-size 1 \
  --enable-mapout-qc \
  --cleanup

test -s result.fasta
test -s result.qc_summary.tsv
test -s result.sample_qc_summary.tsv
test -s result.filtered.fasta
grep -q '^>reference' result.fasta
grep -q '^>query1' result.fasta
grep -q 'query1' result.sample_qc_summary.tsv
