#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/fastq

SRRS=(
  SRR10551665 SRR10551664 SRR10551663   # Early biofilm (IL20-IL22)
  SRR10551662 SRR10551661 SRR10551660   # Thin biofilm  (IL23-IL25)
  SRR10551659 SRR10551658 SRR10551657   # Mature biofilm (IL29-IL31)
)

for srr in "${SRRS[@]}"; do
  echo "=== Downloading ${srr} ==="
  fasterq-dump "${srr}" --split-files -O data/fastq
done

echo "=== Gzipping FASTQ files ==="
gzip -f data/fastq/*.fastq
