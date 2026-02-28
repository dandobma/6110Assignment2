#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/qc

# Run FastQC on all FASTQ files
fastqc data/fastq/*.fastq.gz -o results/qc

# Summarize with MultiQC
multiqc results/qc -o results/qc
