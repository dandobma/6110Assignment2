#!/usr/bin/env bash
set -euo pipefail

INDEX="reference/salmon_index"
FASTQ_DIR="data/fastq"
OUT_DIR="results/salmon"

mkdir -p ${OUT_DIR}

for fq in ${FASTQ_DIR}/*.fastq.gz; do
  sample=$(basename ${fq} .fastq.gz)
  echo "Quantifying ${sample}"

  salmon quant \
    -i ${INDEX} \
    -l A \
    -r ${fq} \
    --validateMappings \
    --gcBias \
    --seqBias \
    --posBias \
    --fldMean 200 \
    --fldSD 20 \
    -o ${OUT_DIR}/${sample}
done
