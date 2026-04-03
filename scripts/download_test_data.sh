#!/usr/bin/env bash
set -euo pipefail

ILLUMINA_ACC="ERR10879100"
ONT_ACC="ERR10879072"
SAMPLE_ID="sample01"

BASE_DIR="$(pwd)/data"

ILLUMINA_DIR="${BASE_DIR}/illumina/${SAMPLE_ID}"
ONT_DIR="${BASE_DIR}/ont/${SAMPLE_ID}"
SRA_CACHE="${BASE_DIR}/sra_cache"

mkdir -p "$ILLUMINA_DIR" "$ONT_DIR" "$SRA_CACHE"

echo "Downloading Illumina paired-end data: $ILLUMINA_ACC"
docker run --rm \
  -v "$ILLUMINA_DIR:/out" \
  -v "$SRA_CACHE:/ncbi/public/sra" \
  ncbi/sra-tools \
  fastq-dump "$ILLUMINA_ACC" \
    --split-files \
    --gzip \
    -O /out

echo "Downloading ONT data: $ONT_ACC"
docker run --rm \
  -v "$ONT_DIR:/out" \
  -v "$SRA_CACHE:/ncbi/public/sra" \
  ncbi/sra-tools \
  fastq-dump "$ONT_ACC" \
    --gzip \
    -O /out

echo "Download complete."