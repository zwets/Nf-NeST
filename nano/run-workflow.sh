#!/bin/bash

export LC_ALL="C"
set -euo pipefail

cd $(dirname $0)

# Set to the absolute directory path containing the barcode* directories
INPUT_DIR=/project/bright-malaria/data/fastq_pass

# We write output to this directory (must be absolute)
OUTPUT_DIR="$PWD/output"

mkdir -p "$OUTPUT_DIR" &&
exec docker run --rm -ti \
    -v "$INPUT_DIR:/data/inputs:ro" \
    -v "$OUTPUT_DIR:/data/output" \
    -v "$PWD:/data/nano:ro" \
    -v "$(realpath "../pyscripts"):/data/pyscripts" \
    supark87/nfnest:latest \
    /data/nano/helper-script.sh
