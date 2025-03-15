#!/bin/bash

export LC_ALL="C"
set -euo pipefail

ln -sf /data/nano/nfNeST_nano.nf /data/
./nextflow run nfNeST_nano.nf -c /data/nano/nextflow.config -with-report /data/output/test_output.html
