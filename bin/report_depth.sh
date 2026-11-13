#!/usr/bin/env bash
# Generate total and strand-specific depth with identical samtools options.
set -euo pipefail
bam=$1
sample=$2
destination=$3
options=${4:-}
mkdir -p "$destination"
python3 "$OPTIONS_RUNNER" depth "$options" -- samtools depth -aa -J "$bam" > "$destination/$sample.depth"
samtools view -b -F 16 "$bam" | python3 "$OPTIONS_RUNNER" depth "$options" -- samtools depth -aa -J - > "$destination/$sample.fwd.depth"
samtools view -b -f 16 "$bam" | python3 "$OPTIONS_RUNNER" depth "$options" -- samtools depth -aa -J - > "$destination/$sample.rev.depth"
