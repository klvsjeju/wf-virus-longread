#!/usr/bin/env bash
# Probe the same container images used for analysis, without installing packages.
set -euo pipefail
output=$1
shift
printf 'Tool\tVersion\n' > "$output"
for tool in "$@"; do
 case "$tool" in
 fastcat|bgzip|samtools|bcftools|minimap2) command=("$tool" --version);;
 mini_align) command=(mini_align -h);;
 NanoPlot) command=(NanoPlot --version);;
 bamCoverage) command=(bamCoverage --version);;
 Flye) command=(flye --version);;
 Clair3) command=(run_clair3.sh --version);;
 snpEff) command=(snpEff -version);;
 SnpSift) command=(SnpSift -version);;
 Python) command=(python3 --version);;
 Java) command=(java -version);;
 *) printf '%s\tNot detected\n' "$tool" >> "$output"; continue;;
 esac
 # Some version/help commands return nonzero even when printing their version.
 raw=$(timeout 15 "${command[@]}" 2>&1 || true)
 version=$(printf '%s\n' "$raw" | sed -nE 's/.*[vV]ersion[[:space:]:=vV"]*([0-9]+\.[0-9]+[^[:space:]",;]*).*/\1/p' | head -n 1)
 if [[ -z "$version" ]]; then
   version=$(printf '%s\n' "$raw" | sed -nE 's/^([[:alnum:]_-]+[[:space:]]+)?(\([^)]*\)[[:space:]]+)?[vV]?([0-9]+\.[0-9]+[^[:space:]",;]*).*$/\3/p' | head -n 1)
 fi
 if [[ -z "$version" && "$tool" == mini_align && "$raw" == *'Align fastq'* ]]; then version='Unversioned alignment wrapper'; fi
 printf '%s\t%s\n' "$tool" "${version:-Not detected}" >> "$output"
done
