#!/usr/bin/env bash
set -euo pipefail
mkdir -p result/alignment result/clair3
bam="result/alignment/${SAMPLE}_sorted.bam"
secondary="$(samtools view -c -f 0x100 "$bam")"
supplementary="$(samtools view -c -f 0x800 "$bam")"
printf 'alignment_type\tcount\nsecondary\t%s\nsupplementary\t%s\n' "$secondary" "$supplementary" > result/alignment/filter_check.tsv
test "$secondary" = 0 && test "$supplementary" = 0
python3 "$OPTIONS_RUNNER" depth "$MASK_DEPTH_OPTS_B64" -- samtools depth -aa -J "$bam" > "result/alignment/${SAMPLE}.depth"
bcftools norm -f "$REF" -m -any result/clair3/merge_output.vcf.gz | bcftools sort -Oz -o result/clair3/norm.vcf.gz
bcftools index -f result/clair3/norm.vcf.gz
bcftools view -f PASS result/clair3/norm.vcf.gz -Oz -o result/clair3/norm.pass.vcf.gz
bcftools index -f result/clair3/norm.pass.vcf.gz
mask=()
if (( MASK_DEPTH > 0 )); then
 awk -v d="$MASK_DEPTH" 'BEGIN {OFS="\t"} $3<d {print $1,$2-1,$2}' "result/alignment/${SAMPLE}.depth" > result/clair3/low_coverage.bed
 mask=(-m result/clair3/low_coverage.bed)
fi
python3 "$OPTIONS_RUNNER" consensus "$CONSENSUS_OPTS_B64" -- bcftools consensus -f "$REF" -s "$SAMPLE" "${mask[@]}" result/clair3/norm.pass.vcf.gz |
 awk -v s="$SAMPLE" '/^>/ {n++;print ">" s "|contig" n;next} {print}' > "result/clair3/${SAMPLE}_consensus.fa"
