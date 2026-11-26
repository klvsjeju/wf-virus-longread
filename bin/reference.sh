#!/usr/bin/env bash
set -euo pipefail
mkdir -p result/alignment result/varcalling
# Prefix caller VCFs and their existing indexes before publishing the result tree.
while IFS= read -r -d '' vcf; do
 name="${vcf##*/}"
 [[ "$name" == "${SAMPLE}_"* ]] && continue
 target="${vcf%/*}/${SAMPLE}_${name}"
 mv "$vcf" "$target"
 for suffix in .tbi .csi; do
  if [[ -f "${vcf}${suffix}" ]]; then mv "${vcf}${suffix}" "${target}${suffix}"; fi
 done
done < <(find result/varcalling -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) -print0)
bam="result/alignment/${SAMPLE}_sorted.bam"
secondary="$(samtools view -c -f 0x100 "$bam")"
supplementary="$(samtools view -c -f 0x800 "$bam")"
printf 'alignment_type\tcount\nsecondary\t%s\nsupplementary\t%s\n' "$secondary" "$supplementary" > result/alignment/filter_check.tsv
test "$secondary" = 0 && test "$supplementary" = 0
bash "$(dirname "$OPTIONS_RUNNER")/report_depth.sh" "$bam" "$SAMPLE" result/alignment "$MASK_DEPTH_OPTS_B64"
# QUAL/DP selection applies to the normalized VCF, regardless of caller FILTER status.
variant_filter="QUAL!=\".\" && FORMAT/DP!=\".\" && QUAL>=${MIN_VARIANT_QUAL:-20} && FORMAT/DP>=${MIN_VARIANT_DEPTH:-30}"
bcftools norm -f "$REF" -m -any "result/varcalling/${SAMPLE}_merge_output.vcf.gz" | bcftools sort -Ou | \
 bcftools view -i "$variant_filter" -Oz -o "result/varcalling/${SAMPLE}_norm.vcf.gz"
bcftools index -f "result/varcalling/${SAMPLE}_norm.vcf.gz"
# Consensus additionally requires caller PASS and an inclusive allele-fraction threshold.
consensus_filter="FORMAT/AF!=\".\" && FORMAT/AF>=${MIN_VARIANT_AF:-0.5}"
bcftools view -f PASS -i "$consensus_filter" "result/varcalling/${SAMPLE}_norm.vcf.gz" -Oz -o "result/varcalling/${SAMPLE}_norm.cons.vcf.gz"
bcftools index -f "result/varcalling/${SAMPLE}_norm.cons.vcf.gz"
bcftools stats "result/varcalling/${SAMPLE}_norm.vcf.gz" > result/varcalling/variants.stats
mask=()
if (( MASK_DEPTH > 0 )); then
 awk -v d="$MASK_DEPTH" 'BEGIN {OFS="\t"} $3<d {print $1,$2-1,$2}' "result/alignment/${SAMPLE}.depth" > result/varcalling/low_coverage.bed
 mask=(-m result/varcalling/low_coverage.bed)
fi
python3 "$OPTIONS_RUNNER" consensus "$CONSENSUS_OPTS_B64" -- bcftools consensus -f "$REF" -s "$SAMPLE" "${mask[@]}" "result/varcalling/${SAMPLE}_norm.cons.vcf.gz" |
 awk -v s="$SAMPLE" '/^>/ {n++;print ">" s (n>1 ? "|contig" n : "");next} {print}' > "result/varcalling/${SAMPLE}_consensus.fa"
# Calling has finished; omit pileup/full-alignment intermediates from published copies.
find result/varcalling -type f \( \
 -name "${SAMPLE}_pileup.vcf*" -o -name "${SAMPLE}_pileup_*.vcf*" -o \
 -name "${SAMPLE}_full_alignment.vcf*" -o -name "${SAMPLE}_full_alignment_*.vcf*" \
 \) -delete
