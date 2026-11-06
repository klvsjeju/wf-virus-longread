#!/usr/bin/env bash
set -euo pipefail
mkdir -p result/alignment result/clair3 result/annot
bam="result/alignment/${SAMPLE}_sorted.bam"
model_path="$(cat selected_model.txt)"
minimap2 -ax map-ont --secondary=no -t "$THREADS" "$MMI" "$READS" |
 samtools view -u -F 0x900 -q "$MAPQ" - |
 samtools sort -@ "$THREADS" -o "$bam" -
samtools index -@ "$THREADS" "$bam"
secondary="$(samtools view -c -f 0x100 "$bam")"
supplementary="$(samtools view -c -f 0x800 "$bam")"
printf 'alignment_type\tcount\nsecondary\t%s\nsupplementary\t%s\n' "$secondary" "$supplementary" > result/alignment/filter_check.tsv
test "$secondary" = 0 && test "$supplementary" = 0
samtools flagstat "$bam" > "result/alignment/${SAMPLE}_flagstat.txt"
bamCoverage -b "$bam" -o "result/alignment/${SAMPLE}.bedgraph" --outFileFormat bedgraph -p "$THREADS"
python3 "$OPTIONS_RUNNER" depth "$MASK_DEPTH_OPTS_B64" -- samtools depth -aa -J "$bam" > "result/alignment/${SAMPLE}.depth"
run_clair3.sh --bam_fn="$PWD/$bam" --ref_fn="$PWD/$REF" --model_path="$model_path" \
 --sample_name="$SAMPLE" --threads="$THREADS" --platform=ont --output="$PWD/result/clair3" \
 --chunk_size=25000 --include_all_ctgs --haploid_sensitive --no_phasing_for_fa \
 --min_coverage="$MIN_COVERAGE" --enable_long_indel --remove_intermediate_dir
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
snpEff ann -c snpEff.config -dataDir "$PWD/snpeff_data" -ud 0 -no-intergenic "$SNPEFF_DB" result/clair3/norm.vcf.gz > result/annot/annot.vcf
SnpSift extractFields -s ';' -e '.' result/annot/annot.vcf POS REF ALT \
 'ANN[*].GENE' 'ANN[*].IMPACT' 'ANN[*].EFFECT' 'GEN[0].GT' 'GEN[0].DP' 'GEN[0].AF' 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' |
 awk -v s="$SAMPLE" 'BEGIN {OFS="\t"} NR==1 {sub(/^#/,"",$0);print "Sample",$0;next} {print s,$0}' > result/annot/var.tsv
