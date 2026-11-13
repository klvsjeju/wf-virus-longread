nextflow.enable.dsl = 2

process INGEST {
 label 'wf_common'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}/metadata" }, mode: 'copy', pattern: '{basecaller.json,fastcat_stats}'
 input:
 tuple val(sample), path(reads, stageAs: 'chunks??/*')
 output:
 tuple val(sample), path('reads.fastq.gz'), path('basecaller.json'), emit: reads
 tuple val(sample), path('fastcat_stats'), emit: stats
 script:
 def read_files = reads instanceof List ? reads : [reads]
 def min_length = params.reference_based_assembly ? 0 : params.min_read_length
 def length_filter = min_length ? "-a ${min_length}" : ''
 def filtered_metadata = min_length ? '--filtered-reads reads.fastq.gz' : ''
 def check_reads = min_length ? "python3 -c \"import gzip,sys; sys.exit(0 if gzip.open('reads.fastq.gz','rb').read(1) else 'No reads meet min_read_length=${min_length} bp; lower --min_read_length or use 0.')\"" : ''
 """
 mkdir fastcat_stats
 fastcat ${length_filter} -s '${sample}' -f fastcat_stats/per-file-stats.tsv \
   -i fastcat_stats/per-file-runids.tsv -l fastcat_stats/per-file-basecallers.tsv \
   --histograms histograms ${read_files.collect { "'${it}'" }.join(' ')} | bgzip -@ ${task.cpus} > reads.fastq.gz
 mv histograms/* fastcat_stats/
 ${check_reads}
 python3 '${projectDir}/bin/inputs.py' fastcat-metadata --stats fastcat_stats --metadata basecaller.json \
   --override '${params.override_basecaller_cfg ?: ''}' ${filtered_metadata} ${read_files.collect { "'${it}'" }.join(' ')}
 """
}

process QC {
 label 'nanoplot'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}/qc" }, mode: 'copy'
 input:
 tuple val(sample), path(reads), path(metadata)
 output:
 tuple val(sample), path('NanoPlot')
 script:
 """
 NanoPlot --fastq '${reads}' -t ${task.cpus} -o NanoPlot -p '${sample}_' --tsv_stats --only-report
 """
}

process INDEX_REFERENCE {
 label 'wf_bacterial'
 input:
 path source, stageAs: 'input.fa'
 output:
 tuple path('reference.fa'), path('reference.fa.fai'), path('reference.fa.map-ont.mmi')
 script:
 """
 cp input.fa reference.fa
 mini_align -X -r reference.fa -t ${task.cpus}
 """
}

process PREPARE_CLAIR3_MODEL {
 label 'clair3'
 tag { basecaller }
 input:
 val basecaller
 path custom_root, stageAs: 'custom_source'
 path mapping, stageAs: 'mapping.tsv'
 path catalog, stageAs: 'catalog.json'
 output:
 tuple val(basecaller), path('clair3_model'), path('model_selection.json')
 script:
 def custom = params.model_dir ? "--custom-root '${custom_root}'" : ''
 """
 python3 '${projectDir}/bin/clair3_models.py' --basecaller '${basecaller}' --catalog '${catalog}' --mapping '${mapping}' ${custom}
 """
}

process ALIGN_READS {
 label 'wf_bacterial'
 tag { sample }
 input:
 tuple val(sample), path(reads), path(metadata)
 tuple path(ref), path(fai), path(mmi)
 output:
 tuple val(sample), path('alignment')
 script:
 """
 mkdir alignment
 mini_align -i '${reads}' -r '${ref}' -p raw -t ${task.cpus} -m
 samtools view -b -F 0x904 -q ${params.min_mapq} -o 'alignment/${sample}_sorted.bam' raw.bam
 samtools index -@ ${task.cpus} 'alignment/${sample}_sorted.bam'
 samtools flagstat 'alignment/${sample}_sorted.bam' > 'alignment/${sample}_flagstat.txt'
 """
}

process BAM_COVERAGE {
 label 'deeptools'
 tag { sample }
 input:
 tuple val(sample), path(alignment)
 output:
 tuple val(sample), path("${sample}.bedgraph")
 script:
 """
 bamCoverage -b '${alignment}/${sample}_sorted.bam' -o '${sample}.bedgraph' --outFileFormat bedgraph -p ${task.cpus}
 """
}

process CALL_CLAIR3 {
 label 'clair3'
 tag { sample }
 input:
 tuple val(sample), path(alignment), path(model_root, stageAs: 'models'), path(provenance)
 tuple path(ref), path(fai), path(mmi)
 output:
 tuple val(sample), path('clair3'), path('model_selection.json', includeInputs: true)
 script:
 """
 run_clair3.sh --bam_fn="\$PWD/${alignment}/${sample}_sorted.bam" --ref_fn="\$PWD/${ref}" \
   --model_path="\$(realpath '${model_root}')" --sample_name='${sample}' --threads=${task.cpus} \
   --platform=ont --output="\$PWD/clair3" --chunk_size=25000 --include_all_ctgs \
   --haploid_sensitive --no_phasing_for_fa --min_coverage=${params.min_coverage} \
   --enable_long_indel --remove_intermediate_dir
 """
}

process REFERENCE_ANALYSIS {
 label 'wf_bacterial'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}" }, mode: 'copy', enabled: !params.run_annotation,
   saveAs: { name -> name == 'result' ? null : name.replaceFirst(/^result\//, '') }
 input:
 tuple val(sample), path(alignment), path(coverage), path(clair3), path(provenance)
 tuple path(ref), path(fai), path(mmi)
 output:
 tuple val(sample), path('result/*', type: 'any'), emit: results
 script:
 """
 mkdir result
 cp -aL --reflink=auto '${alignment}' '${clair3}' '${provenance}' result/
 cp '${coverage}' result/alignment/
 export SAMPLE='${sample}' REF='${ref}' MASK_DEPTH=${params.mask_depth}
 export OPTIONS_RUNNER='${projectDir}/bin/run_options.py'
 export MASK_DEPTH_OPTS_B64='${(params.mask_depth_opts ?: '').toString().bytes.encodeBase64().toString()}'
 export CONSENSUS_OPTS_B64='${(params.consensus_opts ?: '').toString().bytes.encodeBase64().toString()}'
 bash '${projectDir}/bin/reference.sh'
 """
}

process SNPEFF_ANNOTATE {
 label 'snpeff'
 tag { sample }
 input:
 tuple val(sample), path(variants)
 path config, stageAs: 'snpEff.config'
 path data, stageAs: 'snpeff_data'
 output:
 tuple val(sample), path('annot.vcf')
 script:
 """
 snpEff ann -c snpEff.config -dataDir "\$PWD/snpeff_data" -ud 0 -no-intergenic '${params.snpeff_db}' '${variants}' > annot.vcf
 """
}

process SNPSIFT_TABLE {
 label 'snpsift'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}" }, mode: 'copy',
   saveAs: { name -> name == 'result' ? null : name.replaceFirst(/^result\//, '') }
 input:
 tuple val(sample), path(source, stageAs: 'source_result/*'), path(annotated)
 output:
 tuple val(sample), path('result/*', type: 'any'), emit: results
 script:
 """
 cp -aL --reflink=auto source_result result
 mkdir -p result/annot
 cp '${annotated}' 'result/annot/${sample}_annot.vcf'
 SnpSift extractFields -s ';' -e '.' 'result/annot/${sample}_annot.vcf' POS REF ALT \
   'ANN[*].GENE' 'ANN[*].IMPACT' 'ANN[*].EFFECT' 'GEN[0].GT' 'GEN[0].DP' 'GEN[0].AF' 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' | \
   awk -v s='${sample}' 'BEGIN {OFS="\t"} NR==1 {sub(/^#/, ""); print "Sample",\$0; next} {print s,\$0}' > result/annot/var.tsv
 """
}

process DE_NOVO_CONTIGS {
 label 'wf_bacterial'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}" }, mode: 'copy'
 input:
 tuple val(sample), path(reads), path(metadata)
 output:
 tuple val(sample), path("${sample}.contigs.fasta"), path('flye')
 script:
 """
 python3 '${projectDir}/bin/run_options.py' flye '${(params.flye_opts ?: '').toString().bytes.encodeBase64().toString()}' -- flye --nano-hq '${reads}' --genome-size ${params.flye_genome_size} --threads ${task.cpus} --out-dir flye
 test -s flye/assembly.fasta
 cp flye/assembly.fasta '${sample}.contigs.fasta'
 """
}

process DE_NOVO_COVERAGE {
 label 'wf_bacterial'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}" }, mode: 'copy'
 input:
 tuple val(sample), path(reads), path(metadata), path(contigs), path(flye)
 output:
 tuple val(sample), path('coverage')
 script:
 """
 cp '${contigs}' assembly.fa
 mini_align -i '${reads}' -r assembly.fa -p mapped -t ${task.cpus} -m
 samtools view -b -F 0x904 -q ${params.min_mapq} -o reporting.bam mapped.bam
 samtools index reporting.bam
 export OPTIONS_RUNNER='${projectDir}/bin/run_options.py'
 bash '${projectDir}/bin/report_depth.sh' reporting.bam '${sample}' coverage
 """
}

process SAMPLE_SUMMARY {
 label 'wf_common'
 tag { sample }
 input:
 tuple val(sample), path(metadata), path(qc), path(analysis, stageAs: 'analysis/*'), path(stats), path(flye), path(coverage)
 output:
 path 'summary.json'
 script:
 def mode = params.reference_based_assembly ? 'reference' : 'denovo'
 def analysis_path = params.reference_based_assembly ? 'analysis' : analysis
 """
 python3 '${projectDir}/bin/report.py' sample '${sample}' '${mode}' '${metadata}' '${analysis_path}' --stats '${stats}' --flye '${flye}' --coverage '${coverage}' --mask-depth ${params.mask_depth}
 """
}

process WORKFLOW_REPORT {
 label 'wf_common'
 publishDir params.out_dir, mode: 'copy'
 input:
 path summaries, stageAs: 'samples??/*'
 val missing_samples
 val sample_details
 output:
 path 'wf-virus-longread-report.html'
 path 'results.json'
 script:
 def settings = [workflow_version: workflow.manifest.version,
   missing_samples: missing_samples,
   sample_metadata: sample_details,
   run_name: workflow.runName, session_id: workflow.sessionId.toString(),
   reference_based_assembly: params.reference_based_assembly, reference: params.reference,
   run_annotation: params.run_annotation,
   min_mapq: params.min_mapq, min_coverage: params.min_coverage,
   mask_depth: params.mask_depth, mask_depth_opts: params.mask_depth_opts,
   consensus_opts: params.consensus_opts, flye_opts: params.flye_opts,
   min_read_length: params.min_read_length,
   effective_min_read_length: params.reference_based_assembly ? 0 : params.min_read_length,
   flye_read_type: 'nano-hq', flye_genome_size: params.flye_genome_size,
   override_basecaller_cfg: params.override_basecaller_cfg,
   container_override: params.container, containers: params.images, threads: params.threads]
 def summary_files = summaries instanceof List ? summaries : [summaries]
 def encoded = groovy.json.JsonOutput.toJson(settings).bytes.encodeBase64().toString()
 """
 python3 '${projectDir}/bin/report.py' aggregate '${encoded}' ${summary_files.collect { "'${it}'" }.join(' ')}
 """
}

workflow {
 if (params.help) { log.info 'See README.md. Required: --fastq. Default mode: de novo contigs only.'; return }
 if (!params.fastq) error 'Select --fastq (FASTQ file or directory)'
 if (!params.out_dir) error 'Unable to resolve output directory; specify --out_dir'
 ['threads','flye_genome_size','min_mapq','min_coverage','mask_depth','min_read_length'].each { k ->
   if (!(params[k].toString() ==~ /[0-9]+/)) error "--${k} must be an integer"
 }
 if ((params.threads as int) < 1 || (params.flye_genome_size as int) < 1) error 'threads and genome size must be positive'
 if (params.override_basecaller_cfg && !(params.override_basecaller_cfg ==~ /[A-Za-z0-9_.@+-]+/)) error 'Invalid basecaller name'
 def root = file(params.fastq, checkIfExists: true)
 def isFq = { p -> p.isFile() && (p.name ==~ /.*\.(fastq|fq)(\.gz)?/) }
 def samples = []
 def multiplex = false
 if (isFq(root)) samples << [SampleNames.singleName(params.sample, root.simpleName), [root]]
 else if (root.isDirectory()) {
   def direct = root.listFiles().findAll { isFq(it) }.sort { it.name }
   def dirs = root.listFiles().findAll { it.isDirectory() && it.toAbsolutePath().normalize() != java.nio.file.Paths.get(params.out_dir.toString()).toAbsolutePath().normalize() }.sort { it.name }.collect { dir ->
     [dir.baseName, dir.listFiles().findAll { isFq(it) }.sort { it.name }]
   }.findAll { it[1] && (params.analyse_unclassified || it[0] != 'unclassified') }
   if (direct && dirs) error 'Do not mix top-level FASTQ and sample folders'
   if (direct) samples << [SampleNames.singleName(params.sample, root.baseName), direct]
   else {
     samples.addAll(dirs)
     multiplex = true
   }
 }
 if (!samples) error 'No FASTQ files found (supports one subdirectory level)'
 def missing_samples = Channel.value([])
 def sample_details = Channel.value(samples.collectEntries { [(it[0]): [barcode: null]] })
 if (multiplex) {
   def original = samples
   def records = params.sample_sheet ?
     Channel.fromPath(params.sample_sheet, checkIfExists: true).splitCsv(header: true, quote: '"').toList() :
     Channel.value([])
   def selection = records.map { sheet ->
     def result = SampleNames.selectDirectory(original, params.sample_sheet ? sheet : null, params.sample)
     result.warnings.each { log.warn it }
     if (!result.samples && !result.missing) error 'No input samples match the selected sample name'
     result
   }
   samples = selection.flatMap { it.samples }
   missing_samples = selection.map { it.missing }
   sample_details = selection.map { it.metadata }
 } else samples = Channel.fromList(samples)
 def validated = samples.map { sample, reads ->
   if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Invalid sample name: ${sample}"
   if (sample.toLowerCase() in ['execution', 'results.json', 'wf-virus-longread-report.html']) {
     error "Sample name ${sample} conflicts with a shared workflow output; choose another sample name or alias"
   }
   reads.each { p -> if (p.toString().contains("'") || p.toString().contains('\n')) error 'Quotes/newlines in paths are unsupported' }
   tuple(sample, reads)
 }.toList().flatMap { rows ->
   if (rows.collect { it[0] }.toSet().size() != rows.size()) error 'Duplicate sample names'
   rows
 }
 INGEST(validated)
 QC(INGEST.out.reads)
 def analysis_results
 if (params.reference_based_assembly) {
   ['reference'].each { k -> if (!params[k]) error "Reference mode requires --${k}" }
   if (params.run_annotation) {
     ['snpeff_config','snpeff_data','snpeff_db'].each { k -> if (!params[k]) error "Annotation requires --${k}" }
     if (!(params.snpeff_db ==~ /[A-Za-z0-9_.-]+/)) error 'Invalid snpEff database ID'
   }
   def ref = file(params.reference, checkIfExists: true)
   if (params.model_map && !params.model_dir) error '--model_map requires --model_dir; leave both unset for automatic models'
   def models = file(params.model_dir ?: "${projectDir}/model_map.tsv", checkIfExists: true)
   def mapping = file(params.model_map ?: "${projectDir}/model_map.tsv", checkIfExists: true)
   if (!ref.isFile()) error 'reference must be a FASTA file'
   if (params.model_dir && !models.isDirectory()) error 'model_dir must be a directory'
   INDEX_REFERENCE(ref)
   def keyed_reads = INGEST.out.reads.map { sample, reads, metadata ->
     def basecaller = new groovy.json.JsonSlurper().parseText(metadata.text).selected ?: 'unknown'
     if (!(basecaller ==~ /[A-Za-z0-9_.@+-]+/)) error 'Invalid basecaller name in FASTQ metadata'
     tuple(basecaller, sample, reads, metadata)
   }
   PREPARE_CLAIR3_MODEL(keyed_reads.map { key, sample, reads, metadata -> key }.unique(), models, mapping,
     file("${projectDir}/clair3_model_catalog.json", checkIfExists: true))
   def sample_models = keyed_reads.map { key, sample, reads, metadata -> tuple(key, sample) }
     .combine(PREPARE_CLAIR3_MODEL.out, by: 0)
     .map { key, sample, checkpoint, provenance -> tuple(sample, checkpoint, provenance) }
   ALIGN_READS(INGEST.out.reads, INDEX_REFERENCE.out)
   BAM_COVERAGE(ALIGN_READS.out)
   CALL_CLAIR3(ALIGN_READS.out.join(sample_models), INDEX_REFERENCE.out)
   REFERENCE_ANALYSIS(ALIGN_READS.out.join(BAM_COVERAGE.out).join(CALL_CLAIR3.out), INDEX_REFERENCE.out)
   // Reference mode has no Flye directory; use an existing inert file for the path input.
   analysis_results = REFERENCE_ANALYSIS.out.results.map { sample, result -> tuple(sample, result, file("${projectDir}/model_map.tsv"), file("${projectDir}/clair3_model_catalog.json")) }
   if (params.run_annotation) {
     def config = file(params.snpeff_config, checkIfExists: true)
     def data = file(params.snpeff_data, checkIfExists: true)
     if (!data.isDirectory()) error 'snpeff_data must be a directory'
     if (!file("${data}/${params.snpeff_db}/snpEffectPredictor.bin").exists()) error 'snpEff DB binary not found'
     SNPEFF_ANNOTATE(REFERENCE_ANALYSIS.out.results.map { sample, result ->
       tuple(sample, result.find { it.name == 'clair3' }.resolve("${sample}_norm.vcf.gz"))
     }, config, data)
     SNPSIFT_TABLE(REFERENCE_ANALYSIS.out.results.join(SNPEFF_ANNOTATE.out))
     analysis_results = SNPSIFT_TABLE.out.results.map { sample, result -> tuple(sample, result, file("${projectDir}/model_map.tsv"), file("${projectDir}/clair3_model_catalog.json")) }
   }
 } else {
   DE_NOVO_CONTIGS(INGEST.out.reads)
   DE_NOVO_COVERAGE(INGEST.out.reads.join(DE_NOVO_CONTIGS.out))
   analysis_results = DE_NOVO_CONTIGS.out.join(DE_NOVO_COVERAGE.out)
 }
 metadata_results = INGEST.out.reads.map { sample, reads, metadata -> tuple(sample, metadata) }
 SAMPLE_SUMMARY(metadata_results.join(QC.out).join(analysis_results).join(INGEST.out.stats)
   .map { sample, metadata, qc, analysis, flye, coverage, stats -> tuple(sample, metadata, qc, analysis, stats, flye, coverage) })
 WORKFLOW_REPORT(SAMPLE_SUMMARY.out.toList(), missing_samples, sample_details)
}
