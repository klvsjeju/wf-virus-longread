nextflow.enable.dsl = 2

process INGEST {
 tag sample
 publishDir "${params.out_dir}/metadata/${sample}", mode: 'copy', pattern: 'basecaller.json'
 input:
 tuple val(sample), path(reads, stageAs: 'chunks??/*', arity: '1..*')
 output:
 tuple val(sample), path('reads.fastq.gz'), path('basecaller.json')
 script:
 """
 python3 '${projectDir}/bin/inputs.py' merge --output reads.fastq.gz --metadata basecaller.json --override '${params.override_basecaller_cfg ?: ''}' ${reads.collect { "'${it}'" }.join(' ')}
 """
}

process QC {
 tag sample
 publishDir "${params.out_dir}/qc/${sample}", mode: 'copy'
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
 input:
 path source, stageAs: 'input.fa'
 output:
 tuple path('reference.fa'), path('reference.fa.fai'), path('reference.mmi')
 script:
 """
 cp input.fa reference.fa
 samtools faidx reference.fa
 minimap2 -x map-ont -d reference.mmi reference.fa
 """
}

process REFERENCE_ANALYSIS {
 tag sample
 publishDir "${params.out_dir}/reference/${sample}", mode: 'copy', pattern: 'result'
 input:
 tuple val(sample), path(reads), path(metadata)
 tuple path(ref), path(fai), path(mmi)
 path model_root, stageAs: 'models'
 path model_map, stageAs: 'mapping.tsv'
 path config, stageAs: 'snpEff.config'
 path data, stageAs: 'snpeff_data'
 output:
 tuple val(sample), path('result')
 script:
 """
 python3 '${projectDir}/bin/inputs.py' model --root models --mapping mapping.tsv --metadata '${metadata}' > selected_model.txt
 export SAMPLE='${sample}' READS='${reads}' REF='${ref}' MMI='${mmi}'
 export THREADS=${task.cpus} MAPQ=${params.min_mapq} MIN_COVERAGE=${params.min_coverage} MASK_DEPTH=${params.mask_depth}
 export OPTIONS_RUNNER='${projectDir}/bin/run_options.py'
 export MASK_DEPTH_OPTS_B64='${(params.mask_depth_opts ?: '').toString().bytes.encodeBase64().toString()}'
 export CONSENSUS_OPTS_B64='${(params.consensus_opts ?: '').toString().bytes.encodeBase64().toString()}'
 export SNPEFF_DB='${params.snpeff_db}'
 bash '${projectDir}/bin/reference.sh'
 """
}

process DE_NOVO_CONTIGS {
 tag sample
 publishDir "${params.out_dir}/denovo/${sample}", mode: 'copy'
 input:
 tuple val(sample), path(reads), path(metadata)
 output:
 tuple val(sample), path("${sample}.contigs.fasta"), path('flye')
 script:
 """
 python3 '${projectDir}/bin/run_options.py' flye '${(params.flye_opts ?: '').toString().bytes.encodeBase64().toString()}' -- flye --${params.flye_read_type} '${reads}' --genome-size ${params.flye_genome_size} --threads ${task.cpus} --out-dir flye
 test -s flye/assembly.fasta
 cp flye/assembly.fasta '${sample}.contigs.fasta'
 """
}

process SAMPLE_SUMMARY {
 tag sample
 input:
 tuple val(sample), path(metadata), path(qc), path(analysis)
 output:
 path 'summary.json'
 script:
 def mode = params.reference_based_assembly ? 'reference' : 'denovo'
 """
 python3 '${projectDir}/bin/report.py' sample '${sample}' '${mode}' '${metadata}' '${analysis}'
 """
}

process WORKFLOW_REPORT {
 publishDir params.out_dir, mode: 'copy'
 input:
 path summaries, stageAs: 'samples??/*'
 output:
 path 'wf-vzv-report.html'
 path 'results.json'
 script:
 def settings = [workflow_version: workflow.manifest.version,
   run_name: workflow.runName, session_id: workflow.sessionId.toString(),
   reference_based_assembly: params.reference_based_assembly, reference: params.reference,
   min_mapq: params.min_mapq, min_coverage: params.min_coverage,
   mask_depth: params.mask_depth, mask_depth_opts: params.mask_depth_opts,
   consensus_opts: params.consensus_opts, flye_opts: params.flye_opts,
   flye_read_type: params.flye_read_type, flye_genome_size: params.flye_genome_size,
   override_basecaller_cfg: params.override_basecaller_cfg,
   container: params.container, threads: params.threads]
 def encoded = groovy.json.JsonOutput.toJson(settings).bytes.encodeBase64().toString()
 """
 python3 '${projectDir}/bin/report.py' aggregate '${encoded}' ${summaries.collect { "'${it}'" }.join(' ')}
 """
}

workflow {
 if (params.help) { log.info 'See README.md. Required: --fastq. Default mode: de novo contigs only.'; return }
 if (!params.fastq) error 'Select --fastq (FASTQ file or directory)'
 ['threads','flye_genome_size','min_mapq','min_coverage','mask_depth'].each { k ->
   if (!(params[k].toString() ==~ /[0-9]+/)) error "--${k} must be an integer"
 }
 if ((params.threads as int) < 1 || (params.flye_genome_size as int) < 1) error 'threads and genome size must be positive'
 if (!(params.flye_read_type in ['nano-hq','nano-raw'])) error 'Invalid flye_read_type'
 if (params.override_basecaller_cfg && !(params.override_basecaller_cfg ==~ /[A-Za-z0-9_.@+-]+/)) error 'Invalid basecaller name'
 def root = file(params.fastq, checkIfExists: true)
 def isFq = { p -> p.isFile() && (p.name ==~ /.*\.(fastq|fq)(\.gz)?/) }
 def samples = []
 if (isFq(root)) samples << [params.sample, [root]]
 else if (root.isDirectory()) {
   def direct = root.listFiles().findAll { isFq(it) }.sort { it.name }
   def dirs = root.listFiles().findAll { it.isDirectory() }.sort { it.name }.collect { dir ->
     [dir.name, dir.listFiles().findAll { isFq(it) }.sort { it.name }]
   }.findAll { it[1] && (params.analyse_unclassified || it[0] != 'unclassified') }
   if (direct && dirs) error 'Do not mix top-level FASTQ and sample folders'
   if (direct) samples << [params.sample, direct]
   else samples.addAll(dirs)
 }
 if (!samples) error 'No FASTQ files found (supports one subdirectory level)'
 if (params.sample_sheet) {
   if (!root.isDirectory() || root.listFiles().any { isFq(it) }) error 'sample_sheet requires sample/barcode subdirectories'
   def rows = Channel.fromPath(params.sample_sheet, checkIfExists: true).splitCsv(header: true).toList()
   // Perform CSV handling in a dataflow operator; keep sample ordering deterministic.
   def original = samples
   def renamed = rows.map { records ->
     if (!records || !records[0].containsKey('barcode') || !records[0].containsKey('alias')) error 'CSV needs barcode,alias columns'
     def barcodes = records.collect { it.barcode }
     if (barcodes.toSet().size() != barcodes.size()) error 'Duplicate barcode in CSV'
     if (barcodes.toSet() != original.collect { it[0] }.toSet()) error 'CSV must exactly match selected barcode folders'
     def mapping = records.collectEntries { [(it.barcode): it.alias] }
     original.collect { row -> [mapping[row[0]], row[1]] }
   }
   samples = renamed.flatMap { it }
 } else samples = Channel.fromList(samples)
 def validated = samples.map { sample, reads ->
   if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Invalid sample name: ${sample}"
   reads.each { p -> if (p.toString().contains("'") || p.toString().contains('\n')) error 'Quotes/newlines in paths are unsupported' }
   tuple(sample, reads)
 }.toList().flatMap { rows ->
   if (rows.collect { it[0] }.toSet().size() != rows.size()) error 'Duplicate sample names'
   rows
 }
 INGEST(validated)
 QC(INGEST.out)
 def analysis_results
 if (params.reference_based_assembly) {
   ['reference','model_dir','snpeff_config','snpeff_data','snpeff_db'].each { k -> if (!params[k]) error "Reference mode requires --${k}" }
   if (!(params.snpeff_db ==~ /[A-Za-z0-9_.-]+/)) error 'Invalid snpEff database ID'
   def ref = file(params.reference, checkIfExists: true)
   def models = file(params.model_dir, checkIfExists: true)
   def config = file(params.snpeff_config, checkIfExists: true)
   def data = file(params.snpeff_data, checkIfExists: true)
   def mapping = file(params.model_map ?: "${projectDir}/model_map.tsv", checkIfExists: true)
   if (!ref.isFile()) error 'reference must be a FASTA file'
   if (!models.isDirectory() || !data.isDirectory()) error 'model_dir and snpeff_data must be directories'
   if (!file("${data}/${params.snpeff_db}/snpEffectPredictor.bin").exists()) error 'snpEff DB binary not found'
   INDEX_REFERENCE(ref)
   REFERENCE_ANALYSIS(INGEST.out, INDEX_REFERENCE.out, models, mapping, config, data)
   analysis_results = REFERENCE_ANALYSIS.out
 } else {
   DE_NOVO_CONTIGS(INGEST.out)
   analysis_results = DE_NOVO_CONTIGS.out.map { sample, contigs, flye -> tuple(sample, contigs) }
 }
 metadata_results = INGEST.out.map { sample, reads, metadata -> tuple(sample, metadata) }
 SAMPLE_SUMMARY(metadata_results.join(QC.out).join(analysis_results))
 WORKFLOW_REPORT(SAMPLE_SUMMARY.out.collect())
}
