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
 tuple val(sample), path('varcalling'), path('model_selection.json', includeInputs: true)
 script:
 """
 run_clair3.sh --bam_fn="\$PWD/${alignment}/${sample}_sorted.bam" --ref_fn="\$PWD/${ref}" \
   --model_path="\$(realpath '${model_root}')" --sample_name='${sample}' --threads=${task.cpus} \
   --platform=ont --output="\$PWD/varcalling" --chunk_size=25000 --include_all_ctgs \
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
 tuple val(sample), path(alignment), path(coverage), path(varcalling), path(provenance)
 tuple path(ref), path(fai), path(mmi)
 output:
 tuple val(sample), path('result/*', type: 'any'), emit: results
 script:
 """
 mkdir result
 cp -aL --reflink=auto '${alignment}' '${varcalling}' '${provenance}' result/
 cp '${coverage}' result/alignment/
 export SAMPLE='${sample}' REF='${ref}' MASK_DEPTH=${params.mask_depth}
 export MIN_VARIANT_QUAL=${params.min_variant_qual} MIN_VARIANT_DEPTH=${params.min_variant_depth} MIN_VARIANT_AF=${params.min_variant_af}
 export OPTIONS_RUNNER='${projectDir}/bin/run_options.py'
 export MASK_DEPTH_OPTS_B64='${(params.mask_depth_opts ?: '').toString().bytes.encodeBase64().toString()}'
 export CONSENSUS_OPTS_B64='${(params.consensus_opts ?: '').toString().bytes.encodeBase64().toString()}'
 bash '${projectDir}/bin/reference.sh'
 """
}

process PREPARE_SNPEFF_DB {
 label 'snpeff'
 publishDir { "${params.out_dir}/annotation_database" }, mode: 'copy', pattern: 'prepared_db'
 input:
 path reference, stageAs: 'reference.fa'
 path gff, stageAs: 'annotation.gff'
 path cache_root, stageAs: 'db_cache'
 output:
 path 'prepared_db'
 script:
 def gff_arg = params.annotation_gff ? "--gff '${gff}'" : ''
 """
 python3 '${projectDir}/bin/snpeff_db.py' --reference '${reference}' ${gff_arg} --cache '${cache_root}' --output prepared_db
 """
}

process SNPEFF_ANNOTATE {
 label 'snpeff'
 tag { sample }
 input:
 tuple val(sample), path(variants)
 tuple path(config, stageAs: 'snpEff.config'), path(data, stageAs: 'snpeff_data'), val(db_id), val(provenance_b64)
 output:
 tuple val(sample), path('annot.vcf'), path('annotation_provenance.json')
 script:
 """
 snpEff ann -c snpEff.config -dataDir "\$PWD/snpeff_data" -ud 0 -no-intergenic '${db_id}' '${variants}' > annot.vcf
 printf '%s' '${provenance_b64}' | base64 -d > annotation_provenance.json
 """
}

process SNPSIFT_TABLE {
 label 'snpsift'
 tag { sample }
 input:
 tuple val(sample), path(source, stageAs: 'source_result/*'), path(annotated), path(annotation_provenance)
 output:
 tuple val(sample), path('result/*', type: 'any'), emit: results
 script:
 """
 cp -aL --reflink=auto source_result result
 mkdir -p result/annot
 cp '${annotation_provenance}' result/annot/annotation_provenance.json
 cp '${annotated}' 'result/annot/${sample}_annot.vcf'
 SnpSift extractFields -s ';' -e '.' 'result/annot/${sample}_annot.vcf' POS REF ALT \
   'ANN[*].GENE' 'ANN[*].IMPACT' 'ANN[*].EFFECT' 'GEN[0].GT' 'GEN[0].DP' 'GEN[0].AF' 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' | \
   awk -v s='${sample}' 'BEGIN {OFS="\t"} NR==1 {sub(/^#/, ""); print "Sample",\$0; next} {print s,\$0}' > result/annot/var.tsv
 """
}

process ANNOTATION_WORKBOOK {
 label 'wf_common'
 tag { sample }
 publishDir { "${params.out_dir}/${sample}" }, mode: 'copy',
   saveAs: { name -> name == 'result' ? null : name.replaceFirst(/^result\//, '') }
 input:
 tuple val(sample), path(source, stageAs: 'source_result/*')
 output:
 tuple val(sample), path('result/*', type: 'any'), emit: results
 script:
 """
 cp -aL --reflink=auto source_result result
 python3 '${projectDir}/bin/annotation_workbook.py' result/annot/var.tsv result/annot/var.xlsx
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

process TOOL_VERSIONS {
 container { image }
 tag { group }
 input:
 tuple val(group), val(image), val(tools)
 output:
 path "${group}.versions.tsv"
 script:
 """
 bash '${projectDir}/bin/tool_versions.sh' '${group}.versions.tsv' ${tools.join(' ')}
 """
}

process WORKFLOW_REPORT {
 label 'wf_common'
 publishDir params.out_dir, mode: 'copy'
 input:
 path summaries, stageAs: 'samples??/*'
 val missing_samples
 val sample_details
 path tool_versions, stageAs: 'versions/*'
 output:
 path 'wf-virus-longread-report.html'
 path 'results.json'
 script:
 def settings = [workflow_version: workflow.manifest.version, nextflow_version: nextflow.version.toString(),
   missing_samples: missing_samples,
   sample_metadata: sample_details,
   run_name: workflow.runName, session_id: workflow.sessionId.toString(),
   reference_based_assembly: params.reference_based_assembly, reference: params.reference,
   run_annotation: params.run_annotation, annotation_gff: params.annotation_gff, snpeff_cache: params.snpeff_cache,
   min_mapq: params.min_mapq, min_coverage: params.min_coverage,
   min_variant_qual: params.min_variant_qual, min_variant_depth: params.min_variant_depth, min_variant_af: params.min_variant_af,
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
 python3 '${projectDir}/bin/report.py' aggregate '${encoded}' --versions-dir versions ${summary_files.collect { "'${it}'" }.join(' ')}
 """
}

workflow {
 if (params.help) { log.info 'See README.md. Required: --fastq. Default mode: de novo contigs only.'; return }
 if (!params.fastq) error 'Select --fastq (FASTQ file or directory)'
 if (!params.out_dir) error 'Unable to resolve output directory; specify --out_dir'
 ['threads','flye_genome_size','min_mapq','min_coverage','mask_depth','min_read_length','min_variant_qual','min_variant_depth'].each { k ->
   if (!(params[k].toString() ==~ /[0-9]+/)) error "--${k} must be an integer"
 }
 if (!(params.min_variant_af.toString() ==~ /[0-9]+(\.[0-9]+)?/) || (params.min_variant_af as double) > 1) error '--min_variant_af must be a fraction from 0 to 1'
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
     def resources
     if (params.snpeff_config || params.snpeff_data || params.snpeff_db) {
       ['snpeff_config','snpeff_data','snpeff_db'].each { k -> if (!params[k]) error "Legacy annotation requires all three snpeff_config, snpeff_data, snpeff_db parameters" }
       if (params.annotation_gff) error 'Choose annotation_gff automatic mode or legacy snpEff resources, not both'
       if (!(params.snpeff_db ==~ /[A-Za-z0-9_.-]+/)) error 'Invalid snpEff database ID'
       def config = file(params.snpeff_config, checkIfExists: true)
       def data = file(params.snpeff_data, checkIfExists: true)
       if (!data.isDirectory() || !file("${data}/${params.snpeff_db}/snpEffectPredictor.bin").exists()) error 'Legacy snpEff DB binary not found'
       def metadata = [mode:'legacy_manual', db_id:params.snpeff_db,
         validation:[identity_note:'Legacy manually selected database; reference sequence identity has not been automatically verified.', independent_cds_check:false, independent_protein_check:false]]
       def encoded = groovy.json.JsonOutput.toJson(metadata).bytes.encodeBase64().toString()
       resources = Channel.value(tuple(config, data, params.snpeff_db, encoded))
     } else {
       def gff = file(params.annotation_gff ?: "${projectDir}/model_map.tsv", checkIfExists:true)
       if (params.annotation_gff && !gff.isFile()) error 'annotation_gff must be a GFF3 file'
       def cache = file(params.snpeff_cache ?: "${params.out_dir}/annotation_database")
       if (!cache.exists()) cache.mkdirs()
       if (!cache.isDirectory()) error 'snpeff_cache must be a writable directory'
       PREPARE_SNPEFF_DB(ref, gff, cache)
       // Static preparation inputs produce a value channel shared by all samples.
       resources = PREPARE_SNPEFF_DB.out.map { bundle ->
         def metadata = new groovy.json.JsonSlurper().parseText(bundle.resolve('metadata.json').text)
         tuple(bundle.resolve('snpEff.config'), bundle.resolve('data'), metadata.db_id,
           groovy.json.JsonOutput.toJson(metadata).bytes.encodeBase64().toString())
       }
     }
     SNPEFF_ANNOTATE(REFERENCE_ANALYSIS.out.results.map { sample, result ->
       tuple(sample, result.find { it.name == 'varcalling' }.resolve("${sample}_norm.vcf.gz"))
     }, resources)
     SNPSIFT_TABLE(REFERENCE_ANALYSIS.out.results.join(SNPEFF_ANNOTATE.out))
     ANNOTATION_WORKBOOK(SNPSIFT_TABLE.out.results)
     analysis_results = ANNOTATION_WORKBOOK.out.results.map { sample, result -> tuple(sample, result, file("${projectDir}/model_map.tsv"), file("${projectDir}/clair3_model_catalog.json")) }
   }
 } else {
   DE_NOVO_CONTIGS(INGEST.out.reads)
   DE_NOVO_COVERAGE(INGEST.out.reads.join(DE_NOVO_CONTIGS.out))
   analysis_results = DE_NOVO_CONTIGS.out.join(DE_NOVO_COVERAGE.out)
 }
 metadata_results = INGEST.out.reads.map { sample, reads, metadata -> tuple(sample, metadata) }
 SAMPLE_SUMMARY(metadata_results.join(QC.out).join(analysis_results).join(INGEST.out.stats)
   .map { sample, metadata, qc, analysis, flye, coverage, stats -> tuple(sample, metadata, qc, analysis, stats, flye, coverage) })
 def version_jobs = [
   ['common', params.images.wf_common, ['fastcat','bgzip','Python']],
   ['mapping', params.images.wf_bacterial, ['mini_align','minimap2','samtools','bcftools']],
   ['nanoplot', params.container ?: params.images.nanoplot, ['NanoPlot']]]
 if (params.reference_based_assembly) {
   version_jobs += [['clair3', params.container ?: params.images.clair3, ['Clair3']],
                    ['coverage', params.container ?: params.images.deeptools, ['bamCoverage']]]
   if (params.run_annotation) version_jobs += [
     ['snpeff', params.container ?: params.images.snpeff, ['snpEff','Java']],
     ['snpsift', params.container ?: params.images.snpsift, ['SnpSift']]]
 } else {
   version_jobs += [['flye', params.images.wf_bacterial, ['Flye']]]
 }
 TOOL_VERSIONS(Channel.fromList(version_jobs))
 WORKFLOW_REPORT(SAMPLE_SUMMARY.out.toList(), missing_samples, sample_details, TOOL_VERSIONS.out.toList())
}
