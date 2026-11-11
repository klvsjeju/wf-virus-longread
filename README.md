# wf-virus-longread

A Nextflow DSL2 workflow for VZV reference-based analysis or de novo contig assembly from FASTQ reads.

## EPI2ME inputs

1. **FASTQ file or folder**: Select a FASTQ/FQ file (plain or gzip compressed), a folder of FASTQ files for one sample, or a parent folder containing barcode/sample subfolders. One subdirectory level is supported. Leave `sample` empty for an automatic input-derived name, or provide a name for a single input / a filter for multiplex input. An optional `sample_sheet` CSV maps barcode or alias folders to sample names as described below.
2. **Reference FASTA**: Upload/select your own uncompressed FASTA file. Enable **Reference-based assembly** to run reference analysis. No reference is bundled.
3. **Manual basecalling model**: Leave unset to detect the model from fastcat statistics. If model counts are incomplete, the input helper also checks explicit model tags and versioned models in `RG:Z:` read-group tags. Metadata is saved to `<sample>/metadata/basecaller.json`. A manual selection takes priority over detection. Mixed detected models cause an error unless a manual override is supplied. De novo assembly also accepts FASTQ files without model tags. This setting identifies the model used for existing reads; it does not run basecalling.
4. **De novo**: Disable Reference-based assembly to run QC and Flye. Outputs include `<sample>/<sample>.contigs.fasta` and the Flye output folder. No downstream Medaka/Clair3 polishing or reference-based consensus step is performed. Flye's own assembly and internal polishing steps are included.

Reference mode performs alignment, Clair3 variant calling, consensus generation, and snpEff annotation. Enable optional `run_annotation` only when snpEff annotation is needed; only then are `snpeff_config`, `snpeff_data`, and `snpeff_db` required. Clair3 checkpoints are selected automatically from the detected or manually specified basecalling model. No checkpoint directory or mapping TSV input is shown in EPI2ME.

## CLI execution

```bash
nextflow run . -profile standard --fastq /path/to/reads --out_dir output
nextflow run . -profile standard -params-file params.reference.example.json
```

Update paths in the example JSON files for your environment. When out_dir is unset, results are written to analysis_result inside the selected FASTQ directory, or beside the selected FASTQ file. The standard profile pulls the public images listed below when they are absent from the Docker cache.

Input UI formats (`path`, `file-path`) and the model selection list follow the [wf-bacterial-genomes schema](https://github.com/epi2me-labs/wf-bacterial-genomes/blob/master/nextflow_schema.json). Its Medaka stage has not been added to this workflow's de novo branch.

See VALIDATION.txt for the validation scope.

## Sample names and sample sheets (0.3.3)

Naming and folder matching follow the installed wf-bacterial-genomes v2.0.2 ingress rules:

| Input | `sample` empty | `sample` specified |
| --- | --- | --- |
| One FASTQ file | File simpleName (all extensions removed; `isolate.fastq.gz` becomes `isolate`) | Use the supplied name |
| Folder containing FASTQ files directly | Folder baseName | Use the supplied name |
| Parent folder with sample subfolders, no sheet | Use each subfolder baseName | Select only the matching subfolder baseName |
| Parent folder with sample subfolders and sheet | Match folders against `barcode` or `alias`; use the sheet alias | Select the row matching the supplied barcode or alias; use its alias |

Spaces in final names become underscores. The existing filename safety check still requires names to start with a letter or digit and contain only letters, digits, underscores, periods or hyphens. Duplicate final names are rejected. Names come from input paths or the sheet, not FASTQ read headers. A sheet does not rename single-file or flat-folder input.

The CSV requires unique `barcode,alias` columns. Barcodes must start with `barcode`, followed by at least two digits, with the same length in every row. Aliases must be nonempty and must not start with `barcode`. Optional `type` values are `test_sample`, `positive_control`, `negative_control` and `no_template_control`; optional `analysis_group` must be populated for every row.

Folders absent from the sheet are ignored with a warning. Sheet entries without matching FASTQ folders are recorded in the HTML report and `results.json` as missing input; they are not analysed. If only missing entries remain selected, the workflow produces a report without analysis tasks. Supplying both a barcode folder and its alias folder for the same sheet row is an error. Without a sheet, `sample` must match the original folder baseName; with a sheet it must match the original barcode or alias, before spaces are replaced.

Reference result names use the resolved sample alias: `<sample>_norm.pass.vcf.gz` (PASS variants), `<sample>_norm.vcf.gz` (all normalized variants), `<sample>_merge_output.vcf.gz` (caller VCF), optional `<sample>_annot.vcf`, and `<sample>_consensus.fa`. VCF indexes retain the matching sample prefix. Pileup/full-alignment VCFs and their indexes, including nested per-chunk VCFs, are omitted from published results after calling finishes; original calling task outputs in the Nextflow work cache are retained. BAM, depth, de novo contig FASTA and NanoPlot output prefixes also use the resolved alias. A manually entered literal `sample` remains a literal name; leave it empty for input-derived naming or provide the actual name / sample sheet alias.

For a single reference sequence, the consensus FASTA header is exactly `><sample>`. With multiple reference sequences, the first header is `><sample>` and subsequent records use `><sample>|contig2`, `><sample>|contig3`, etc. to keep identifiers unique. Only FASTA identifiers change; sequence generation and VCF reference chromosome identifiers are preserved.

## Manual additional options

Enter additional argument strings under Advanced options in EPI2ME. Leave them empty to retain the existing behavior.

| Setting | Command | Example |
| --- | --- | --- |
| `flye_opts` | Flye in de novo mode | `--meta --iterations 2 --asm-coverage 50` |
| `mask_depth_opts` | samtools depth for mask calculation | `-q 15 -Q 20` |
| `consensus_opts` | bcftools consensus in reference mode | `--mask-with N --mark-del -` |

The mask BED contains positions with depth below `mask_depth`. For example, `mask_depth=30` with `mask_depth_opts="-q 15 -Q 20"` masks positions with depth below 30 after requiring base quality of at least 15 and mapping quality of at least 20. Setting `mask_depth=0` disables BED mask generation and application. Additional depth options also affect the saved `.depth` file. The default `-aa -J` flags remain enabled. Depth options cannot restore reads already removed by `min_mapq` during alignment. Restricting depth calculation with `-r` or `-b` leaves positions outside that range untreated by this BED mask.

```json
{
  "mask_depth": 30,
  "mask_depth_opts": "-q 15 -Q 20",
  "consensus_opts": "--mask-with N --mark-del -",
  "flye_opts": "--meta --iterations 2"
}
```

On the CLI, quote the entire option string:

```bash
nextflow run . --fastq /path/to/reads --flye_opts '--meta --iterations 2'
```

Workflow-managed options, including input/output paths, reference, sample, mask paths, threads, and genome size, cannot be duplicated in additional argument strings. Use the dedicated settings instead. Quoted values are passed as individual arguments. Shell commands, pipes, redirects, and variable expansion are not executed. File paths in additional options are not automatically staged into Nextflow task directories or containers. Unsupported tool options cause the relevant tool to fail.

Option documentation: [samtools depth](https://www.htslib.org/doc/samtools-depth.html), [bcftools consensus](https://samtools.github.io/bcftools/bcftools.html#consensus), [Flye](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md).

## Result folder layout

All per-sample outputs are grouped directly under `<out_dir>/<sample>/`:

```text
<out_dir>/
  <sample>/
    alignment/             # reference mode
    clair3/                # reference mode: VCFs, consensus, mask, statistics
    annot/                 # optional snpEff annotation
    model_selection.json   # reference mode
    qc/NanoPlot/
    metadata/fastcat_stats/
    metadata/basecaller.json
    <sample>.contigs.fasta  # de novo mode
    flye/                  # de novo mode
  wf-virus-longread-report.html
  results.json
  execution/
```

Each additional sample has its own folder. The names `execution`, `results.json` and `wf-virus-longread-report.html` are reserved for shared outputs (case-insensitive). This replaces the older `reference/<sample>/result/`, `qc/<sample>/`, `metadata/<sample>/` and `denovo/<sample>/` nesting. Existing analysis output directories are not migrated automatically.

`output_definition.json` registers the output root, each sample directory, and its analysis/QC/metadata/tool subfolders using the reference workflow's `inode/directory` MIME type. This lets EPI2ME treat them as folders for its native Open in Explorer action. The action is provided by Desktop, not by HTML report links; absent optional folders are not produced solely for display.

## EPI2ME results

Details / Reports / Files / Logs are tabs provided by EPI2ME Desktop. The workflow does not need to create folders with those names.

| Tab | Content provided by this workflow |
| --- | --- |
| Details | Run status and input settings collected by the app. Additional run information and options are also recorded in results.json and the HTML report. |
| Reports | Top-level wf-virus-longread-report.html using the reference workflow's EPI2ME Labs design, sample dropdowns, interactive plots and searchable tables. Includes sample/read summaries, assembly or variant/coverage/consensus statistics, optional snpEff results, models and file links. NanoPlot and Nextflow execution reports are generated separately. |
| Files | Contigs, consensus, mask BED, BAM/BAI, VCF, annotations, depth, QC, metadata, and execution outputs described in output_definition.json. |
| Logs | Nextflow/process logs collected by EPI2ME. Task execution records are also saved to execution/trace.tsv. |

The aggregate report waits for all per-sample analysis, metadata, and QC results. If analysis fails, the aggregate report may not be generated; check Logs for the cause. `results.json` uses this workflow's own summary format rather than the official workflows' internal results/checkpoints format. `output_definition.json` follows the output description format used by official workflows.

The report can be opened without external web resources. Preserve the output folder structure to use links to result files. If the EPI2ME embedded viewer restricts file links, open the files through the Files tab. Report discovery and tab rendering across Desktop versions have not yet been tested.

The report uses the same `ezcharts.LabsReport`, `Tabs`, `DataTable`, fastcat QC and `EZChart` components as installed wf-bacterial-genomes v2.0.2. It runs inside the existing public wf-common image; no custom image build or new launch setting is required. The workflow version remains 0.3.3.

- Sample summary retains barcode/alias metadata and lists missing input samples.
- Read quality, length and base yield plots reuse fastcat histograms without rereading FASTQ. Median read length, yield and N50 use histogram bins. The histogram yield plot handles samples where the longest bin contributes at least half of the yield, and the quality axis includes Q scores above 30.
- De novo shows Flye contig length, coverage, circular/repeat flags and multiplicity. Reference mode shows full normalized VCF counts before PASS filtering, Ts/Tv, complementary substitution heatmaps, indel length counts, coverage tables/plots and consensus/mask links.
- Coverage statistics use all assessed positions from the existing samtools depth output, including zero-depth positions. Region-restricting depth options change which positions are assessed. Plots retain at most 500 mean-depth bins per reference, displaying the 20 longest references; full depth files and all-reference statistics remain available. Strand-specific coverage is not computed.
- Optional snpEff tables show impact/effect counts and at most 500 annotation records per sample, with a link to the full TSV. Counts include each distinct impact/effect once per annotated record. Annotation and report variant statistics both use the full normalized VCF, including non-PASS records. Consensus uses the PASS VCF. Consensus options can select a subset of PASS variants.
- Software versions lists the actual report runtime packages. The Execution section lists configured analysis image tags and run identifiers. Results retain the workflow's own JSON schema; they do not implement the reference workflow's internal checkpoint schema.


## Automatic Clair3 models (0.3.0)

Clair3 is pinned to 1.2.0, which uses TensorFlow checkpoints. Clair3 2.x PyTorch checkpoints are incompatible with this workflow. The model catalog records exact model download URLs from Rerio revision `198f90f33434437001ad0d023e3ca912cdbca47c`.

Reference runs prepare one checkpoint per distinct selected basecaller and share it among matching samples. The workflow checks for matching preinstalled models in `/opt/clair3-models`, `/opt/models` and `/opt/conda/bin/models`, then downloads the required archive from the official Oxford Nanopore model CDN if absent. Network access is required for downloads. Completed preparation tasks can be reused with Nextflow `-resume` when their work directory/cache remain available; this is not a permanent cross-run download cache.

Automatic exact matches:

- R10.4.1 E8.2 400 bps HAC/SUP: v4.0.0, v4.1.0, v4.2.0, v4.3.0, v5.0.0, v5.2.0.
- R10.4.1 E8.2 260 bps HAC/SUP: v4.0.0, v4.1.0.

Unsupported or unknown basecaller names fail with guidance rather than silently using a different model. The manual basecaller field still includes older choices usable with custom checkpoints. These are upstream chemistry/basecaller matches, not VZV-specific accuracy validation. Some upstream models are marked deprecated. See [Rerio](https://github.com/nanoporetech/rerio) for model information and license terms.

Model identity, source and file hashes are saved in `<sample>/model_selection.json` and included in the aggregate report/results.json. Downloaded archives are hashed for provenance; unless the catalog supplies an expected sha256, the recorded hash is not an independent integrity guarantee. Archives are fetched over HTTPS and checked for unsafe paths/links and nonempty TensorFlow checkpoint files.

Advanced CLI users may still supply `--model_dir` for a compatible custom checkpoint or `--model_dir` plus `--model_map` for a custom library. These options have been removed from the EPI2ME input schema. Neither is required for automatic selection.

The standard Docker profile uses the public Clair3 1.2.0 image. Older Word guides describe 0.2.0.

## EPI2ME launch form organization

The workflow input groups appear in this order: Input options, Sample options, Advanced option, Reference annotation, Miscellaneous option. Sample options contains sample_sheet followed by sample. Advanced option contains the manual basecaller choice and Flye/alignment/masking/consensus settings. Miscellaneous option contains threads and out_dir. Task memory is fixed at 4 GB; the Docker image is configured internally and is not shown as a launch field. The final Nextflow configuration area is provided by EPI2ME Desktop rather than duplicated in the workflow schema. Actual rendering and availability of this built-in area depend on the Desktop version.


## Memory, container and default output

Each task requests a fixed 4 GB of memory. The memory and Docker image fields are removed from the EPI2ME form. The standard profile uses a public image per tool group. The CLI `--container` override applies to NanoPlot, Clair3, deepTools, snpEff and SnpSift. Common tasks use the two reference workflow images.

When out_dir is omitted or empty:
- FASTQ file /data/reads/sample.fastq.gz -> /data/reads/analysis_result
- FASTQ directory /data/reads -> /data/reads/analysis_result

The input location must be writable and visible to the execution environment. You can still select an explicit Output directory. Different runs using the same input folder share this default destination and may overwrite/mix published files; specify separate output directories to retain separate runs. Nextflow creates execution outputs and publishes analysis files there during execution. EPI2ME may supply its own output directory; if so that explicit value takes priority. Final Desktop output handling has not been tested.


## Optional reference annotation

`run_annotation` defaults to false. Reference mode still generates alignment, depth, variants, optional mask and consensus, but does not run snpEff/SnpSift or create annotation outputs. Enable Run reference annotation to additionally produce annot/<sample>_annot.vcf and annot/var.tsv. Supply snpeff_config (snpEff configuration file defining the database), snpeff_data (database root containing DB_ID/snpEffectPredictor.bin), and snpeff_db (matching database ID). The database must match the input reference coordinates and contig names. De novo does not run this annotation stage.


## Fixed Flye read type

De novo assembly always runs Flye with `--nano-hq`, following wf-bacterial-genomes. The Flye read type field and configurable flye_read_type parameter have been removed. Reports record the fixed nano-hq value. Additional Flye options remain available, but read-type flags cannot be supplied through flye_opts.


## Consensus allele-frequency example

To apply only variants with AF strictly greater than 50%, enter this in Additional consensus options:

```text
--include 'FORMAT/AF>0.5'
```

AF is expressed as a fraction (0.5 means 50%), so use 0.5 rather than 50. This expression requires a numeric FORMAT/AF field in the Clair3 VCF; check the VCF header if using a different caller/version. A site at exactly 50% is excluded. This filters variants applied to the consensus, not the published VCF. Excluded variants leave the reference sequence at those positions unless masked separately. To include exactly 50%, use FORMAT/AF>=0.5 instead.


## Workflow name

The workflow is named wf-virus-longread. Its manifest expects the GitHub repository https://github.com/klvsjeju/wf-virus-longread. Rename the GitHub repository or publish the project at that URL before importing under the new name. The report is wf-virus-longread-report.html. This is a naming change; VZV-oriented genome-size defaults and reference resources have not been generalized automatically.


## Nextflow 23.04.2 compatibility fix (0.3.1)

Task tags and sample-specific publish directories use dynamic closures so task inputs are resolved at execution time. The newer input arity option was removed, and single staged paths are explicitly wrapped as lists before building read/report arguments. Directory discovery uses the Nextflow Path directly rather than calling File.toPath(). The minimum Nextflow version is now 23.04.2. Workflow version 0.3.1 continues to use the wf-virus-longread:0.3.0 analysis image; no image dependency change is required for this script fix.

Verified with Nextflow 23.04.2 and a temporary Java 17 runtime: two-sample de novo, two-sample reference with a shared custom checkpoint and annotation disabled, and single-file de novo with the default analysis_result output. External bioinformatics commands were replaced by test doubles. These runs validate process/channel/file wiring and reporting, not real Flye/Clair3 accuracy or container execution. The reference model preparation helper, FASTQ merge and report code executed normally against synthetic inputs. Python unit tests also passed (13 tests).


## Container dependency conflict fix

DeepTools 3.5.6 requires NumPy 2, while Clair3 1.2.0 uses TensorFlow requiring NumPy below 2. Use deepTools 3.5.5, numpy>=1.23.5,<2 and CPU TensorFlow 2.15.0. The complete linux-64 dependency set passed a Conda libmamba dry-run; see dependency-solve-summary.json. These pins apply when building the optional combined analysis image with Dockerfile.


## FASTQ read-group model detection fix

The input helper accepts explicit basecall_model_version_id, basecall_model and model_version_id key=value tags, and now also extracts versioned DNA/RNA basecaller names embedded in RG:Z:runid_model_barcode read-group tags (or RG=...). Explicit valid model tags take priority over the read-group fallback. Unknown placeholder values do not count as a model. Barcode and modified-base suffixes are excluded from the versioned basecaller name. Model detection does not guess from read sequences, quality or filenames. The Python header parser runs as a fallback after fastcat when model metadata is incomplete.

basecaller.json now records records_scanned, records_with_model and model_detection_sources. The RG fix was verified with synthetic plain/gzip FASTQ, explicit-tag precedence, mixed-model detection and missing metadata. It has not yet been checked against the user's failing FASTQ. The current release also uses fastcat and parallel bgzip for ingestion.


## Shared steps and public images (0.3.2)

The shared FASTQ ingestion, reference alignment and de novo assembly commands follow the installed wf-bacterial-genomes v2.0.2 release. Ingestion uses fastcat and threaded bgzip, writes file-level statistics under <sample>/metadata/fastcat_stats, and selects the basecaller from those statistics with the header fallback described above. Alignment uses mini_align with the map-ont preset and retains the workflow MAPQ filter. Flye uses the reference workflow's nano-hq input mode. The bacterial workflow's read length, DUST, polishing and annotation settings are specific to that workflow.

The standard Docker profile selects these public images by task label:

| Step | Image |
| --- | --- |
| FASTQ ingestion and reporting | ontresearch/wf-common:shafdd79f8e4a6faad77513c36f623693977b92b08e |
| Index, alignment, Flye and consensus postprocessing | ontresearch/wf-bacterial-genomes:shab005a1bd7d4bb439bd864fecf5223e2741a53631 |
| Clair3 model preparation and calling | hkubal/clair3:v1.2.0 |
| Read QC | quay.io/biocontainers/nanoplot:1.43.0--pyhdfd78af_1 |
| BAM coverage | quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 |
| Optional annotation | quay.io/biocontainers/snpeff:5.4.0c--hdfd78af_0 |
| Optional annotation table | quay.io/biocontainers/snpsift:5.2--hdfd78af_0 |

Docker pulls each image when its first task runs. To update an EPI2ME installation, install the complete updated workflow tree so main.nf, nextflow.config, bin and the model catalog come from the same revision. Use the standard profile. The local profile requires the listed tools, including fastcat, bgzip and mini_align, on PATH.


The public snpEff 5.4.0c image supports the existing Dumas DB, whose stored database version is 5.4. A temporary database built with snpEff 5.2 was also successfully annotated with this image. The original pOka DB binary was unavailable during this check.

The 0.3.2 validation used real public containers on the EPI2ME WSL Docker daemon. A single synthetic sample completed reference alignment, variant calling, masking, consensus, annotation and reporting. The same synthetic dataset completed de novo Flye assembly and reporting. See VALIDATION.txt for test scope and the two-sample channel test.
