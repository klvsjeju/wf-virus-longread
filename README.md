# wf-virus-longread

A Nextflow DSL2 workflow for VZV reference-based analysis or de novo contig assembly from FASTQ reads.

## EPI2ME inputs

1. **FASTQ file or folder**: Select a FASTQ/FQ file (plain or gzip compressed), a folder of FASTQ files for one sample, or a parent folder containing barcode/sample subfolders. One subdirectory level is supported. Use `sample` to name a single file or a single-sample folder. An optional `sample_sheet` CSV must contain `barcode,alias` columns and match all selected sample folders.
2. **Reference FASTA**: Upload/select your own uncompressed FASTA file. Enable **Reference-based assembly** to run reference analysis. No reference is bundled.
3. **Manual basecalling model**: Leave unset to detect the model from the FASTQ header tags `basecall_model_version_id`, `basecall_model`, or `model_version_id`. All read headers are inspected, and metadata is saved to `metadata/<sample>/basecaller.json`. A manual selection takes priority over detection. Mixed detected models cause an error unless a manual override is supplied. De novo assembly also accepts FASTQ files without model tags. This setting identifies the model used for existing reads; it does not run basecalling.
4. **De novo**: Disable Reference-based assembly to run QC and Flye. Outputs include `denovo/<sample>/<sample>.contigs.fasta` and the Flye output folder. No downstream Medaka/Clair3 polishing or reference-based consensus step is performed. Flye's own assembly and internal polishing steps are included.

Reference mode performs alignment, Clair3 variant calling, consensus generation, and snpEff annotation. Enable optional `run_annotation` only when snpEff annotation is needed; only then are `snpeff_config`, `snpeff_data`, and `snpeff_db` required. Clair3 checkpoints are selected automatically from the detected or manually specified basecalling model. No checkpoint directory or mapping TSV input is shown in EPI2ME.

## CLI execution

```bash
docker build -t wf-virus-longread:0.3.0 .
nextflow run . -profile standard --fastq /path/to/reads --out_dir output
nextflow run . -profile standard -params-file params.reference.example.json
```

Update paths in the example JSON files for your environment. When out_dir is unset, results are written to analysis_result inside the selected FASTQ directory, or beside the selected FASTQ file. Running from GitHub also requires building `wf-virus-longread:0.3.0` locally or specifying an accessible image with `--container`.

Input UI formats (`path`, `file-path`) and the model selection list follow the [wf-bacterial-genomes schema](https://github.com/epi2me-labs/wf-bacterial-genomes/blob/master/nextflow_schema.json). Its Medaka stage has not been added to this workflow's de novo branch.

The EPI2ME Desktop input UI and complete analysis have not yet been tested. See VALIDATION.txt for the validation scope.

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

## EPI2ME results

Details / Reports / Files / Logs are tabs provided by EPI2ME Desktop. The workflow does not need to create folders with those names.

| Tab | Content provided by this workflow |
| --- | --- |
| Details | Run status and input settings collected by the app. Additional run information and options are also recorded in results.json and the HTML report. |
| Reports | Top-level wf-virus-longread-report.html with per-sample models, contig/consensus counts and lengths, and links to QC and result files. NanoPlot and Nextflow execution reports are generated separately. |
| Files | Contigs, consensus, mask BED, BAM/BAI, VCF, annotations, depth, QC, metadata, and execution outputs described in output_definition.json. |
| Logs | Nextflow/process logs collected by EPI2ME. Task execution records are also saved to execution/trace.tsv. |

The aggregate report waits for all per-sample analysis, metadata, and QC results. If analysis fails, the aggregate report may not be generated; check Logs for the cause. `results.json` uses this workflow's own summary format rather than the official workflows' internal results/checkpoints format. `output_definition.json` follows the output description format used by official workflows.

The report can be opened without external web resources. Preserve the output folder structure to use links to result files. If the EPI2ME embedded viewer restricts file links, open the files through the Files tab. Report discovery and tab rendering across Desktop versions have not yet been tested.


## Automatic Clair3 models (0.3.0)

Clair3 is pinned to 1.2.0, which uses TensorFlow checkpoints. Clair3 2.x PyTorch checkpoints are incompatible with this workflow. The model catalog records exact model download URLs from Rerio revision `198f90f33434437001ad0d023e3ca912cdbca47c`.

Reference runs prepare one checkpoint per distinct selected basecaller and share it among matching samples. The workflow checks for matching preinstalled models in `/opt/clair3-models` and `/opt/conda/bin/models`, then downloads the required archive from the official Oxford Nanopore model CDN if absent. Network access is required for downloads. Completed preparation tasks can be reused with Nextflow `-resume` when their work directory/cache remain available; this is not a permanent cross-run download cache.

Automatic exact matches:

- R10.4.1 E8.2 400 bps HAC/SUP: v4.0.0, v4.1.0, v4.2.0, v4.3.0, v5.0.0, v5.2.0.
- R10.4.1 E8.2 260 bps HAC/SUP: v4.0.0, v4.1.0.

Unsupported or unknown basecaller names fail with guidance rather than silently using a different model. The manual basecaller field still includes older choices usable with custom checkpoints. These are upstream chemistry/basecaller matches, not VZV-specific accuracy validation. Some upstream models are marked deprecated. See [Rerio](https://github.com/nanoporetech/rerio) for model information and license terms.

Model identity, source and file hashes are saved in `reference/<sample>/result/model_selection.json` and included in the aggregate report/results.json. Downloaded archives are hashed for provenance; unless the catalog supplies an expected sha256, the recorded hash is not an independent integrity guarantee. Archives are fetched over HTTPS and checked for unsafe paths/links and nonempty TensorFlow checkpoint files.

Advanced CLI users may still supply `--model_dir` for a compatible custom checkpoint or `--model_dir` plus `--model_map` for a custom library. These options have been removed from the EPI2ME input schema. Neither is required for automatic selection.

After updating from 0.2.0, build the new image with `docker build -t wf-virus-longread:0.3.0 .`, make it available to the Docker engine used by EPI2ME, push the source changes and update the installed workflow. GitHub import does not build the image. Older Word guides describe 0.2.0 and do not include these changes.

## EPI2ME launch form organization

The workflow input groups appear in this order: Input options, Sample options, Advanced option, Reference annotation, Miscellaneous option. Sample options contains sample_sheet followed by sample. Advanced option contains the manual basecaller choice and Flye/alignment/masking/consensus settings. Miscellaneous option contains threads and out_dir. Task memory is fixed at 4 GB; the Docker image is configured internally and is not shown as a launch field. The final Nextflow configuration area is provided by EPI2ME Desktop rather than duplicated in the workflow schema. Actual rendering and availability of this built-in area depend on the Desktop version.


## Memory, container and default output

Each task requests a fixed 4 GB of memory. The memory and Docker image fields are removed from the EPI2ME form. The container remains necessary for the standard Docker profile; its configured image must be available to the Docker engine. An advanced CLI --container override remains supported, while --memory does not change the fixed process memory.

When out_dir is omitted or empty:
- FASTQ file /data/reads/sample.fastq.gz -> /data/reads/analysis_result
- FASTQ directory /data/reads -> /data/reads/analysis_result

The input location must be writable and visible to the execution environment. You can still select an explicit Output directory. Different runs using the same input folder share this default destination and may overwrite/mix published files; specify separate output directories to retain separate runs. Nextflow creates execution outputs and publishes analysis files there during execution. EPI2ME may supply its own output directory; if so that explicit value takes priority. Final Desktop output handling has not been tested.


## Optional reference annotation

`run_annotation` defaults to false. Reference mode still generates alignment, depth, variants, optional mask and consensus, but does not run snpEff/SnpSift or create annotation outputs. Enable Run reference annotation to additionally produce annot/annot.vcf and annot/var.tsv. Supply snpeff_config (snpEff configuration file defining the database), snpeff_data (database root containing DB_ID/snpEffectPredictor.bin), and snpeff_db (matching database ID). The database must match the input reference coordinates and contig names. De novo does not run this annotation stage.


## Fixed Flye read type

De novo assembly always runs Flye with `--nano-hq`, following wf-bacterial-genomes. The Flye read type field and configurable flye_read_type parameter have been removed. Reports record the fixed nano-hq value. Additional Flye options remain available, but read-type flags cannot be supplied through flye_opts.


## Consensus allele-frequency example

To apply only variants with AF strictly greater than 50%, enter this in Additional consensus options:

```text
--include 'FORMAT/AF>0.5'
```

AF is expressed as a fraction (0.5 means 50%), so use 0.5 rather than 50. This expression requires a numeric FORMAT/AF field in the Clair3 VCF; check the VCF header if using a different caller/version. A site at exactly 50% is excluded. This filters variants applied to the consensus, not the published VCF. Excluded variants leave the reference sequence at those positions unless masked separately. To include exactly 50%, use FORMAT/AF>=0.5 instead.


## Workflow name

The workflow is named wf-virus-longread. Its manifest expects the GitHub repository https://github.com/klvsjeju/wf-virus-longread. Rename the GitHub repository or publish the project at that URL before importing under the new name. The report is wf-virus-longread-report.html and the default container is wf-virus-longread:0.3.0. This is a naming change; VZV-oriented genome-size defaults and reference resources have not been generalized automatically.
