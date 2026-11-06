# wf-vzv

A Nextflow DSL2 workflow for VZV reference-based analysis or de novo contig assembly from FASTQ reads.

## EPI2ME inputs

1. **FASTQ file or folder**: Select a FASTQ/FQ file (plain or gzip compressed), a folder of FASTQ files for one sample, or a parent folder containing barcode/sample subfolders. One subdirectory level is supported. Use `sample` to name a single file or a single-sample folder. An optional `sample_sheet` CSV must contain `barcode,alias` columns and match all selected sample folders.
2. **Reference FASTA**: Upload/select your own uncompressed FASTA file. Enable **Reference-based assembly** to run reference analysis. No reference is bundled.
3. **Manual basecalling model**: Leave unset to detect the model from the FASTQ header tags `basecall_model_version_id`, `basecall_model`, or `model_version_id`. All read headers are inspected, and metadata is saved to `metadata/<sample>/basecaller.json`. A manual selection takes priority over detection. Mixed detected models cause an error unless a manual override is supplied. De novo assembly also accepts FASTQ files without model tags. This setting identifies the model used for existing reads; it does not run basecalling.
4. **De novo**: Disable Reference-based assembly to run QC and Flye. Outputs include `denovo/<sample>/<sample>.contigs.fasta` and the Flye output folder. No downstream Medaka/Clair3 polishing or reference-based consensus step is performed. Flye's own assembly and internal polishing steps are included.

Reference mode performs alignment, Clair3 variant calling, consensus generation, and snpEff annotation. In addition to the reference, it requires `model_dir`, `snpeff_config`, `snpeff_data`, and `snpeff_db`. The basecalling model name and Clair3 checkpoints are separate settings. Supply a checkpoint folder directly, or provide a checkpoint library with a verified `model_map` TSV. Models are not downloaded automatically.

## CLI execution

```bash
docker build -t wf-vzv:0.2.0 .
nextflow run . -profile standard --fastq /path/to/reads --out_dir output
nextflow run . -profile standard -params-file params.reference.example.json
```

Update paths in the example JSON files for your environment. The default output directory is `output`. Running from GitHub also requires building `wf-vzv:0.2.0` locally or specifying an accessible image with `--container`.

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
| Reports | Top-level wf-vzv-report.html with per-sample models, contig/consensus counts and lengths, and links to QC and result files. NanoPlot and Nextflow execution reports are generated separately. |
| Files | Contigs, consensus, mask BED, BAM/BAI, VCF, annotations, depth, QC, metadata, and execution outputs described in output_definition.json. |
| Logs | Nextflow/process logs collected by EPI2ME. Task execution records are also saved to execution/trace.tsv. |

The aggregate report waits for all per-sample analysis, metadata, and QC results. If analysis fails, the aggregate report may not be generated; check Logs for the cause. `results.json` uses this workflow's own summary format rather than the official workflows' internal results/checkpoints format. `output_definition.json` follows the output description format used by official workflows.

The report can be opened without external web resources. Preserve the output folder structure to use links to result files. If the EPI2ME embedded viewer restricts file links, open the files through the Files tab. Report discovery and tab rendering across Desktop versions have not yet been tested.
