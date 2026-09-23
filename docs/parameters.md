# All-In-One Pipeline Parameters

This page documents the user-facing parameters defined by the All-In-One pipeline. Parameters are grouped by workflow area to make the reference easier to navigate.

Command-line pipeline parameters use two hyphens, for example `--run_qc false`. Nextflow engine options such as `-profile`, `-c`, and `-params-file` use one hyphen.

> **Back to the main page:** [All-In-One README](../README.md)

## Parameter groups

- [Input and output](#input-and-output)
- [Sequencing mode](#sequencing-mode)
- [Quality control](#quality-control)
- [Host and rRNA removal](#host-and-rrna-removal)
- [Assembly](#assembly)
- [Nanopore polishing and CLC](#nanopore-polishing-and-clc)
- [Assembly validation and AutoCycler](#assembly-validation-and-autocycler)
- [Sequence search and viral analysis](#sequence-search-and-viral-analysis)
- [Read taxonomic classification](#read-taxonomic-classification)
- [Virulence-factor classifier](#virulence-factor-classifier)
- [Contig characterization](#contig-characterization)
- [Automated read mapping](#automated-read-mapping)
- [Krona rendering](#krona-rendering)
- [Advanced and developer options](#advanced-and-developer-options)

## Input and output

Input sample sheet, project naming, and output locations.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--samplesheet`](#samplesheet) | string | `samplesheet.csv` | CSV sample sheet with the columns `sample_id`, `fastq_1`, `fastq_2`, and `long_read`. |
| [`--project_id`](#project_id) | string | `YOU_FORGOT_TO_SET_THIS` | Project identifier used in output paths and report names. |
| [`--outdir`](#outdir) | string | `${launchDir}/results` | Directory where pipeline results and execution reports are written. |
| [`--work_dir`](#work_dir) | string | `${launchDir}/work` | Nextflow work directory. |
| [`--cleanup`](#cleanup) | boolean / null | `null` | Whether Nextflow should remove intermediate work files after successful completion. |

<a id="samplesheet"></a>
### `--samplesheet`

CSV sample sheet with the columns `sample_id`, `fastq_1`, `fastq_2`, and `long_read`.

- Type: `string`
- Default: `samplesheet.csv`

<a id="project_id"></a>
### `--project_id`

Project identifier used in output paths and report names.

- Type: `string`
- Default: `YOU_FORGOT_TO_SET_THIS`

<a id="outdir"></a>
### `--outdir`

Directory where pipeline results and execution reports are written.

- Type: `string`
- Default: `${launchDir}/results`

<a id="work_dir"></a>
### `--work_dir`

Nextflow work directory.

- Type: `string`
- Default: `${launchDir}/work`

<a id="cleanup"></a>
### `--cleanup`

Whether Nextflow should remove intermediate work files after successful completion.

- Type: `boolean / null`
- Default: `null`

## Sequencing mode

Select the sequencing-data mode represented by the sample sheet.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--shortreads`](#shortreads) | boolean | `false` | Run in paired short-read mode. Set only one sequencing-mode switch to true. |
| [`--longreads`](#longreads) | boolean | `false` | Run in long-read mode. Set only one sequencing-mode switch to true. |
| [`--hybrid`](#hybrid) | boolean | `false` | Run with paired short reads and long reads. Set only one sequencing-mode switch to true. |

<a id="shortreads"></a>
### `--shortreads`

Run in paired short-read mode. Set only one sequencing-mode switch to true.

- Type: `boolean`
- Default: `false`

<a id="longreads"></a>
### `--longreads`

Run in long-read mode. Set only one sequencing-mode switch to true.

- Type: `boolean`
- Default: `false`

<a id="hybrid"></a>
### `--hybrid`

Run with paired short reads and long reads. Set only one sequencing-mode switch to true.

- Type: `boolean`
- Default: `false`

## Quality control

Read trimming, FastQC, NanoPlot, MultiQC, and reporting options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--run_qc`](#run_qc) | boolean | `true` | Run the primary read quality-control workflow. |
| [`--run_qc_stats`](#run_qc_stats) | boolean | `true` | Generate FastQC, NanoPlot, MultiQC, and summary statistics where applicable. |
| [`--trim`](#trim) | boolean | `true` | Trim and filter reads before downstream analysis. |
| [`--fastp_opts`](#fastp_opts) | string | `-q 20 -e 20 --cut_front --cut_tail -w 16 -l 50` | Additional command-line options passed to fastp for paired short reads. |
| [`--fastp_long_opts`](#fastp_long_opts) | string | `` | Additional command-line options passed to fastp-long. |
| [`--quality_phread_fastp_long`](#quality_phread_fastp_long) | integer | `8` | Minimum qualified Phred score used by fastp-long. |
| [`--exercise_report`](#exercise_report) | boolean | `false` | Create the exercise-style final report. |
| [`--agnostic_read_analysis`](#agnostic_read_analysis) | boolean | `true` | Enable analyses that do not require an expected organism assignment. |

<a id="run_qc"></a>
### `--run_qc`

Run the primary read quality-control workflow.

- Type: `boolean`
- Default: `true`

<a id="run_qc_stats"></a>
### `--run_qc_stats`

Generate FastQC, NanoPlot, MultiQC, and summary statistics where applicable.

- Type: `boolean`
- Default: `true`

<a id="trim"></a>
### `--trim`

Trim and filter reads before downstream analysis.

- Type: `boolean`
- Default: `true`

<a id="fastp_opts"></a>
### `--fastp_opts`

Additional command-line options passed to fastp for paired short reads.

- Type: `string`
- Default: `-q 20 -e 20 --cut_front --cut_tail -w 16 -l 50`

<a id="fastp_long_opts"></a>
### `--fastp_long_opts`

Additional command-line options passed to fastp-long.

- Type: `string`
- Default: ``

<a id="quality_phread_fastp_long"></a>
### `--quality_phread_fastp_long`

Minimum qualified Phred score used by fastp-long.

- Type: `integer`
- Default: `8`

<a id="exercise_report"></a>
### `--exercise_report`

Create the exercise-style final report.

- Type: `boolean`
- Default: `false`

<a id="agnostic_read_analysis"></a>
### `--agnostic_read_analysis`

Enable analyses that do not require an expected organism assignment.

- Type: `boolean`
- Default: `true`

## Host and rRNA removal

Host-depletion and ribosomal-RNA filtering options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--map2reads`](#map2reads) | boolean | `false` | Map reads against the host reference and retain host-removed reads. |
| [`--host_db`](#host_db) | string | `` | Host reference FASTA or preconfigured host reference used for read removal. |
| [`--remove_rRNA_reads`](#remove_rRNA_reads) | boolean | `false` | Remove ribosomal-RNA reads before downstream analysis. |
| [`--bbmap_args`](#bbmap_args) | string | `k=13 usemodulo=f rebuild=f interleaved=auto fastareadlen=500 unpigz=f touppercase=t ` | Additional BBMap settings for rRNA filtering. |
| [`--rRNA_method`](#rRNA_method) | string | `bbmap_minimap2` | rRNA-removal implementation: `bbmap_minimap2` or `ribodetector`. Allowed values: `bbmap_minimap2`, `ribodetector`. |
| [`--ribodetector_mode`](#ribodetector_mode) | string | `cpu` | RiboDetector execution mode: `cpu` or `gpu`. Allowed values: `cpu`, `gpu`. |
| [`--ribodetector_short_len`](#ribodetector_short_len) | integer | `100` | Read-length setting supplied to RiboDetector for short reads. |
| [`--ribodetector_long_len`](#ribodetector_long_len) | integer | `100` | Read-length setting supplied to RiboDetector for long reads. |
| [`--ribodetector_ensure`](#ribodetector_ensure) | string | `rrna` | RiboDetector class that must be retained or reported. |
| [`--ribodetector_chunk_size`](#ribodetector_chunk_size) | integer | `256` | Number of sequences processed per RiboDetector chunk. |
| [`--ribodetector_cpu_cluster_options`](#ribodetector_cpu_cluster_options) | string | `` | Additional scheduler options for CPU RiboDetector jobs. |
| [`--ribodetector_gpu_cluster_options`](#ribodetector_gpu_cluster_options) | string | `--gres=gpu:1` | Additional scheduler options for GPU RiboDetector jobs. |
| [`--ribodetector_gpu_memory`](#ribodetector_gpu_memory) | integer | `4` | GPU-memory setting supplied to RiboDetector. |

<a id="map2reads"></a>
### `--map2reads`

Map reads against the host reference and retain host-removed reads.

- Type: `boolean`
- Default: `false`

<a id="host_db"></a>
### `--host_db`

Host reference FASTA or preconfigured host reference used for read removal.

- Type: `string`
- Default: ``

<a id="remove_rRNA_reads"></a>
### `--remove_rRNA_reads`

Remove ribosomal-RNA reads before downstream analysis.

- Type: `boolean`
- Default: `false`

<a id="bbmap_args"></a>
### `--bbmap_args`

Additional BBMap settings for rRNA filtering.

- Type: `string`
- Default: `k=13 usemodulo=f rebuild=f interleaved=auto fastareadlen=500 unpigz=f touppercase=t `

<a id="rRNA_method"></a>
### `--rRNA_method`

rRNA-removal implementation: `bbmap_minimap2` or `ribodetector`.

- Type: `string`
- Default: `bbmap_minimap2`
- Allowed values: `bbmap_minimap2`, `ribodetector`

<a id="ribodetector_mode"></a>
### `--ribodetector_mode`

RiboDetector execution mode: `cpu` or `gpu`.

- Type: `string`
- Default: `cpu`
- Allowed values: `cpu`, `gpu`

<a id="ribodetector_short_len"></a>
### `--ribodetector_short_len`

Read-length setting supplied to RiboDetector for short reads.

- Type: `integer`
- Default: `100`

<a id="ribodetector_long_len"></a>
### `--ribodetector_long_len`

Read-length setting supplied to RiboDetector for long reads.

- Type: `integer`
- Default: `100`

<a id="ribodetector_ensure"></a>
### `--ribodetector_ensure`

RiboDetector class that must be retained or reported.

- Type: `string`
- Default: `rrna`

<a id="ribodetector_chunk_size"></a>
### `--ribodetector_chunk_size`

Number of sequences processed per RiboDetector chunk.

- Type: `integer`
- Default: `256`

<a id="ribodetector_cpu_cluster_options"></a>
### `--ribodetector_cpu_cluster_options`

Additional scheduler options for CPU RiboDetector jobs.

- Type: `string`
- Default: ``

<a id="ribodetector_gpu_cluster_options"></a>
### `--ribodetector_gpu_cluster_options`

Additional scheduler options for GPU RiboDetector jobs.

- Type: `string`
- Default: `--gres=gpu:1`

<a id="ribodetector_gpu_memory"></a>
### `--ribodetector_gpu_memory`

GPU-memory setting supplied to RiboDetector.

- Type: `integer`
- Default: `4`

## Assembly

Short-read, long-read, hybrid, and metagenomic assembly options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--run_assembly`](#run_assembly) | boolean | `false` | Run one or more enabled assembly workflows. |
| [`--subassembly`](#subassembly) | boolean | `false` | Run the subassembly workflow. |
| [`--spades`](#spades) | boolean | `false` | Run SPAdes assembly. |
| [`--metaspades`](#metaspades) | boolean | `false` | Run metaSPAdes assembly. |
| [`--plasmidspades`](#plasmidspades) | boolean | `false` | Run plasmidSPAdes assembly. |
| [`--unicycler`](#unicycler) | boolean | `false` | Run Unicycler assembly. |
| [`--unicycler_mode`](#unicycler_mode) | string | `normal` | Unicycler bridging mode, such as `conservative`, `normal`, or `bold`. Allowed values: `conservative`, `normal`, `bold`. |
| [`--dragonflye`](#dragonflye) | boolean | `false` | Run Dragonflye long-read assembly. |
| [`--dragonflye_isolate`](#dragonflye_isolate) | boolean | `false` | Use isolate-oriented Dragonflye behavior. |
| [`--raven`](#raven) | boolean | `false` | Run Raven long-read assembly. |
| [`--myloasm`](#myloasm) | boolean | `false` | Run Myloasm long-read metagenome assembly. |
| [`--myloasm_hifi`](#myloasm_hifi) | boolean | `false` | Treat Myloasm input as HiFi reads. |
| [`--metaspades_default_cpus`](#metaspades_default_cpus) | integer | `254` | Default CPU request for metaSPAdes. |
| [`--metaspades_default_mem`](#metaspades_default_mem) | string | `500 GB` | Default memory request for metaSPAdes. |
| [`--metaSPAdes`](#metaSPAdes) | string | `spades.py --only-assembler --tmp-dir ${params.tmp_dir} --meta` | metaSPAdes command template. |
| [`--myloasm_cpus`](#myloasm_cpus) | integer | `50` | CPU request for Myloasm. |
| [`--myloasm_memory`](#myloasm_memory) | string | `256 GB` | Memory request for Myloasm. |
| [`--myloasm_time`](#myloasm_time) | string | `72h` | Wall-time request for Myloasm. |
| [`--myloasm_opts`](#myloasm_opts) | string | `--clean-dir` | Additional Myloasm command-line options. |
| [`--dragonflye_opts`](#dragonflye_opts) | string | `--trim --trimopts "--discard_middle" --depth 100 --opts "--iterations 5 --meta" --racon 3 --polypolish 5 --polypolish_careful --nanohq` | Additional Dragonflye command-line options. |
| [`--raven_opts`](#raven_opts) | string | `--trim --trimopts "--discard_middle" --depth 0 --racon 3` | Additional Raven command-line options. |
| [`--dragonflye_min_quality`](#dragonflye_min_quality) | integer | `8` | Minimum long-read quality used by Dragonflye-associated filtering. |
| [`--use_gsize`](#use_gsize) | boolean | `false` | Provide an explicit genome-size estimate to long-read assembly. |
| [`--gsize`](#gsize) | string | `` | Expected genome size, for example `5m`. |
| [`--threads`](#threads) | integer | `8` | Generic thread count used by legacy assembly steps. |
| [`--memory`](#memory) | integer | `8` | Generic memory value used by legacy assembly steps. |
| [`--medaka`](#medaka) | boolean | `false` | Run Medaka consensus polishing where supported. |
| [`--medaka_model`](#medaka_model) | string | `` | Medaka model name. |

<a id="run_assembly"></a>
### `--run_assembly`

Run one or more enabled assembly workflows.

- Type: `boolean`
- Default: `false`

<a id="subassembly"></a>
### `--subassembly`

Run the subassembly workflow.

- Type: `boolean`
- Default: `false`

<a id="spades"></a>
### `--spades`

Run SPAdes assembly.

- Type: `boolean`
- Default: `false`

<a id="metaspades"></a>
### `--metaspades`

Run metaSPAdes assembly.

- Type: `boolean`
- Default: `false`

<a id="plasmidspades"></a>
### `--plasmidspades`

Run plasmidSPAdes assembly.

- Type: `boolean`
- Default: `false`

<a id="unicycler"></a>
### `--unicycler`

Run Unicycler assembly.

- Type: `boolean`
- Default: `false`

<a id="unicycler_mode"></a>
### `--unicycler_mode`

Unicycler bridging mode, such as `conservative`, `normal`, or `bold`.

- Type: `string`
- Default: `normal`
- Allowed values: `conservative`, `normal`, `bold`

<a id="dragonflye"></a>
### `--dragonflye`

Run Dragonflye long-read assembly.

- Type: `boolean`
- Default: `false`

<a id="dragonflye_isolate"></a>
### `--dragonflye_isolate`

Use isolate-oriented Dragonflye behavior.

- Type: `boolean`
- Default: `false`

<a id="raven"></a>
### `--raven`

Run Raven long-read assembly.

- Type: `boolean`
- Default: `false`

<a id="myloasm"></a>
### `--myloasm`

Run Myloasm long-read metagenome assembly.

- Type: `boolean`
- Default: `false`

<a id="myloasm_hifi"></a>
### `--myloasm_hifi`

Treat Myloasm input as HiFi reads.

- Type: `boolean`
- Default: `false`

<a id="metaspades_default_cpus"></a>
### `--metaspades_default_cpus`

Default CPU request for metaSPAdes.

- Type: `integer`
- Default: `254`

<a id="metaspades_default_mem"></a>
### `--metaspades_default_mem`

Default memory request for metaSPAdes.

- Type: `string`
- Default: `500 GB`

<a id="metaSPAdes"></a>
### `--metaSPAdes`

metaSPAdes command template.

- Type: `string`
- Default: `spades.py --only-assembler --tmp-dir ${params.tmp_dir} --meta`

<a id="myloasm_cpus"></a>
### `--myloasm_cpus`

CPU request for Myloasm.

- Type: `integer`
- Default: `50`

<a id="myloasm_memory"></a>
### `--myloasm_memory`

Memory request for Myloasm.

- Type: `string`
- Default: `256 GB`

<a id="myloasm_time"></a>
### `--myloasm_time`

Wall-time request for Myloasm.

- Type: `string`
- Default: `72h`

<a id="myloasm_opts"></a>
### `--myloasm_opts`

Additional Myloasm command-line options.

- Type: `string`
- Default: `--clean-dir`

<a id="dragonflye_opts"></a>
### `--dragonflye_opts`

Additional Dragonflye command-line options.

- Type: `string`
- Default: `--trim --trimopts "--discard_middle" --depth 100 --opts "--iterations 5 --meta" --racon 3 --polypolish 5 --polypolish_careful --nanohq`

<a id="raven_opts"></a>
### `--raven_opts`

Additional Raven command-line options.

- Type: `string`
- Default: `--trim --trimopts "--discard_middle" --depth 0 --racon 3`

<a id="dragonflye_min_quality"></a>
### `--dragonflye_min_quality`

Minimum long-read quality used by Dragonflye-associated filtering.

- Type: `integer`
- Default: `8`

<a id="use_gsize"></a>
### `--use_gsize`

Provide an explicit genome-size estimate to long-read assembly.

- Type: `boolean`
- Default: `false`

<a id="gsize"></a>
### `--gsize`

Expected genome size, for example `5m`.

- Type: `string`
- Default: ``

<a id="threads"></a>
### `--threads`

Generic thread count used by legacy assembly steps.

- Type: `integer`
- Default: `8`

<a id="memory"></a>
### `--memory`

Generic memory value used by legacy assembly steps.

- Type: `integer`
- Default: `8`

<a id="medaka"></a>
### `--medaka`

Run Medaka consensus polishing where supported.

- Type: `boolean`
- Default: `false`

<a id="medaka_model"></a>
### `--medaka_model`

Medaka model name.

- Type: `string`
- Default: ``

## Nanopore polishing and CLC

Dorado, variant-calling, Nanopore, and CLC Genomics Server options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--bam_file`](#bam_file) | string | `` | Input BAM file or BAM location for Dorado workflows. |
| [`--dorado_polish`](#dorado_polish) | boolean | `false` | Run Dorado polishing. |
| [`--run_variant`](#run_variant) | boolean | `false` | Run the variant-calling branch. |
| [`--nanopore_assembly`](#nanopore_assembly) | boolean | `false` | Enable the specialized Nanopore assembly branch. |
| [`--clc`](#clc) | boolean | `false` | Run assembly using CLC Genomics Server. |
| [`--clc_grid`](#clc_grid) | string | `slurm` | CLC execution backend or grid identifier. |
| [`--clc_min_length`](#clc_min_length) | integer | `1000` | Minimum contig length reported by CLC assembly. |
| [`--clc_create_assembly_graph`](#clc_create_assembly_graph) | boolean | `true` | Request a CLC assembly graph. |
| [`--clc_create_report`](#clc_create_report) | boolean | `true` | Request a CLC assembly report. |
| [`--clc_keep_circular_contigs_under_len_threshold`](#clc_keep_circular_contigs_under_len_threshold) | boolean | `false` | Retain circular CLC contigs shorter than the minimum-length threshold. |
| [`--clc_report_export_format`](#clc_report_export_format) | string | `export_pdf` | CLC assembly-report export format. |
| [`--clc_graph_export_format`](#clc_graph_export_format) | string | `` | CLC assembly-graph export format. |

<a id="bam_file"></a>
### `--bam_file`

Input BAM file or BAM location for Dorado workflows.

- Type: `string`
- Default: ``

<a id="dorado_polish"></a>
### `--dorado_polish`

Run Dorado polishing.

- Type: `boolean`
- Default: `false`

<a id="run_variant"></a>
### `--run_variant`

Run the variant-calling branch.

- Type: `boolean`
- Default: `false`

<a id="nanopore_assembly"></a>
### `--nanopore_assembly`

Enable the specialized Nanopore assembly branch.

- Type: `boolean`
- Default: `false`

<a id="clc"></a>
### `--clc`

Run assembly using CLC Genomics Server.

- Type: `boolean`
- Default: `false`

<a id="clc_grid"></a>
### `--clc_grid`

CLC execution backend or grid identifier.

- Type: `string`
- Default: `slurm`

<a id="clc_min_length"></a>
### `--clc_min_length`

Minimum contig length reported by CLC assembly.

- Type: `integer`
- Default: `1000`

<a id="clc_create_assembly_graph"></a>
### `--clc_create_assembly_graph`

Request a CLC assembly graph.

- Type: `boolean`
- Default: `true`

<a id="clc_create_report"></a>
### `--clc_create_report`

Request a CLC assembly report.

- Type: `boolean`
- Default: `true`

<a id="clc_keep_circular_contigs_under_len_threshold"></a>
### `--clc_keep_circular_contigs_under_len_threshold`

Retain circular CLC contigs shorter than the minimum-length threshold.

- Type: `boolean`
- Default: `false`

<a id="clc_report_export_format"></a>
### `--clc_report_export_format`

CLC assembly-report export format.

- Type: `string`
- Default: `export_pdf`

<a id="clc_graph_export_format"></a>
### `--clc_graph_export_format`

CLC assembly-graph export format.

- Type: `string`
- Default: ``

## Assembly validation and AutoCycler

CheckM, CheckV, consensus assembly, and long-read subsampling options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--run_checkm`](#run_checkm) | boolean | `true` | Run CheckM, CheckM2, and CheckV assembly validation where configured. |
| [`--checkm_time`](#checkm_time) | string | `36h` | CheckM wall-time request. |
| [`--checkm_cpus`](#checkm_cpus) | integer | `32` | CheckM CPU request. |
| [`--checkm_mem`](#checkm_mem) | string | `64GB` | CheckM memory request. |
| [`--run_autocycler`](#run_autocycler) | boolean | `false` | Run the AutoCycler consensus-assembly workflow. |
| [`--ac_enable_contig_filtering`](#ac_enable_contig_filtering) | boolean | `true` | Filter AutoCycler input contigs before clustering. |
| [`--ac_min_contig_len`](#ac_min_contig_len) | integer | `1000` | Minimum AutoCycler contig length. |
| [`--ac_max_contigs_per_assembly`](#ac_max_contigs_per_assembly) | integer | `25` | Maximum contigs retained from each assembly for AutoCycler. |
| [`--skip_subsample`](#skip_subsample) | boolean | `true` | Skip AutoCycler read subsampling. |
| [`--subsample_count`](#subsample_count) | integer | `6` | Number of long-read subsamples to create. |
| [`--min_read_depth`](#min_read_depth) | integer | `25` | Minimum read depth used during subsampling. |
| [`--autocycler_assemblers`](#autocycler_assemblers) | string | `flye,miniasm,necat,myloasm` | Comma-separated assembler set used by AutoCycler. |
| [`--autocycler_subsamples`](#autocycler_subsamples) | integer | `4` | Number of subsamples passed through the AutoCycler assembler set. |
| [`--use_flye`](#use_flye) | boolean | `true` | Include Flye in AutoCycler assembly generation. |
| [`--use_miniasm`](#use_miniasm) | boolean | `false` | Include Miniasm in AutoCycler assembly generation. |
| [`--use_necat`](#use_necat) | boolean | `false` | Include NECAT in AutoCycler assembly generation. |
| [`--use_myloasm`](#use_myloasm) | boolean | `true` | Include Myloasm in AutoCycler assembly generation. |

<a id="run_checkm"></a>
### `--run_checkm`

Run CheckM, CheckM2, and CheckV assembly validation where configured.

- Type: `boolean`
- Default: `true`

<a id="checkm_time"></a>
### `--checkm_time`

CheckM wall-time request.

- Type: `string`
- Default: `36h`

<a id="checkm_cpus"></a>
### `--checkm_cpus`

CheckM CPU request.

- Type: `integer`
- Default: `32`

<a id="checkm_mem"></a>
### `--checkm_mem`

CheckM memory request.

- Type: `string`
- Default: `64GB`

<a id="run_autocycler"></a>
### `--run_autocycler`

Run the AutoCycler consensus-assembly workflow.

- Type: `boolean`
- Default: `false`

<a id="ac_enable_contig_filtering"></a>
### `--ac_enable_contig_filtering`

Filter AutoCycler input contigs before clustering.

- Type: `boolean`
- Default: `true`

<a id="ac_min_contig_len"></a>
### `--ac_min_contig_len`

Minimum AutoCycler contig length.

- Type: `integer`
- Default: `1000`

<a id="ac_max_contigs_per_assembly"></a>
### `--ac_max_contigs_per_assembly`

Maximum contigs retained from each assembly for AutoCycler.

- Type: `integer`
- Default: `25`

<a id="skip_subsample"></a>
### `--skip_subsample`

Skip AutoCycler read subsampling.

- Type: `boolean`
- Default: `true`

<a id="subsample_count"></a>
### `--subsample_count`

Number of long-read subsamples to create.

- Type: `integer`
- Default: `6`

<a id="min_read_depth"></a>
### `--min_read_depth`

Minimum read depth used during subsampling.

- Type: `integer`
- Default: `25`

<a id="autocycler_assemblers"></a>
### `--autocycler_assemblers`

Comma-separated assembler set used by AutoCycler.

- Type: `string`
- Default: `flye,miniasm,necat,myloasm`

<a id="autocycler_subsamples"></a>
### `--autocycler_subsamples`

Number of subsamples passed through the AutoCycler assembler set.

- Type: `integer`
- Default: `4`

<a id="use_flye"></a>
### `--use_flye`

Include Flye in AutoCycler assembly generation.

- Type: `boolean`
- Default: `true`

<a id="use_miniasm"></a>
### `--use_miniasm`

Include Miniasm in AutoCycler assembly generation.

- Type: `boolean`
- Default: `false`

<a id="use_necat"></a>
### `--use_necat`

Include NECAT in AutoCycler assembly generation.

- Type: `boolean`
- Default: `false`

<a id="use_myloasm"></a>
### `--use_myloasm`

Include Myloasm in AutoCycler assembly generation.

- Type: `boolean`
- Default: `true`

## Sequence search and viral analysis

BLASTX, DIAMOND, MMseqs2, Hecatomb, and VirusSeeker options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--run_blastx`](#run_blastx) | boolean | `false` | Run DIAMOND BLASTX analysis. |
| [`--diamond_args`](#diamond_args) | string | `--block-size 20 --iterate faster --index-chunks 1 --tmpdir ${params.tmp_dir}` | Additional DIAMOND search options. |
| [`--run_mmseqs`](#run_mmseqs) | boolean | `false` | Run MMseqs2 nucleotide/protein search workflows. |
| [`--mmseqs_search_max_seqs`](#mmseqs_search_max_seqs) | integer | `20` | Maximum number of MMseqs2 target sequences retained per query. |
| [`--mmseqs_search_opts`](#mmseqs_search_opts) | string | `--max-seqs 25 -s 7.0 -e 1.0E-8 --search-type 3 --split 0 --cov-mode 2 --cov 0.15` | Additional options passed to `mmseqs search`. |
| [`--mmseqs_createdb_opts`](#mmseqs_createdb_opts) | string | `--createdb-mode 1` | Additional options passed to `mmseqs createdb`. |
| [`--mmseqs_format`](#mmseqs_format) | string | `query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen` | Column list passed to MMseqs2 `convertalis`. |
| [`--mmseqs_convertalis_opts`](#mmseqs_convertalis_opts) | string | `` | Additional options passed to `mmseqs convertalis`. |
| [`--mmseqs_cpus`](#mmseqs_cpus) | integer | `254` | CPU request for MMseqs2 searches. |
| [`--mmseqs_mem`](#mmseqs_mem) | integer | `960` | Memory request for MMseqs2 searches. |
| [`--db_BP`](#db_BP) | string | `` | Optional database-size value used by search calculations. |
| [`--hecatomb_config`](#hecatomb_config) | string | `` | Optional Hecatomb configuration file. |
| [`--vs`](#vs) | boolean | `false` | Run the VirusSeeker workflow. |

<a id="run_blastx"></a>
### `--run_blastx`

Run DIAMOND BLASTX analysis.

- Type: `boolean`
- Default: `false`

<a id="diamond_args"></a>
### `--diamond_args`

Additional DIAMOND search options.

- Type: `string`
- Default: `--block-size 20 --iterate faster --index-chunks 1 --tmpdir ${params.tmp_dir}`

<a id="run_mmseqs"></a>
### `--run_mmseqs`

Run MMseqs2 nucleotide/protein search workflows.

- Type: `boolean`
- Default: `false`

<a id="mmseqs_search_max_seqs"></a>
### `--mmseqs_search_max_seqs`

Maximum number of MMseqs2 target sequences retained per query.

- Type: `integer`
- Default: `20`

<a id="mmseqs_search_opts"></a>
### `--mmseqs_search_opts`

Additional options passed to `mmseqs search`.

- Type: `string`
- Default: `--max-seqs 25 -s 7.0 -e 1.0E-8 --search-type 3 --split 0 --cov-mode 2 --cov 0.15`

<a id="mmseqs_createdb_opts"></a>
### `--mmseqs_createdb_opts`

Additional options passed to `mmseqs createdb`.

- Type: `string`
- Default: `--createdb-mode 1`

<a id="mmseqs_format"></a>
### `--mmseqs_format`

Column list passed to MMseqs2 `convertalis`.

- Type: `string`
- Default: `query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen`

<a id="mmseqs_convertalis_opts"></a>
### `--mmseqs_convertalis_opts`

Additional options passed to `mmseqs convertalis`.

- Type: `string`
- Default: ``

<a id="mmseqs_cpus"></a>
### `--mmseqs_cpus`

CPU request for MMseqs2 searches.

- Type: `integer`
- Default: `254`

<a id="mmseqs_mem"></a>
### `--mmseqs_mem`

Memory request for MMseqs2 searches.

- Type: `integer`
- Default: `960`

<a id="db_BP"></a>
### `--db_BP`

Optional database-size value used by search calculations.

- Type: `string`
- Default: ``

<a id="hecatomb_config"></a>
### `--hecatomb_config`

Optional Hecatomb configuration file.

- Type: `string`
- Default: ``

<a id="vs"></a>
### `--vs`

Run the VirusSeeker workflow.

- Type: `boolean`
- Default: `false`

## Read taxonomic classification

Kraken2, Bracken, Sourmash, Mash, MetaPhlAn, GOTTCHA, and Taxpasta options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--run_reads_taxonomic_classifier`](#run_reads_taxonomic_classifier) | boolean | `false` | Run the combined read-taxonomy workflow. |
| [`--run_kraken2`](#run_kraken2) | boolean | `true` | Run Kraken2 classification. |
| [`--run_sourmash`](#run_sourmash) | boolean | `true` | Run Sourmash sketching and gather. |
| [`--run_mash`](#run_mash) | boolean | `true` | Run Mash distance analysis. |
| [`--run_metaphlan4`](#run_metaphlan4) | boolean | `false` | Run MetaPhlAn 4 classification. |
| [`--run_gottcha`](#run_gottcha) | boolean | `false` | Run GOTTCHA classification. |
| [`--run_taxpasta`](#run_taxpasta) | boolean | `true` | Normalize taxonomy outputs with Taxpasta. |
| [`--make_taxonomic_classifier_summary`](#make_taxonomic_classifier_summary) | boolean | `true` | Create a combined taxonomy summary across enabled classifiers. |
| [`--run_sourmash_kreport`](#run_sourmash_kreport) | boolean | `false` | Convert Sourmash results into Kraken-style reports. |
| [`--sourmash_taxonomy`](#sourmash_taxonomy) | string | `` | Optional Sourmash taxonomy mapping file. |
| [`--sourmash_tax_extra_args`](#sourmash_tax_extra_args) | string | `--ignore-abund` | Additional arguments for Sourmash taxonomy conversion. |
| [`--kraken2_extra_args`](#kraken2_extra_args) | string | `--use-names --confidence 0.05` | Additional Kraken2 command-line options. |
| [`--run_bracken`](#run_bracken) | boolean | `true` | Estimate species abundance from Kraken2 results using Bracken. |
| [`--bracken_level`](#bracken_level) | string | `S` | Bracken taxonomic level, such as `S` for species. Allowed values: `D`, `P`, `C`, `O`, `F`, `G`, `S`. |
| [`--bracken_threshold`](#bracken_threshold) | integer | `1` | Minimum reads required by Bracken. |
| [`--bracken_read_length_short`](#bracken_read_length_short) | integer | `150` | Expected short-read length supplied to Bracken. |
| [`--bracken_read_length_long`](#bracken_read_length_long) | integer | `300` | Expected long-read length supplied to Bracken. |
| [`--bracken_extra_args`](#bracken_extra_args) | string | `` | Additional Bracken command-line options. |
| [`--sourmash_scaled_param`](#sourmash_scaled_param) | string | `k=21,scaled=1000` | Sourmash sketch parameters. |
| [`--sourmash_extra_args`](#sourmash_extra_args) | string | `` | Additional Sourmash command-line options. |
| [`--metaphlan4_extra_args`](#metaphlan4_extra_args) | string | `` | Additional MetaPhlAn 4 command-line options. |
| [`--metaphlan4_db_arg`](#metaphlan4_db_arg) | string | `--db_dir` | Argument name used to supply the MetaPhlAn database directory. |
| [`--gottcha_extra_args`](#gottcha_extra_args) | string | `` | Additional GOTTCHA command-line options. |

<a id="run_reads_taxonomic_classifier"></a>
### `--run_reads_taxonomic_classifier`

Run the combined read-taxonomy workflow.

- Type: `boolean`
- Default: `false`

<a id="run_kraken2"></a>
### `--run_kraken2`

Run Kraken2 classification.

- Type: `boolean`
- Default: `true`

<a id="run_sourmash"></a>
### `--run_sourmash`

Run Sourmash sketching and gather.

- Type: `boolean`
- Default: `true`

<a id="run_mash"></a>
### `--run_mash`

Run Mash distance analysis.

- Type: `boolean`
- Default: `true`

<a id="run_metaphlan4"></a>
### `--run_metaphlan4`

Run MetaPhlAn 4 classification.

- Type: `boolean`
- Default: `false`

<a id="run_gottcha"></a>
### `--run_gottcha`

Run GOTTCHA classification.

- Type: `boolean`
- Default: `false`

<a id="run_taxpasta"></a>
### `--run_taxpasta`

Normalize taxonomy outputs with Taxpasta.

- Type: `boolean`
- Default: `true`

<a id="make_taxonomic_classifier_summary"></a>
### `--make_taxonomic_classifier_summary`

Create a combined taxonomy summary across enabled classifiers.

- Type: `boolean`
- Default: `true`

<a id="run_sourmash_kreport"></a>
### `--run_sourmash_kreport`

Convert Sourmash results into Kraken-style reports.

- Type: `boolean`
- Default: `false`

<a id="sourmash_taxonomy"></a>
### `--sourmash_taxonomy`

Optional Sourmash taxonomy mapping file.

- Type: `string`
- Default: ``

<a id="sourmash_tax_extra_args"></a>
### `--sourmash_tax_extra_args`

Additional arguments for Sourmash taxonomy conversion.

- Type: `string`
- Default: `--ignore-abund`

<a id="kraken2_extra_args"></a>
### `--kraken2_extra_args`

Additional Kraken2 command-line options.

- Type: `string`
- Default: `--use-names --confidence 0.05`

<a id="run_bracken"></a>
### `--run_bracken`

Estimate species abundance from Kraken2 results using Bracken.

- Type: `boolean`
- Default: `true`

<a id="bracken_level"></a>
### `--bracken_level`

Bracken taxonomic level, such as `S` for species.

- Type: `string`
- Default: `S`
- Allowed values: `D`, `P`, `C`, `O`, `F`, `G`, `S`

<a id="bracken_threshold"></a>
### `--bracken_threshold`

Minimum reads required by Bracken.

- Type: `integer`
- Default: `1`

<a id="bracken_read_length_short"></a>
### `--bracken_read_length_short`

Expected short-read length supplied to Bracken.

- Type: `integer`
- Default: `150`

<a id="bracken_read_length_long"></a>
### `--bracken_read_length_long`

Expected long-read length supplied to Bracken.

- Type: `integer`
- Default: `300`

<a id="bracken_extra_args"></a>
### `--bracken_extra_args`

Additional Bracken command-line options.

- Type: `string`
- Default: ``

<a id="sourmash_scaled_param"></a>
### `--sourmash_scaled_param`

Sourmash sketch parameters.

- Type: `string`
- Default: `k=21,scaled=1000`

<a id="sourmash_extra_args"></a>
### `--sourmash_extra_args`

Additional Sourmash command-line options.

- Type: `string`
- Default: ``

<a id="metaphlan4_extra_args"></a>
### `--metaphlan4_extra_args`

Additional MetaPhlAn 4 command-line options.

- Type: `string`
- Default: ``

<a id="metaphlan4_db_arg"></a>
### `--metaphlan4_db_arg`

Argument name used to supply the MetaPhlAn database directory.

- Type: `string`
- Default: `--db_dir`

<a id="gottcha_extra_args"></a>
### `--gottcha_extra_args`

Additional GOTTCHA command-line options.

- Type: `string`
- Default: ``

## Virulence-factor classifier

Machine-learning and VFDB homology settings for VF classification.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--vf_classifier`](#vf_classifier) | boolean | `false` | Run the virulence-factor classifier independently of VirusSeeker. |
| [`--vf_input_type`](#vf_input_type) | string | `nucleotide` | VF-classifier input type: `nucleotide` or `protein`. Allowed values: `nucleotide`, `protein`. |
| [`--vf_mode`](#vf_mode) | string | `isolate` | VF-classifier analysis mode: `isolate` or `metagenome`. Allowed values: `isolate`, `metagenome`. |
| [`--vf_homology`](#vf_homology) | string | `yes` | Enable or disable VFDB homology support using `yes` or `no`. Allowed values: `yes`, `no`. |
| [`--vf_gpu`](#vf_gpu) | boolean | `true` | Expose a GPU to the VF-classifier container and request a GPU from the scheduler. |
| [`--vf_threshold`](#vf_threshold) | number | `0.34` | Model probability threshold used to call virulence factors. |

<a id="vf_classifier"></a>
### `--vf_classifier`

Run the virulence-factor classifier independently of VirusSeeker.

- Type: `boolean`
- Default: `false`

<a id="vf_input_type"></a>
### `--vf_input_type`

VF-classifier input type: `nucleotide` or `protein`.

- Type: `string`
- Default: `nucleotide`
- Allowed values: `nucleotide`, `protein`

<a id="vf_mode"></a>
### `--vf_mode`

VF-classifier analysis mode: `isolate` or `metagenome`.

- Type: `string`
- Default: `isolate`
- Allowed values: `isolate`, `metagenome`

<a id="vf_homology"></a>
### `--vf_homology`

Enable or disable VFDB homology support using `yes` or `no`.

- Type: `string`
- Default: `yes`
- Allowed values: `yes`, `no`

<a id="vf_gpu"></a>
### `--vf_gpu`

Expose a GPU to the VF-classifier container and request a GPU from the scheduler.

- Type: `boolean`
- Default: `true`

<a id="vf_threshold"></a>
### `--vf_threshold`

Model probability threshold used to call virulence factors.

- Type: `number`
- Default: `0.34`

## Contig characterization

Annotation, AMR/VF, typing, plasmid, prophage, and chimeric-contig analyses.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--characterize_contigs`](#characterize_contigs) | boolean | `false` | Run the selected contig-characterization tools. |
| [`--prokka`](#prokka) | boolean | `false` | Annotate contigs with Prokka. |
| [`--busco`](#busco) | boolean | `false` | Evaluate assembly completeness with BUSCO. |
| [`--amr_vf`](#amr_vf) | boolean | `false` | Run AMR and virulence-factor similarity searches. |
| [`--mlst`](#mlst) | boolean | `false` | Run multilocus sequence typing. |
| [`--ge_screen`](#ge_screen) | boolean | `false` | Run the genetic-engineering screening workflow. |
| [`--rgi`](#rgi) | boolean | `false` | Run the CARD Resistance Gene Identifier. |
| [`--amrfinder`](#amrfinder) | boolean | `false` | Run NCBI AMRFinderPlus. |
| [`--chimeric_detection`](#chimeric_detection) | boolean | `false` | Run chimeric-contig detection. |
| [`--blast_contigs`](#blast_contigs) | boolean | `true` | Run nucleotide BLAST analysis of assembled contigs. |
| [`--run_phispy`](#run_phispy) | boolean | `false` | Run PhiSpy prophage prediction. |
| [`--run_mobsuite`](#run_mobsuite) | boolean | `false` | Run MOB-suite plasmid analysis. |
| [`--run_plasme`](#run_plasme) | boolean | `false` | Run PLASMe plasmid prediction. |
| [`--rgi_load_cmd`](#rgi_load_cmd) | string | `` | Optional RGI database-loading command. |
| [`--plasme_mode`](#plasme_mode) | string | `balance` | PLASMe sensitivity/specificity mode. Allowed values: `specific`, `balance`, `sensitive`. |
| [`--desired_rank`](#desired_rank) | string | `genus` | Taxonomic rank targeted by chimeric-contig analysis. |
| [`--max_target_seqs`](#max_target_seqs) | integer | `10` | Maximum BLAST target sequences retained for chimeric-contig analysis. |

<a id="characterize_contigs"></a>
### `--characterize_contigs`

Run the selected contig-characterization tools.

- Type: `boolean`
- Default: `false`

<a id="prokka"></a>
### `--prokka`

Annotate contigs with Prokka.

- Type: `boolean`
- Default: `false`

<a id="busco"></a>
### `--busco`

Evaluate assembly completeness with BUSCO.

- Type: `boolean`
- Default: `false`

<a id="amr_vf"></a>
### `--amr_vf`

Run AMR and virulence-factor similarity searches.

- Type: `boolean`
- Default: `false`

<a id="mlst"></a>
### `--mlst`

Run multilocus sequence typing.

- Type: `boolean`
- Default: `false`

<a id="ge_screen"></a>
### `--ge_screen`

Run the genetic-engineering screening workflow.

- Type: `boolean`
- Default: `false`

<a id="rgi"></a>
### `--rgi`

Run the CARD Resistance Gene Identifier.

- Type: `boolean`
- Default: `false`

<a id="amrfinder"></a>
### `--amrfinder`

Run NCBI AMRFinderPlus.

- Type: `boolean`
- Default: `false`

<a id="chimeric_detection"></a>
### `--chimeric_detection`

Run chimeric-contig detection.

- Type: `boolean`
- Default: `false`

<a id="blast_contigs"></a>
### `--blast_contigs`

Run nucleotide BLAST analysis of assembled contigs.

- Type: `boolean`
- Default: `true`

<a id="run_phispy"></a>
### `--run_phispy`

Run PhiSpy prophage prediction.

- Type: `boolean`
- Default: `false`

<a id="run_mobsuite"></a>
### `--run_mobsuite`

Run MOB-suite plasmid analysis.

- Type: `boolean`
- Default: `false`

<a id="run_plasme"></a>
### `--run_plasme`

Run PLASMe plasmid prediction.

- Type: `boolean`
- Default: `false`

<a id="rgi_load_cmd"></a>
### `--rgi_load_cmd`

Optional RGI database-loading command.

- Type: `string`
- Default: ``

<a id="plasme_mode"></a>
### `--plasme_mode`

PLASMe sensitivity/specificity mode.

- Type: `string`
- Default: `balance`
- Allowed values: `specific`, `balance`, `sensitive`

<a id="desired_rank"></a>
### `--desired_rank`

Taxonomic rank targeted by chimeric-contig analysis.

- Type: `string`
- Default: `genus`

<a id="max_target_seqs"></a>
### `--max_target_seqs`

Maximum BLAST target sequences retained for chimeric-contig analysis.

- Type: `integer`
- Default: `10`

## Automated read mapping

Automatic reference selection, ANI filtering, and mapping options.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--run_readmapping_auto`](#run_readmapping_auto) | boolean | `false` | Automatically select references and map reads. |
| [`--readmapping_auto_cluster_dedup`](#readmapping_auto_cluster_dedup) | boolean | `false` | Cluster and deduplicate automatically selected references. |
| [`--readmapping_auto_ani_threshold`](#readmapping_auto_ani_threshold) | integer | `99` | ANI percentage required when deduplicating selected references. |

<a id="run_readmapping_auto"></a>
### `--run_readmapping_auto`

Automatically select references and map reads.

- Type: `boolean`
- Default: `false`

<a id="readmapping_auto_cluster_dedup"></a>
### `--readmapping_auto_cluster_dedup`

Cluster and deduplicate automatically selected references.

- Type: `boolean`
- Default: `false`

<a id="readmapping_auto_ani_threshold"></a>
### `--readmapping_auto_ani_threshold`

ANI percentage required when deduplicating selected references.

- Type: `integer`
- Default: `99`

## Krona rendering

Browser timing and image dimensions for automated Krona snapshots.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--krona_snapshot_width`](#krona_snapshot_width) | integer | `1600` | Browser viewport width used for Krona screenshots. |
| [`--krona_snapshot_height`](#krona_snapshot_height) | integer | `1600` | Browser viewport height used for Krona screenshots. |
| [`--krona_snapshot_wait_ms`](#krona_snapshot_wait_ms) | integer | `1600` | Delay before attempting a Krona screenshot. |
| [`--krona_snapshot_popup_timeout_ms`](#krona_snapshot_popup_timeout_ms) | integer | `30000` | Maximum time to wait for the Krona snapshot popup. |
| [`--krona_snapshot_download_timeout_ms`](#krona_snapshot_download_timeout_ms) | integer | `30000` | Maximum time to wait for a Krona image download. |
| [`--krona_snapshot_png_width`](#krona_snapshot_png_width) | integer | `1800` | Output PNG width for Krona snapshots. |
| [`--krona_font`](#krona_font) | integer | `15` | Krona label font size. |

<a id="krona_snapshot_width"></a>
### `--krona_snapshot_width`

Browser viewport width used for Krona screenshots.

- Type: `integer`
- Default: `1600`

<a id="krona_snapshot_height"></a>
### `--krona_snapshot_height`

Browser viewport height used for Krona screenshots.

- Type: `integer`
- Default: `1600`

<a id="krona_snapshot_wait_ms"></a>
### `--krona_snapshot_wait_ms`

Delay before attempting a Krona screenshot.

- Type: `integer`
- Default: `1600`

<a id="krona_snapshot_popup_timeout_ms"></a>
### `--krona_snapshot_popup_timeout_ms`

Maximum time to wait for the Krona snapshot popup.

- Type: `integer`
- Default: `30000`

<a id="krona_snapshot_download_timeout_ms"></a>
### `--krona_snapshot_download_timeout_ms`

Maximum time to wait for a Krona image download.

- Type: `integer`
- Default: `30000`

<a id="krona_snapshot_png_width"></a>
### `--krona_snapshot_png_width`

Output PNG width for Krona snapshots.

- Type: `integer`
- Default: `1800`

<a id="krona_font"></a>
### `--krona_font`

Krona label font size.

- Type: `integer`
- Default: `15`

## Advanced and developer options

Internal paths, experimental branches, and generated report identifiers.

| Parameter | Type | Default | Description |
|---|---|---|---|
| [`--test_unmapped_reassembly`](#test_unmapped_reassembly) | boolean | `false` | Enable the experimental unmapped-read analysis branch. |
| [`--scripts`](#scripts) | string | `${projectDir}/scripts` | Directory containing pipeline helper scripts. Normally follows the repository automatically. |
| [`--trace_report_suffix`](#trace_report_suffix) | string | Generated | Timestamp suffix used for trace, timeline, DAG, and execution-report filenames. |

<a id="test_unmapped_reassembly"></a>
### `--test_unmapped_reassembly`

Enable the experimental unmapped-read analysis branch.

- Type: `boolean`
- Default: `false`

<a id="scripts"></a>
### `--scripts`

Directory containing pipeline helper scripts. Normally follows the repository automatically.

- Type: `string`
- Default: `${projectDir}/scripts`

<a id="trace_report_suffix"></a>
### `--trace_report_suffix`

Timestamp suffix used for trace, timeline, DAG, and execution-report filenames.

- Type: `string`
- Default: Generated
