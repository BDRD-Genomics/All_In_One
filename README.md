# All-In-One Sequencing Pipeline

**All-In-One** is a modular [Nextflow](https://www.nextflow.io/) DSL2 workflow for metagenomic and viral sequencing analysis. It supports **paired short-read**, **long-read**, and **hybrid** sequencing data and combines quality control, host/rRNA removal, assembly, taxonomic classification, viral analysis, genome characterization, automated reference mapping, and reporting in a single configurable pipeline.

**Documentation:** Browse the [documentation](docs/index.md), jump to the [parameter reference](docs/parameters.md), or see the [Docker and Apptainer guide](docs/containers.md).

The Markdown parameter reference is available directly in the repository at [`docs/parameters.md`](docs/parameters.md).

---

## Overview

The pipeline is designed so that major analysis stages can be enabled or disabled independently with Nextflow parameters. A run can therefore be configured for anything from basic read QC to a larger workflow containing assembly, taxonomic classification, AMR characterization, viral screening, and reference-based read mapping.

### Major workflow areas

| Workflow area | Included capabilities |
|---|---|
| Quality control | fastp, FastQC, NanoPlot, MultiQC, read-distribution statistics |
| Host / rRNA removal | host read mapping, BBMap/minimap2, RiboDetector |
| Assembly | SPAdes, metaSPAdes, plasmidSPAdes, Unicycler, Dragonflye, Raven, Myloasm |
| Assembly validation | CheckM / CheckM2 / CheckV and AutoCycler workflows |
| Taxonomic classification | Kraken2/Bracken, Sourmash, Mash, MetaPhlAn, GOTTCHA, Taxpasta |
| Sequence search / viral analysis | DIAMOND BLASTX, MMseqs2, VirusSeeker workflows |
| Contig characterization | Prokka, BUSCO, MLST, AMRFinderPlus, RGI, PhiSpy, MOB-suite, PLASMe |
| Reporting | Nextflow trace, timeline, report, DAG, QC summaries, and exercise reports |

---

## Repository layout

```text
All_In_One/
├── main.nf                 # Main DSL2 workflow entry point
├── nextflow.config         # Top-level Nextflow configuration
├── nextflow_schema.json    # Nextflow parameter schema
├── params.config           # User-facing pipeline parameters
├── conf/                   # Profiles, resources, containers, databases, site config
├── workflows/              # Workflows and local DSL2 modules
├── scripts/                # Supporting Python, R, and shell scripts
├── env/                    # Conda environment definitions
├── containers/
│   ├── docker/             # Dockerfiles and container build environments
│   └── apptainer/          # Apptainer/SIF build scripts and offline container setup
├── ICTV/                   # ICTV viral family/genome-size resources
└── docs/
    ├── index.md            # Documentation home
    ├── usage.md            # Usage guide
    ├── parameters.md       # Complete parameter reference
    ├── containers.md       # Docker, Apptainer, and offline container guide
    └── configuration.md    # Installation and site configuration
```

---

## Requirements

At minimum, a system running All-In-One requires:

- Nextflow with DSL2 support
- Java compatible with the installed Nextflow release
- Linux
- One configured software backend:
  - Apptainer
  - Docker
  - Conda/Mamba
- Required reference databases for the workflows being enabled

For cluster execution, the repository includes a **Slurm** profile in `conf/profiles.config`.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/BDRD-Genomics/All_In_One.git
cd All_In_One
```

Create a site-specific configuration from the included template:

```bash
cp conf/site.config.example site.config
```

Edit `site.config` for the local installation. This is where paths such as container directories, database directories, temporary storage, MEGAN, and other installation-specific resources should be defined.

`site.config` is intended to be installation-specific and should not be committed when it contains local paths or credentials.

---

## Input sample sheet

The pipeline reads a CSV sample sheet with the following columns:

```text
sample_id,fastq_1,fastq_2,long_read
```

Use `No_Read` for a sequencing input that is not present for that sample.

### Short-read example

```csv
sample_id,fastq_1,fastq_2,long_read
sample01,/data/sample01_R1.fastq.gz,/data/sample01_R2.fastq.gz,No_Read
```

### Long-read example

```csv
sample_id,fastq_1,fastq_2,long_read
sample01,No_Read,No_Read,/data/sample01.fastq.gz
```

### Hybrid example

```csv
sample_id,fastq_1,fastq_2,long_read
sample01,/data/sample01_R1.fastq.gz,/data/sample01_R2.fastq.gz,/data/sample01.fastq.gz
```

Select the corresponding input mode with one of:

```text
--shortreads true
--longreads true
--hybrid true
```

Only the sequencing mode appropriate for the run should be enabled.

---

## Quick start

### Short-read QC

```bash
nextflow run main.nf \
    -c site.config \
    -profile slurm,apptainer \
    --samplesheet samplesheet.csv \
    --project_id example_project \
    --shortreads true \
    --run_qc true \
    --outdir results
```

### Hybrid assembly

```bash
nextflow run main.nf \
    -c site.config \
    -profile slurm,apptainer \
    --samplesheet samplesheet.csv \
    --project_id example_project \
    --hybrid true \
    --run_qc true \
    --run_assembly true \
    --metaspades true \
    --outdir results
```

### Long-read assembly with Myloasm

```bash
nextflow run main.nf \
    -c site.config \
    -profile slurm,apptainer \
    --samplesheet samplesheet.csv \
    --project_id example_project \
    --longreads true \
    --run_qc true \
    --run_assembly true \
    --myloasm true \
    --outdir results
```

---

## Execution profiles

Execution and software backends are separated so profiles can be combined.

Examples:

```bash
-profile local,apptainer
```

```bash
-profile slurm,apptainer
```

```bash
-profile slurm,conda
```

```bash
-profile local,docker
```

The available execution/backend profiles are defined in [`conf/profiles.config`](conf/profiles.config).

---

## Common workflow switches

The following are some of the primary switches used to control a run:

| Parameter | Purpose |
|---|---|
| `--run_qc` | Run read quality control |
| `--map2reads` | Remove reads mapping to the configured host reference |
| `--remove_rRNA_reads` | Remove rRNA reads |
| `--run_assembly` | Enable genome/metagenome assembly |
| `--run_blastx` | Run DIAMOND BLASTX analysis |
| `--run_mmseqs` | Run MMseqs2 searches |
| `--run_reads_taxonomic_classifier` | Enable read-level taxonomic classification |
| `--characterize_contigs` | Enable contig characterization workflows |
| `--vs` | Enable VirusSeeker workflow components |
| `--run_autocycler` | Enable AutoCycler workflow |

For every option, its default, type, allowed values, and related settings, see the **[complete parameter reference](docs/parameters.md)**.

---

## Configuration

Configuration is separated by purpose under `conf/`:

| File | Purpose |
|---|---|
| `conf/profiles.config` | Nextflow execution and software backend profiles |
| `conf/resources.config` | CPU, memory, time, and process resource configuration |
| `conf/containers.config` | Apptainer/container image locations |
| `conf/conda.config` | Conda environment configuration |
| `conf/databases.config` | Database locations |
| `conf/databases.config.internal` | Internal database configuration |
| `conf/site.config.example` | Template for installation-specific overrides |

User-facing analysis switches and defaults are maintained in [`params.config`](params.config).

---

## Containers

Container build definitions are stored in:

```text
containers/docker/
```

This directory contains the Dockerfiles and Conda environment definitions used to build the pipeline's component images. Runtime container selection is configured through `conf/containers.config` and the selected Nextflow profile.

See the **[Docker and Apptainer documentation](https://bdrd-genomics.github.io/All_In_One/containers/)** for backend selection, image building, SIF conversion, and offline HPC deployment.

---

## Output

By default, pipeline output is written beneath:

```text
results/
```

unless overridden with:

```bash
--outdir /path/to/output
```

Nextflow execution metadata is also generated, including:

- execution report
- timeline report
- trace file
- workflow DAG

Individual analysis workflows write their results below the configured project/output structure.

---

## Parameter reference

The complete parameter documentation is maintained separately so the GitHub landing page stays readable.

### **[View all pipeline parameters →](docs/parameters.md)**

The parameter page is organized into sections for:

- Input and output
- Sequencing mode
- Quality control
- Host and rRNA removal
- Assembly
- Assembly validation and AutoCycler
- Sequence search and viral analysis
- Read taxonomic classification
- Contig characterization
- Krona rendering
- Advanced/developer options

---

## Development

The pipeline uses Nextflow DSL2 with reusable modules under:

```text
workflows/modules/local/
```

When adding or changing a user-facing parameter, update both:

```text
params.config
docs/parameters.md
```

so the documented parameter reference remains synchronized with pipeline behavior.

---

## License
