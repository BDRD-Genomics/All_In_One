# All-In-One Sequencing Pipeline

**All-In-One** is a modular [Nextflow](https://www.nextflow.io/) DSL2 workflow for metagenomic and viral sequencing analysis. It supports **paired short-read**, **long-read**, and **hybrid** sequencing data and combines quality control, host/rRNA removal, assembly, taxonomic classification, viral analysis, genome characterization, and reporting in a single configurable pipeline.

**Documentation:** Browse the [documentation](docs/index.md), jump to the [parameter reference](docs/parameters.md), see the [database setup guide](docs/databases.md), or review the [Docker and Apptainer guide](docs/containers.md).

---

## Overview

The pipeline is designed so that major analysis stages can be enabled or disabled independently with Nextflow parameters. A run can therefore be configured for anything from basic read QC to a larger workflow containing assembly, taxonomic classification, AMR characterization, viral screening, and reference-based read filtering.

### Major workflow areas

| Workflow area                    | Included capabilities                                                    |
| -------------------------------- | ------------------------------------------------------------------------ |
| Quality control                  | fastp, FastQC, NanoPlot, MultiQC, read-distribution statistics           |
| Host / rRNA removal              | host read mapping, BBMap/minimap2, RiboDetector                          |
| Assembly                         | SPAdes, metaSPAdes, plasmidSPAdes, Unicycler, Dragonflye, Raven, Myloasm |
| Assembly validation              | CheckM, CheckM2, CheckV, BUSCO, and AutoCycler workflows                 |
| Taxonomic classification         | Kraken2/Bracken, Sourmash, Mash, MetaPhlAn, GOTTCHA, Taxpasta            |
| Sequence search / viral analysis | BLAST, DIAMOND BLASTX, MMseqs2, VirusSeeker workflows                    |
| Contig characterization          | Prokka, MLST, AMRFinderPlus, RGI, PhiSpy, MOB-suite, PLASMe              |
| Reporting                        | Nextflow trace, timeline, execution report, DAG, and QC summaries        |

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
├── env/                    # Conda/Mamba environment definitions
├── databases/
│   └── manifest.tsv        # Database source and installation manifest
├── containers/
│   ├── docker/             # Docker image build definitions
│   └── apptainer/          # Offline/local Apptainer resources
├── ICTV/                   # ICTV viral family/genome-size resources
└── docs/
    ├── index.md            # Documentation home
    ├── usage.md            # Usage guide
    ├── parameters.md       # Complete parameter reference
    ├── databases.md        # Database installation and offline setup
    ├── containers.md       # Docker, Apptainer, and offline container guide
    └── configuration.md    # Installation and site configuration
```

---

## Requirements

At minimum, a system running All-In-One requires:

* Nextflow with DSL2 support
* Java compatible with the installed Nextflow release
* Linux
* One configured software backend:

  * Apptainer
  * Docker
  * Conda/Mamba
* Required reference databases for the workflows being enabled

For cluster execution, the repository includes a **Slurm** profile in `conf/profiles.config`.

---

## Database setup

All-In-One uses a single database root for shared reference data.

For example:

```bash
export AIO_DATABASE_DIR=/data/All_In_One_databases
```

The same location can be supplied to the pipeline with:

```bash
--database_dir /data/All_In_One_databases
```

Database provisioning is intentionally performed separately from Nextflow execution. This prevents large downloads such as NCBI `nt` and `nr` from unexpectedly starting inside a local or Slurm analysis job.

### Automated database installation

The repository includes:

```text
scripts/download_databases.sh
scripts/verify_databases.sh
databases/manifest.tsv
```

List supported database installers:

```bash
./scripts/download_databases.sh --list
```

Install the standard database set:

```bash
./scripts/download_databases.sh \
    --db-root "$AIO_DATABASE_DIR" \
    --preset standard
```

Install selected databases only:

```bash
./scripts/download_databases.sh \
    --db-root "$AIO_DATABASE_DIR" \
    --db checkm2,checkv,taxdb
```

Install the full database set, including NCBI `nt` and `nr`:

```bash
./scripts/download_databases.sh \
    --db-root "$AIO_DATABASE_DIR" \
    --preset full
```

### Database presets

| Preset     | Databases                                                             |
| ---------- | --------------------------------------------------------------------- |
| `minimal`  | `taxdb`, `checkm2`, `checkv`                                          |
| `standard` | minimal + `checkm`, `busco`, `kraken2`, `metaphlan4`, `amrfinderplus` |
| `full`     | standard + `nt`, `nr`                                                 |

NCBI `nt` and `nr` are intentionally excluded from the standard preset because of their size.

The installer delegates to upstream database installers whenever possible, including `update_blastdb.pl`, CheckM2, CheckV, BUSCO, Kraken2, MetaPhlAn, and AMRFinderPlus.

Verify an installed database root with:

```bash
./scripts/verify_databases.sh \
    --db-root "$AIO_DATABASE_DIR"
```

For database sources, expected directory layout, manual installation, updates, and air-gapped deployment, see:

**[Database setup guide](docs/databases.md)**

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

Edit `site.config` for the local installation. Installation-specific settings such as database locations, temporary storage, container directories, and other site resources can be defined there.

`site.config` is intended to be installation-specific and should not be committed when it contains local paths or credentials.

A typical database root can also be exported directly:

```bash
export AIO_DATABASE_DIR=/data/All_In_One_databases
```

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
    --run_id example_run \
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
    --run_id example_run \
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
    --run_id example_run \
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
-profile local,conda
```

```bash
-profile slurm,conda
```

```bash
-profile local,docker
```

The available execution/backend profiles are defined in [`conf/profiles.config`](conf/profiles.config).

### Conda/Mamba

The Conda profile uses repository-managed environments from:

```text
env/
```

Mamba is used as the environment solver for faster dependency resolution.

For example:

```bash
nextflow run main.nf \
    -c site.config \
    -profile local,conda \
    --samplesheet samplesheet.csv \
    --run_id example_run \
    ...
```

### RiboDetector CPU/GPU support

RiboDetector supports both CPU and GPU execution. GPU support depends on the selected software backend:

| Backend     |  CPU mode |                                                        GPU mode |
| ----------- | --------: | --------------------------------------------------------------: |
| Conda/Mamba | Supported | Supported and tested with the provided CUDA-enabled environment |
| Apptainer   | Supported |                Supported and tested with NVIDIA GPU passthrough |
| Docker      | Supported |              Not supported in the current All-In-One deployment |

Select the RiboDetector execution mode with:

```bash
--ribodetector_mode cpu
```

or:

```bash
--ribodetector_mode gpu
```

For Conda/Mamba GPU execution, `env/environment.ribodetector.yml` provides a CUDA-enabled PyTorch environment. The host system must provide a compatible NVIDIA driver.

For Apptainer GPU execution, NVIDIA GPU passthrough is required. On supported systems, this is provided using Apptainer's `--nv` option.

Docker should be used with RiboDetector CPU mode in the current All-In-One deployment. Use the Conda/Mamba or Apptainer backend when GPU execution is required.

---

## Common workflow switches

The following are some of the primary switches used to control a run:

| Parameter                          | Purpose                                               |
| ---------------------------------- | ----------------------------------------------------- |
| `--run_qc`                         | Run read quality control                              |
| `--map2reads`                      | Remove reads mapping to the configured host reference |
| `--remove_rRNA_reads`              | Remove rRNA reads                                     |
| `--ribodetector_mode`              | Select `cpu` or `gpu` RiboDetector execution          |
| `--run_assembly`                   | Enable genome/metagenome assembly                     |
| `--run_blastx`                     | Run DIAMOND BLASTX analysis                           |
| `--run_mmseqs`                     | Run MMseqs2 searches                                  |
| `--run_reads_taxonomic_classifier` | Enable read-level taxonomic classification            |
| `--characterize_contigs`           | Enable contig characterization workflows              |
| `--vs`                             | Enable VirusSeeker workflow components                |
| `--run_autocycler`                 | Enable AutoCycler workflow                            |

For every option, its default, type, allowed values, and related settings, see the **[complete parameter reference](docs/parameters.md)**.

---

## Configuration

Configuration is separated by purpose under `conf/`:

| File                       | Purpose                                               |
| -------------------------- | ----------------------------------------------------- |
| `conf/profiles.config`     | Nextflow execution and software backend profiles      |
| `conf/resources.config`    | CPU, memory, time, and process resource configuration |
| `conf/containers.config`   | Docker/Apptainer runtime image configuration          |
| `conf/conda.config`        | Conda/Mamba environment configuration                 |
| `conf/databases.config`    | Reference database locations                          |
| `conf/site.config.example` | Template for installation-specific overrides          |

User-facing analysis switches and defaults are maintained in [`params.config`](params.config).

### Project directory vs. run ID

Nextflow's built-in:

```text
projectDir
```

refers to the root directory of the cloned All-In-One repository.

The pipeline parameter:

```text
--run_id
```

identifies an individual analysis run and is used in output directory and report naming.

For example:

```text
projectDir
/home/user/All_In_One

run_id
example_run
```

These values are independent.

---

## Containers

Container build definitions are stored in:

```text
containers/docker/
```

Docker images are published through the project's container registry and can also be consumed by Apptainer using `docker://` image references.

For online Apptainer execution, images can be converted and cached automatically from the Docker registry.

For offline or air-gapped systems, prebuilt `.sif` files can be staged under the configured container directory.

Recommended Apptainer cache locations on systems with limited root filesystem capacity can be configured with:

```bash
export APPTAINER_TMPDIR=/data/tmp/apptainer_tmp
export APPTAINER_CACHEDIR=/data/tmp/apptainer_cache
```

See the **[Docker and Apptainer documentation](docs/containers.md)** for backend selection, image handling, SIF conversion, and offline HPC deployment.

---

## Offline / air-gapped deployment

All-In-One is designed to support disconnected HPC environments.

For offline deployment, prepare the following on an internet-connected system:

1. Clone or package the All-In-One repository.
2. Download required container images or build the required Apptainer `.sif` files.
3. Download the required databases using `scripts/download_databases.sh`.
4. Transfer the repository, containers, and database tree to the offline system.
5. Configure `site.config`, `AIO_DATABASE_DIR`, and the local container directory.

Database archives can be staged with:

```bash
tar -I 'pigz -1' \
    -cf All_In_One_databases.tar.gz \
    -C "$(dirname "$AIO_DATABASE_DIR")" \
    "$(basename "$AIO_DATABASE_DIR")"
```

See:

* [Database setup](docs/databases.md)
* [Container setup](docs/containers.md)
* [Configuration](docs/configuration.md)

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

Individual analysis outputs are organized under the selected run identifier:

```text
<outdir>/<run_id>/
```

For example:

```text
results/example_run/
```

Nextflow execution metadata is also generated, including:

* execution report
* timeline report
* trace file
* workflow DAG

Individual workflows create their own subdirectories beneath the configured run output.

---

## Parameter reference

The complete parameter documentation is maintained separately so the GitHub landing page stays readable.

### **[View all pipeline parameters →](docs/parameters.md)**

The parameter page is organized into sections for:

* Input and output
* Sequencing mode
* Quality control
* Host and rRNA removal
* Assembly
* Assembly validation and AutoCycler
* Sequence search and viral analysis
* Read taxonomic classification
* Contig characterization
* Krona rendering
* Advanced/developer options

---

## Database reference

The complete database installation and deployment documentation is available at:

### **[Database setup and sources →](docs/databases.md)**

This guide includes:

* database installation presets
* official/upstream database sources
* automated download commands
* expected directory layout
* database verification
* update behavior
* manual installation
* offline and air-gapped deployment

---

## Development

The pipeline uses Nextflow DSL2 with reusable modules under:

```text
workflows/modules/local/
```

When adding or changing a user-facing parameter, update:

```text
params.config
nextflow_schema.json
docs/parameters.md
```

When adding or changing a database dependency, update:

```text
conf/databases.config
databases/manifest.tsv
scripts/download_databases.sh
docs/databases.md
```

When adding or changing a software dependency, update the appropriate:

```text
env/environment.<tool>.yml
conf/conda.config
conf/containers.config
```

and confirm the relevant execution profiles remain valid.

---

## License

