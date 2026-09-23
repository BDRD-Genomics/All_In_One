# Containers

All-In-One supports containerized execution with **Apptainer/Singularity** and **Docker**. A Conda/Mamba backend is also available for environments where containers are not appropriate.

[Apptainer](#apptainer){ .md-button .md-button--primary }
[Docker](#docker){ .md-button }
[Configuration](configuration.md){ .md-button }

## Supported backends

| Backend | Typical use |
|---|---|
| **Apptainer / Singularity** | HPC and Slurm clusters, including offline or restricted systems |
| **Docker** | Local workstations, servers, image development, and image publishing |
| **Conda / Mamba** | Alternative software backend when containers are unavailable |

Execution and software backends are selected through Nextflow profiles defined in `conf/profiles.config`.

## Apptainer

Apptainer is the recommended container backend for Slurm/HPC execution.

### Slurm + Apptainer

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

### Local + Apptainer

```bash
nextflow run main.nf \
    -c site.config \
    -profile local,apptainer \
    --samplesheet samplesheet.csv \
    --project_id example_project \
    --shortreads true \
    --run_qc true \
    --outdir results
```

### Local SIF images

All-In-One can use locally stored `.sif` images. This is especially useful on HPC systems that do not have internet access at runtime.

Container mappings are maintained in:

```text
conf/containers.config
```

Installation-specific container roots and overrides should be placed in the local site configuration created from:

```text
conf/site.config.example
```

A local image directory can contain images such as:

```text
containers/apptainer/
├── allinone-aio_qc-1.0.1.sif
├── allinone-amrfinder_plus-1.0.1.sif
├── allinone-bbmap-1.0.1.sif
├── allinone-blast-1.0.1.sif
├── allinone-busco-1.0.1.sif
├── allinone-diamond-1.0.1.sif
├── allinone-fastqc-1.0.1.sif
├── allinone-myloasm-1.0.1.sif
└── ...
```

The exact image names used by the workflow are controlled by `conf/containers.config`.

## Docker

Docker image definitions are stored under:

```text
containers/docker/
```

The directory contains component-specific Dockerfiles and Conda environment definitions. Examples include:

```text
Dockerfile.aio_qc
Dockerfile.amrfinder_plus
Dockerfile.autocycler
Dockerfile.bbmap
Dockerfile.blast
Dockerfile.busco
Dockerfile.checkm-genome
Dockerfile.diamond
Dockerfile.dragonflye
Dockerfile.fastp_long
Dockerfile.fastqc
Dockerfile.kraken2
Dockerfile.mmseqs2
Dockerfile.myloasm
Dockerfile.nanoplot
Dockerfile.prokka
Dockerfile.rgi
Dockerfile.ribodetector
Dockerfile.spades
Dockerfile.unicycler
Dockerfile.vs
```

### Local + Docker

```bash
nextflow run main.nf \
    -c site.config \
    -profile local,docker \
    --samplesheet samplesheet.csv \
    --project_id example_project \
    --shortreads true \
    --run_qc true \
    --outdir results
```

### Container registry

Project images are designed to use the BDRD Genomics GitHub Container Registry namespace:

```text
ghcr.io/bdrd-genomics/
```

For example:

```text
ghcr.io/bdrd-genomics/allinone-aio_qc:1.0.1
ghcr.io/bdrd-genomics/allinone-blast:1.0.1
ghcr.io/bdrd-genomics/allinone-fastqc:1.0.1
```

The authoritative process-to-image mapping remains `conf/containers.config`.

## Building Docker images

The repository includes:

```text
containers/docker/build_all.sh
```

for building the collection of component images.

An individual image can also be built directly. For example:

```bash
cd containers/docker

docker build \
    -f Dockerfile.fastqc \
    -t ghcr.io/bdrd-genomics/allinone-fastqc:1.0.1 \
    .
```

## Converting Docker images to Apptainer

Docker/OCI images can be converted to Apptainer SIF files with `apptainer build`.

```bash
apptainer build \
    allinone-fastqc-1.0.1.sif \
    docker://ghcr.io/bdrd-genomics/allinone-fastqc:1.0.1
```

Another example:

```bash
apptainer build \
    allinone-aio_qc-1.0.1.sif \
    docker://ghcr.io/bdrd-genomics/allinone-aio_qc:1.0.1
```

## Offline HPC workflow

For disconnected or restricted HPC systems, a typical deployment pattern is:

1. Build or pull the Docker/OCI images on an internet-connected system.
2. Convert each image to an Apptainer `.sif` file.
3. Transfer the `.sif` files to the HPC environment.
4. Configure the local image directory in the site/container configuration.
5. Run with `-profile slurm,apptainer`.

This avoids pulling external images during an offline pipeline run.

## How container configuration is separated

All-In-One keeps execution mode and image selection separate:

| File | Purpose |
|---|---|
| `conf/profiles.config` | Selects execution/backend behavior such as `local`, `slurm`, `docker`, `apptainer`, or `conda` |
| `conf/containers.config` | Maps pipeline processes to container images |
| `conf/conda.config` | Defines Conda/Mamba software environments |
| `conf/site.config.example` | Template for installation-specific paths and overrides |

This makes it possible to switch between configurations such as:

```text
-profile local,apptainer
-profile slurm,apptainer
-profile local,docker
-profile slurm,conda
```

without editing workflow modules.

See [Configuration](configuration.md) for site and database configuration details.
