# Usage

## Clone the pipeline

```bash
git clone https://github.com/BDRD-Genomics/All_In_One.git
cd All_In_One
```

## Input sample sheet

The pipeline reads a CSV sample sheet with these columns:

```text
sample_id,fastq_1,fastq_2,long_read
```

Use `No_Read` when a sequencing input is absent for a sample.

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

## Select the sequencing mode

Set the appropriate mode for the input data:

```bash
--shortreads true
--longreads true
--hybrid true
```

Only enable the mode that matches the sample sheet for the run.

## Example run

A typical Slurm + Apptainer invocation is:

```bash
nextflow run main.nf \
    -profile slurm,apptainer \
    --samplesheet samplesheet.csv \
    --project_id example_project \
    --shortreads true \
    --run_qc true
```

For every available workflow option, see the [Parameters](parameters.md) tab.

## Resume a run

Nextflow can reuse completed tasks from the work directory:

```bash
nextflow run main.nf -resume [other options]
```

## Execution backends

The repository contains configuration for multiple software/execution approaches. Available profiles and exact behavior are defined in `conf/profiles.config` and the other files under `conf/`.

Common profile combinations include:

```bash
-profile slurm,apptainer
```

or, where appropriate:

```bash
-profile local,docker
```

Check the [Configuration](configuration.md) tab before running on a new system.
