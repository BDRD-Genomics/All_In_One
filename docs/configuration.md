# Installation and configuration

## Site configuration

Installation-specific locations are kept outside the tracked pipeline defaults.

```bash
cp conf/site.config.example site.config
```

Edit the copied file for the target system. It controls:

- Container directory
- Database root
- Temporary and Conda cache directories
- Dorado installation and polishing models
- MEGAN command-line tools

Pass the file before the `run` command:

```bash
nextflow -c site.config run main.nf ...
```


## Profiles

Execution and software profiles are composable.

| Profile | Purpose |
|---|---|
| `local` | Execute processes on the current machine. |
| `slurm` | Submit processes through Slurm. |
| `apptainer` | Use configured Apptainer images. |
| `conda` | Create and use repository-defined Conda environments. |
| `docker` | Use Docker containers. |

Typical combinations:

```bash
# Cluster with Apptainer
nextflow -c site.config run main.nf -profile slurm,apptainer -params-file params.yml

# Cluster with Conda
nextflow -c site.config run main.nf -profile slurm,conda -params-file params.yml

# Local testing with Conda
nextflow -c site.config run main.nf -profile local,conda -params-file params.yml
```

Do not combine multiple software backends in one run.

## Databases

All_In_One uses a single database root, `database_dir`. The recommended setup is to place all required external databases beneath this directory using the directory structure defined in `conf/databases.config`.

Set the root in `site.config`:

```groovy
params {
    database_dir = '/data/All_In_One_databases'
}
```

Alternatively, set it with an environment variable:

```bash
export AIO_DATABASE_DIR=/data/All_In_One_databases
```

A typical layout is:

```text
/data/All_In_One_databases/
├── amrfinderplus/
├── bbmap/
├── blastdb/
├── busco/
├── checkm/
├── checkm2/
├── checkv/
├── gottcha/
├── kraken2/
├── krona/
├── mash/
├── metaphlan4/
├── mmseqs/
├── mob_suite/
├── PLASMe/
├── sourmash/
└── taxonomy/
```

Only databases required by enabled workflows need to be installed.

### Overriding individual database locations

The paths in `conf/databases.config` are defaults. If an existing database is stored elsewhere, override only that parameter without reorganizing the rest of the database tree.

For example:

```bash
nextflow run main.nf \
    --database_dir /data/All_In_One_databases \
    --checkm2_db /software/checkm2/uniref100.KO.1.dmnd \
    --diamond_dbdir /shared/blast/nr/nr
```

The same overrides may be placed in `site.config`:

```groovy
params {
    database_dir = '/data/All_In_One_databases'

    checkm2_db    = '/software/checkm2/uniref100.KO.1.dmnd'
    diamond_dbdir = '/shared/blast/nr/nr'
}
```

## Containers

Container filenames are defined in `conf/containers.config` and resolved beneath `container_dir`.


## Building the documentation site

The Markdown files render directly on GitHub. A local website matching the grouped parameter-reference layout can be built with MkDocs:

```bash
python -m pip install -r requirements-docs.txt
mkdocs serve
```

Open the address printed by MkDocs. GitHub Pages deployment is configured in `.github/workflows/docs.yml`.

In the GitHub repository settings, open **Pages** and select **GitHub Actions** as the source. A push to `main` that changes the documentation or parameter schema will then rebuild and publish the site.
