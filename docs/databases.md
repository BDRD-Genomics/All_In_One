# Database setup

All-In-One uses a single configurable database root.

```bash
export AIO_DATABASE_DIR=/path/to/databases
```

The same location can be supplied to the pipeline as `--database_dir`.

Database provisioning is intentionally separate from Nextflow execution. Large downloads such as NCBI `nt` and `nr` should not begin unexpectedly inside a Slurm job.

## Automated installer

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

Install only selected databases:

```bash
./scripts/download_databases.sh \
    --db-root "$AIO_DATABASE_DIR" \
    --db checkm2,checkv,taxdb
```

Install the full set including NCBI `nt` and `nr`:

```bash
./scripts/download_databases.sh \
    --db-root "$AIO_DATABASE_DIR" \
    --preset full
```

`nt` and `nr` are intentionally opt-in because of their size.

## Presets

| Preset | Databases |
|---|---|
| `minimal` | `taxdb`, `checkm2`, `checkv` |
| `standard` | minimal + `checkm`, `busco`, `kraken2`, `metaphlan4`, `amrfinderplus` |
| `full` | standard + `nt`, `nr` |

## Database sources and update mechanisms

| Database | Automated method | Official/upstream source |
|---|---|---|
| NCBI `nt` | `update_blastdb.pl --decompress nt` | NCBI BLAST database archive |
| NCBI `nr` | `update_blastdb.pl --decompress nr` | NCBI BLAST database archive |
| NCBI `taxdb` | `update_blastdb.pl --decompress taxdb` | NCBI BLAST database archive |
| CheckM v1 | download and extract official CheckM data archive | CheckM/UQ data archive |
| CheckM2 | `checkm2 database --download --path ...` | CheckM2 |
| CheckV | `checkv download_database ...` | CheckV/NERSC |
| BUSCO | `busco --download prokaryota --download_path ...` | BUSCO |
| Kraken2 | `kraken2-build --standard` | Kraken2 / NCBI |
| MetaPhlAn | `metaphlan --install` | bioBakery MetaPhlAn |
| AMRFinderPlus | `amrfinder -u` / `amrfinder -U` | NCBI AMRFinderPlus |

## Expected layout

```text
$AIO_DATABASE_DIR/
├── blast/
│   ├── taxdb/
│   ├── nt/
│   └── nr/
├── checkm/
├── checkm2/
├── checkv/
├── busco/
├── kraken2/
│   └── standard/
├── metaphlan4/
└── amrfinderplus/
```

## Verify installation

```bash
./scripts/verify_databases.sh --db-root "$AIO_DATABASE_DIR"
```

A missing `nt` or `nr` entry is expected when using `minimal` or `standard`.

## Offline / air-gapped systems

Download databases on an internet-connected host:

```bash
export AIO_DATABASE_DIR=/data/All_In_One_databases

./scripts/download_databases.sh \
    --db-root "$AIO_DATABASE_DIR" \
    --preset standard
```

Archive the completed database tree:

```bash
tar -I 'pigz -1' \
    -cf All_In_One_databases.tar.gz \
    -C "$(dirname "$AIO_DATABASE_DIR")" \
    "$(basename "$AIO_DATABASE_DIR")"
```

Transfer the archive to the offline system, extract it, and point `AIO_DATABASE_DIR` / `--database_dir` at the extracted directory.

For BUSCO offline execution, use locally downloaded lineage datasets and enable BUSCO's offline behavior in the workflow configuration.

## Updating

The installer can be rerun. Upstream tools such as NCBI `update_blastdb.pl`, CheckM2, BUSCO, MetaPhlAn, and AMRFinderPlus handle their own update logic.

For AMRFinderPlus, database compatibility is tied to the AMRFinderPlus software version. After upgrading AMRFinderPlus, rerun the database updater.

For large production installations, record the database version/date used for each analysis and avoid silently replacing databases underneath an active project.

## Manual installation

If automated installation is not appropriate, use the source column in `databases/manifest.tsv` to locate the upstream project or database archive and install into the same directory layout shown above.

The pipeline should consume databases from `--database_dir`; it should not download databases during analysis jobs.
