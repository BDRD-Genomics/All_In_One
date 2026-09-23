# Configuration

All-In-One separates portable pipeline logic from installation-specific configuration. The main configuration files are under `conf/`.

## Configuration files

| File | Purpose |
|---|---|
| `conf/profiles.config` | Execution and software profiles |
| `conf/resources.config` | Process resource definitions |
| `conf/containers.config` | Container image definitions |
| `conf/conda.config` | Conda environment configuration |
| `conf/databases.config` | Database configuration |
| `conf/databases.config.internal` | Internal database configuration |
| `conf/site.config.example` | Template for installation-specific settings |

## Site configuration

Create a local site configuration from the template:

```bash
cp conf/site.config.example site.config
```

Edit `site.config` for the local environment. Installation-specific paths such as container locations, database roots, temporary storage, and external software paths belong here rather than in workflow modules.

Do not commit a site configuration that contains credentials or sensitive local infrastructure information.

## Containers

Container definitions and build recipes are under:

```text
containers/docker/
```

The repository includes individual Dockerfiles and Conda environment definitions for pipeline components. Container mappings used by Nextflow are defined in `conf/containers.config`.

## Conda

Conda environments are stored under `env/`, with additional image-build environment files under `containers/docker/`.

The active backend is selected through the configured Nextflow profile rather than by editing individual workflow modules.

## Databases

Database locations are configured in:

```text
conf/databases.config
conf/databases.config.internal
```

Only databases needed by enabled workflows need to be available for a given run.

## Slurm and resources

Cluster behavior is controlled through `conf/profiles.config` and `conf/resources.config`. Site-specific queue, account, storage, or scheduler settings should remain outside reusable workflow code wherever possible.

## Parameters

Pipeline behavior is controlled through command-line parameters. See the complete [Parameters](parameters.md) reference.
