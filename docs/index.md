# All-In-One Sequencing Pipeline

**All-In-One** is a modular Nextflow DSL2 workflow for infectious-disease sequencing analysis. It supports **paired short-read**, **long-read**, and **hybrid** sequencing data and brings quality control, host/rRNA removal, assembly, taxonomic classification, viral analysis, genome characterization, and reporting into one configurable pipeline.

[View parameters](parameters.md){ .md-button .md-button--primary }

[Usage guide](usage.md){ .md-button }

[Configuration](configuration.md){ .md-button }

## Major workflow areas

| Workflow area                    | Included capabilities                                                         |
| -------------------------------- | ----------------------------------------------------------------------------- |
| Quality control                  | fastp, FastQC, NanoPlot, MultiQC, read-distribution statistics                |
| Host / rRNA removal              | Host read mapping, BBMap/minimap2, RiboDetector                               |
| Assembly                         | SPAdes, metaSPAdes, plasmidSPAdes, Unicycler, Dragonflye, Myloasm, AutoCycler |
| Assembly validation              | CheckM, CheckM2, CheckV                                                       |
| Taxonomic classification         | Kraken2/Bracken, Sourmash, Mash, MetaPhlAn, GOTTCHA2, Taxpasta                |
| Sequence search / viral analysis | DIAMOND BLASTX, BLAST, MMseqs2, VirusSeeker workflows                         |
| Contig characterization          | Prokka, BUSCO, MLST, AMRFinderPlus, RGI, PhiSpy, MOB-suite, PLASMe            |
| Specialized analysis             | Chimeric-contig detection                                                     |
| Reporting                        | Nextflow trace, timeline, execution report, DAG, and QC summaries             |

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

## Documentation

Use the tabs above to move between the usage guide, complete parameter reference, container documentation, and installation/configuration guidance.

All_In_One supports Docker, Apptainer, Conda, and offline execution profiles. See the [container documentation](containers.md) for information about online container retrieval, local Apptainer SIF files, and air-gapped deployments.

For the source code, issues, and release history, visit the [GitHub repository](https://github.com/BDRD-Genomics/All_In_One).
