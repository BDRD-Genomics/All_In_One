# All-In-One Sequencing Pipeline

**All-In-One** is a modular Nextflow DSL2 workflow for infectious-disease sequencing analysis. It supports **paired short-read**, **long-read**, and **hybrid** sequencing data and brings quality control, host/rRNA removal, assembly, taxonomic classification, viral analysis, genome characterization, automated reference mapping, and reporting into one configurable pipeline.

[View parameters](parameters.md){ .md-button .md-button--primary }
[Usage guide](usage.md){ .md-button }
[Configuration](configuration.md){ .md-button }

## Major workflow areas

| Workflow area | Included capabilities |
|---|---|
| Quality control | fastp, FastQC, NanoPlot, MultiQC, read-distribution statistics |
| Host / rRNA removal | host read mapping, BBMap/minimap2, RiboDetector |
| Assembly | SPAdes, metaSPAdes, plasmidSPAdes, Unicycler, Dragonflye, Raven, Myloasm, CLC |
| Assembly validation | CheckM / CheckM2 / CheckV and AutoCycler workflows |
| Taxonomic classification | Kraken2/Bracken, Sourmash, Mash, MetaPhlAn, GOTTCHA, Taxpasta |
| Sequence search / viral analysis | DIAMOND BLASTX, MMseqs2, VirusSeeker workflows |
| Contig characterization | Prokka, BUSCO, MLST, AMRFinderPlus, RGI, PhiSpy, MOB-suite, PLASMe |
| Specialized analysis | VF classifier, chimeric-contig detection, automated read mapping |
| Reporting | Nextflow trace, timeline, report, DAG, QC summaries, and exercise reports |

## Repository layout

```text
All_In_One/
├── main.nf
├── nextflow.config
├── params.config
├── conf/
├── workflows/
├── scripts/
├── env/
├── containers/docker/
├── ICTV/
└── docs/
```

## Documentation

Use the tabs above to move between the usage guide, the complete parameter reference, and installation/configuration guidance.

For the source code, issues, and release history, visit the [GitHub repository](https://github.com/BDRD-Genomics/All_In_One).
