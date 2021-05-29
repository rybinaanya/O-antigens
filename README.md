# Identification and comparison of somatic antigen structures of symbiotic and pathogenic bacteria from Morganellaceae family

## Introduction

Analysis of the functional organization of the O-antigen operons makes it possible to identify genes specific for each O-serogroup and use them both for molecular typing of strains and for identifying pathways of bacterial evolution.

The aim of this project is to identify and compare candidate O-antigen operons in bacteria of the Morganellaceae family with different levels of virulence.

Goals of the project:
* to find genes of O-antigen biosynthesis in the Providencia and Xenorhabdus species analyzing the literature;
* to analyze the quality of all Providencia assemblies available at NCBI (204 assemblies) and select the best one for each species;
* to create and test the pipeline for identifying, visualizing O-antigen operons;
* to compare candidate O-antigen operons in selected Providencia and Xenorhabdus species;

## Programs used in the study

This project was performed with:
* python v3.7
* QUAST v5.0.2
* Prokka v1.14
* PGAP 
* [Operon-Mapper](https://biocomputo.ibt.unam.mx/operon_mapper/)
* [BPROM](http://www.softberry.com/berry.phtml?topic=bprom&group=programs&subgroup=gfindb)
* [FindTerm](http://www.softberry.com/berry.phtml?topic=findterm&group=programs&subgroup=gfindb)
* MAFFT
* MEGAX v10.2.4

Python packeges:
* Biopython v1.78
* dna-features-viewer v3.0.3
* matplotlib v3.3.4
* numpy v1.20.1

## Workflow

1. Download all Providencia assemblies
NCBI link: https://www.ncbi.nlm.nih.gov/assembly/?term=Providencia

2. Check assemblies statistics
```bash
# Running QUAST
quast ${directory_with_assemblies}/*.gz
```

3. Select better assemblies, based on statistics:
* number of contigs
* N50
* L50

Within the framework of this project, a jupyter notebook was written, in which the quality of all assemblies was compared and the best ones were selected for further analysis.

4. Identify operon boundaries with Operon-Mapper

5. Validate operon boundaries with BPROM and FindTerm
Sequense extraction was performed using [Python script](https://github.com/rybinaanya/O-antigens/blob/main/extract_sequence.py)

6. Visualize O-antigen operons
The project used two annotations(Prokka and PGAP), which have some differences in the gb file structure. Therefore, to visualize the operons, we created scripts that use the files of each annotation.

All scripts and examples of figures available [here](https://github.com/rybinaanya/O-antigens/tree/main/operon_visualization)

7. Codon-based test of neutrality
