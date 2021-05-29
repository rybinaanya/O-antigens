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
>Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler. QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075. doi: 10.1093/bioinformatics/btt086
* Prokka v1.14
>Torsten Seemann, Prokka: rapid prokaryotic genome annotation, Bioinformatics, Volume 30, Issue 14, 15 July 2014, Pages 2068–2069, doi: 10.1093/bioinformatics/btu153
* PGAP 
>Tatusova T, DiCuccio M, Badretdin A, Chetvernin V, Nawrocki EP, Zaslavsky L, Lomsadze A, Pruitt KD, Borodovsky M, Ostell J. NCBI prokaryotic genome annotation pipeline. Nucleic Acids Res. 2016 Aug 19;44(14):6614-24. doi: 10.1093/nar/gkw569
* [Operon-Mapper](https://biocomputo.ibt.unam.mx/operon_mapper/)
>Blanca Taboada, Karel Estrada, Ricardo Ciria, Enrique Merino, Operon-mapper: a web server for precise operon identification in bacterial and archaeal genomes, Bioinformatics, Volume 34, Issue 23, 01 December 2018, Pages 4118–4120, doi: 10.1093/bioinformatics/bty496
* [BPROM](http://www.softberry.com/berry.phtml?topic=bprom&group=programs&subgroup=gfindb)
>V. Solovyev & A Salamov (2011) Automatic Annotation of Microbial Genomes and Metagenomic Sequences. In Metagenomics and its Applications in Agriculture, Biomedicine and Environmental Studies (Ed. R.W. Li), Nova Science Publishers, p. 61-78
* [FindTerm](http://www.softberry.com/berry.phtml?topic=findterm&group=programs&subgroup=gfindb)
>V. Solovyev, A Salamov (2011) Automatic Annotation of Microbial Genomes and Metagenomic Sequences. In Metagenomics and its Applications in Agriculture, Biomedicine and Environmental Studies (Ed. R.W. Li), Nova Science Publishers, p. 61-78
* MAFFT
>Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30(14):3059-3066. doi:10.1093/nar/gkf436
* MEGAX v10.2.4
>Sudhir Kumar, Glen Stecher, Michael Li, Christina Knyaz, Koichiro Tamura, MEGA X: Molecular Evolutionary Genetics Analysis across Computing Platforms, Molecular Biology and Evolution, Volume 35, Issue 6, June 2018, Pages 1547–1549, doi: /10.1093/molbev/msy096

Python packeges:
* Biopython v1.78
>Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
* dna-features-viewer v3.0.3
>Valentin Zulkower, Susan Rosser, DNA Features Viewer: a sequence annotation formatting and plotting library for Python, Bioinformatics, Volume 36, Issue 15, 1 August 2020, Pages 4350–4352, doi: 10.1093/bioinformatics/btaa213
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

## Results

![Z test wec](/Z_test/Z_test_wec_operon.jpeg)
![Z test gal](/Z_test/Z_test_gal_operon.jpeg)


## Authors:

* Anna Churkina, Almazov National Medical Research Centre, Saint Petersburg, Russia
* Anna Rybina, Skolkovo Institute of Science and Technology, Moscow, Russia

Supervisors:
Polina Kuchur, Aleksey Komissarov
Applied Genomics Laboratory, SCAMT Institute, ITMO University, 191002 Saint-Petersburg, Russia
