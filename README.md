# Identification and comparison of somatic antigen structures of symbiotic and pathogenic bacteria from Morganellaceae family

## Introduction

Analysis of the functional organization of the O-antigen operons makes it possible to identify genes specific for each O-serogroup and use them both for molecular typing of strains and for identifying pathways of bacterial evolution.

This project aims to identify and compare candidate O-antigen operons in bacteria of the Morganellaceae family with different levels of virulence.

Goals of the project:
* to find genes of O-antigen biosynthesis in the Providencia and Xenorhabdus species analyzing the literature;
* to analyze the quality of all Providencia assemblies available at NCBI (204 assemblies) and select the best one for each species;
* to create and test the pipeline for identifying and visualizing O-antigen operons;
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
* [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/)
>Wheeler DL, Church DM, Federhen S, et al. Database resources of the National Center for Biotechnology. Nucleic Acids Res. 2003;31(1):28-33. doi:10.1093/nar/gkg033
* MAFFT v7.475
>Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30(14):3059-3066. doi:10.1093/nar/gkf436
* MEGAX v10.2.4
>Sudhir Kumar, Glen Stecher, Michael Li, Christina Knyaz, Koichiro Tamura, MEGA X: Molecular Evolutionary Genetics Analysis across Computing Platforms, Molecular Biology and Evolution, Volume 35, Issue 6, June 2018, Pages 1547–1549, doi: /10.1093/molbev/msy096

Python packages:
* Biopython v1.78
>Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
* dna-features-viewer v3.0.3
>Valentin Zulkower, Susan Rosser, DNA Features Viewer: a sequence annotation formatting and plotting library for Python, Bioinformatics, Volume 36, Issue 15, 1 August 2020, Pages 4350–4352, doi: 10.1093/bioinformatics/btaa213
* matplotlib v3.3.4
* numpy v1.20.1
* pandas v1.2.3

## Workflow

1. Download all Providencia assemblies
NCBI link: https://www.ncbi.nlm.nih.gov/assembly/?term=Providencia

2. Check assemblies statistics
using the follwoing command:
```bash
# Running QUAST
quast ${directory_with_assemblies}/*.gz
```

or a script (example below) to save output in a certain directory:
```
#!/bin/bash

path_downloaded_genomes=/home/rybina/O_antigens/ncbi_dataset/data/
path_out=/home/rybina/O_antigens/quast_Providencia_ncbi/

quast.py $(ls ${path_downloaded_genomes}/*.1/*fna) -o ${path_out} -t 10

```

3. Select better assemblies, based on statistics:
* number of contigs
* Total length (>= 50000 bp)
* N50
* L50

Within the framework of this project, a jupyter notebook was written, in which the quality of all assemblies was compared and the best ones were selected for further analysis.

4. Annotation

* Prokka
```bash
for file in ${directory_with_assemblies}/*.fna; do prokka --prefix "${file%.fna}"; done
```

5. Identify operon boundaries with Operon-Mapper

6. Validate operon boundaries with BPROM and FindTerm

Sequense extraction was performed using [Python script](https://github.com/rybinaanya/O-antigens/blob/main/extract_sequence.py)

7. Visualize O-antigen operons

All scripts and examples of figures available [here](https://github.com/rybinaanya/O-antigens/tree/main/operon_visualization)

8. Codon-based test of neutrality

This stage of analysis includes:
* creation of files for each of the genes of the conserved operons. Using [this script](https://github.com/rybinaanya/O-antigens/blob/main/record_multifasta.py), you can extract gene sequences by gene coordinates from a gff file for further use;
* gene alignment
* computation of Z-test statistics of neutral evolution in MEGAX with parameters: Nei-Gojobori (Jukes-Cantor) model, 1000 bootstraps

## Results

1. Out of 204 Providencia assemblies currently available at NCBI, we selected only seven complete genomes and two genomes of contig/scaffold assembly level that demonstrated the best quality. In total, we analyzed six distinct Providencia species.

2. We have designed a unique pipeline for identifying candidate O-antigen operons. The key features of our pipeline are: employing several genome annotation tools, additional prediction and validation of operon boundaries, and manual curation and annotation of unknown genes. Our pipeline allowed us to detect about 20 O-antigen genes that were not previously described in Providencia and Xenorhabdus species. 

An example of visualization an operon obtained using a custom script:

```python 
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator


def draw_region_by_coordinates(gff_file, start_, end_):
    biopython_translator = BiopythonTranslator()
    biopython_translator.ignored_features_types = ['CDS']
    biopython_translator.label_fields = ["gene", "product"]
    graphic_record = biopython_translator.translate_record(gff_file)

    graphic_record.crop((start_, end_)).plot()

draw_region_by_coordinates("GCA_016618195.1_ASM1661819v1_genomic.gff", 4103300, 4116180)
```

![Z test wec](/operon_visualization/img/PGAP_gff.png)


3. Both Providencia and Xenorhabdus species have O-antigen conserved operon involved in the nucleotide (UDP- or dTDP-) sugar biosynthesis (wecB, wecC, rffG, rfbA, rffC, rffA), glucosyl to lipid transfer (wecA, wecF, wecG), and O-antigen processing (wzzE,  wzxE, wxyE).

|    Gene    | P-value |       Z(dN−dS)             | Selection |                                           
| :---------:|:-------:| :-------------------------:|:---------:|
| wecA       | 0.0     |       -8.95                |negative   |
| wzzE       | 0.0     |       -9.25                |negative   |
| wecB       | 0.0     |       -9.35                |negative   |
| wecC       | 0.0     |       -10.32               |negative   |
| rffG       | 0.0     |       -10.47               |negative   |
| rffH       | 0.0     |       -6.84                |negative   |
| wecD       | 0.0     |       -6.47                |negative   |
| wecE       | 0.0     |       -9.71                |negative   |
| wzxE       | 0.0     |       -10.59               |negative   |
| wecF       | 0.0     |       -10.18               |negative   |
| wzyE       | 0.0     |       -5.88                |negative   |
| wecG       | 0.0     |       -6.83                |negative   |

|    Gene    | P-value |       Z(dN−dS)             | Selection |                
| :---------:|:-------:| :-------------------------:|:---------:|
| galE       | 0.0     |       -8.15                |negative   |
| galT       | 0.0     |       -6.20                |negative   |
| galK       | 0.0     |       -6.33                |negative   |
| galM       | 0.0     |       -8.22                |negative   |
| galR       | 0.0     |       -7.53                |negative   |

4. We did not observe a correlation between lifestyle and O-antigen operon organization as we initially expected. 

## Authors:

* Anna Churkina, Almazov National Medical Research Centre, Saint Petersburg, Russia
* Anna Rybina, Skolkovo Institute of Science and Technology, Moscow, Russia

Supervisors:
Polina Kuchur, Aleksey Komissarov
Applied Genomics Laboratory, SCAMT Institute, ITMO University, 191002 Saint-Petersburg, Russia
