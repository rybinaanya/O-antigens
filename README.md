# Identification and comparison of somatic antigen structures of symbiotic and pathogenic bacteria from Morganellaceae family

## Introduction

Analysis of the functional organization of the O-antigen operons makes it possible to identify genes specific for each O-serogroup and use them both for molecular typing of strains and for identifying pathways of bacterial evolution.

This project aims to identify and compare candidate O-antigen operons in bacteria of the Morganellaceae family with different levels of virulence.

Goals of the project:
* to find genes of O-antigen biosynthesis in the _Providencia_ and _Xenorhabdus_ species analyzing the literature;
* to analyze the quality of all _Providencia_ assemblies available at NCBI (204 assemblies) and select the best one for each species;
* to create and test the pipeline for identifying and visualizing O-antigen operons;
* to compare candidate O-antigen operons in selected _Providencia_ and _Xenorhabdus_ species;

## Programs used in the study

This project was performed with:
* python v3.7
* QUAST v5.0.2 (Gurevich A. et al., 2013)
* Prokka v1.14 (Seemann T., 2014)
* PGAP (Tatusova T. et al., 2016)
* [Operon-Mapper](https://biocomputo.ibt.unam.mx/operon_mapper/) (Blanca T. et al., 2018)
* [BPROM](http://www.softberry.com/berry.phtml?topic=bprom&group=programs&subgroup=gfindb) (Solovyev V. et al., 2011)
* [FindTerm](http://www.softberry.com/berry.phtml?topic=findterm&group=programs&subgroup=gfindb) (Solovyev V. et al., 2011)
* [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/) (Wheeler DL. et al., 2003)
* MAFFT v7.475 (Katoh K., et al., 2002)
* MEGAX v10.2.4 (Sudhir K. et al., 2018)

Python packages:
* Biopython v1.78 (Cock PA. et al., 2009)
* dna-features-viewer v3.0.3 (Zulkowe V. et al., 2020)
* matplotlib v3.3.4
* numpy v1.20.1
* pandas v1.2.3

## Workflow

__1. Download all _Providencia_ assemblies__

In this project we worked with all _Providencia_ assemblies__ available to 04.03.2021 (204 assemblies)

__2. Check assemblies statistics__

Use the follwoing command:
```bash
# Running QUAST
quast ${directory_with_assemblies}/*.gz
```

or to save output in a certain directory, a script (example below):
```bash
#!/bin/bash

path_downloaded_genomes=/home/rybina/O_antigens/ncbi_dataset/data/
path_out=/home/rybina/O_antigens/quast_Providencia_ncbi/

quast.py $(ls ${path_downloaded_genomes}/*.1/*fna) -o ${path_out} -t 10

```

__3. Select the best assemblies, based on statistics:__
* number of contigs
* Total length (>= 50000 bp)
* N50
* L50

Within the framework of this project, we selected assembly of the best quality for each species (see a [jupyter notebook](https://github.com/rybinaanya/O-antigens/blob/main/select_assemblies.ipynb)).

### Summary of selected assemblies

|   Assembly ID   | Organism Name                            |      Link                                                          |                                           
| :--------------:|:--------------------------------------:| :-------------------------------------------------------------------:|
| GCF_001874625.1 | Providencia rettgeri strain RB151        | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_001874625.1/) |       
| GCF_010669105.1 | Providencia stuartii strain MF1          | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_010669105.1/) |    
| GCF_000259175.1 | Providencia stuartii MRSN 2154           | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000259175.1/) |     
| GCA_016618195.1 | Providencia vermicola strain LLDRA6      | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCA_016618195.1/) |     
| GCF_013702185.1 | Providencia alcalifaciens strain 1701003 | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_013702185.1/) |     
| GCF_900455075.1 | Providencia rustigianii strain NCTC11802 | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_900455075.1/) |    
| GCF_900635875.1 | Providencia rustigianii strain NCTC6933  | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_900635875.1/) |
| GCF_014652175.1 | Providencia thailandensis KCTC 23281     | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_014652175.1/) |  
| GCA_013702025.1 | Providencia rettgeri strain 2055         | [Assambly](https://www.ncbi.nlm.nih.gov/assembly/GCF_013702025.1/) |    

__4. Annotation__

* Prokka
```bash
for file in ${directory_with_assemblies}/*.fna; do prokka --prefix "${file%.fna}"; done
```

or for specifying genus and name of output folder, run the command (you are in the working directory already):
```bash
for f in $(ls); do prokka --outdir  ${f}/${f}_prokka --force --genus Providencia --prefix ${f} ${f}/*.fna; done
```

__5. Identify operon boundaries with Operon-Mapper__

Upload genome FASTA file to the web-server and submit with default options. After job completion, download output compressed folder that includes list of predicted operons with their gene compostion, predicted ORFs coordinates and their DNA sequences and other files.

__6. Search for candidate O-antigen operons__

For PGAP annotation, [Python script](https://github.com/rybinaanya/O-antigens/blob/main/operon_search_pgap.py) was used.
Script for Prokka annotation is deposited [here](https://github.com/rybinaanya/O-antigens/blob/main/operon_search_prokka.py).

__7. Validate boundaries of candidate O-antigen operons with BPROM and FindTerm__

Sequence extraction was performed using [Python script](https://github.com/rybinaanya/O-antigens/blob/main/extract_sequence.py)

__8. Visualize O-antigen operons__

All scripts and examples of figures available [here](https://github.com/rybinaanya/O-antigens/tree/main/operon_visualization).

Example of operon search and operon visualization in _Xenorhabdus bovienii_ str. CS03 could be found [here](https://github.com/rybinaanya/O-antigens/blob/main/Xenorhabdus_bovienii_pgap_example.ipynb). 

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


Sometimes annotation of a gene was undefined (as in the picture above: name for one gene is absent). In that case, we extracted DNA fragment covering the respective gene and its upstream and downstream 300-nucleotide region using [Python script](https://github.com/rybinaanya/O-antigens/blob/main/extract_sequence.py) and clarified annotation of resulting sequence using BLAST and ORFfinder.


__9. Codon-based test of neutrality__

This stage of analysis includes:
* creating multi-FASTA files for each gene from the conserved operons. Using [this script](https://github.com/rybinaanya/O-antigens/blob/main/record_multifasta.py), you can extract gene sequences by gene coordinates from a PROKKA gff file. To get sequences of conserved operon genes from PGAP gff annotation file, use [this script](https://github.com/rybinaanya/O-antigens/blob/main/get_conserved_operon_fasta_from_pgap.py). Resulting multi-FASTA files could be found [here](https://github.com/rybinaanya/O-antigens/tree/main/Z_test) 
* running multiple sequence alignment on gene multi-FASTA files using MAFFT:
```
for f in *.fasta; do mafft ${f} > ${f%%.*}.mafft.fa; done
```
* computing the Z-test statistics of neutral evolution in MEGAX: syn-nonsynonymous substitution type and Nei-Gojobori (Jukes-Cantor) model were applied on 1000 bootstraps with pairwise deletion as gaps/missing data treatment


## Results

1. Out of 204 _Providencia_ assemblies currently available at NCBI, we selected only seven complete genomes and two genomes of contig/scaffold assembly level that demonstrated the best quality. In total, we analyzed six distinct _Providencia_ species.

2. We have designed a unique pipeline for identifying candidate O-antigen operons. The key features of our pipeline are: employing several genome annotation tools, additional prediction and validation of operon boundaries, and manual curation and annotation of unknown genes. Our pipeline allowed us to detect 23 and 12 genes in variable and conserved candidate O-antigen operons, respectively. Among these genes, only seven (_galE_, _wxz_, _wzc_, _wza_, _ugd_, _rmlA_, and _wecA_)  were previously described as O-antigen ones in _Providencia_ and _Xenorhabdus_ species.  

3. Both Providencia and Xenorhabdus species have O-antigen conserved operon involved in the nucleotide (UDP- or dTDP-) sugar biosynthesis (_wecB_, _wecC_, _rffG_, _rfbA_, _rffC_, _rffA_), glucosyl to lipid transfer (_wecA_, _wecF_, _wecG_), and O-antigen processing (_wzzE_,  _wzxE_, _wzyE_). 

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

Almost all _Providencia_ used in the study have _galETKMR_ operon (except for _Providencia alcalifaciens_ strain 1701003) with the same genetic order and gene content.

|    Gene    | P-value |       Z(dN−dS)             | Selection |                
| :---------:|:-------:| :-------------------------:|:---------:|
| galE       | 0.0     |       -8.15                |negative   |
| galT       | 0.0     |       -6.20                |negative   |
| galK       | 0.0     |       -6.33                |negative   |
| galM       | 0.0     |       -8.22                |negative   |
| galR       | 0.0     |       -7.53                |negative   |

4. We did not observe a correlation between lifestyle and O-antigen operon organization as we initially expected. 

For more details on results, check the [slides](https://github.com/rybinaanya/O-antigens/blob/main/O-antigens_presentation.pdf)

## Authors:

* Anna Churkina, Almazov National Medical Research Centre, Saint Petersburg, Russia
* Anna Rybina, Skolkovo Institute of Science and Technology, Moscow, Russia

Supervisors:
* Polina Kuchur, Aleksey Komissarov, 
Applied Genomics Laboratory, SCAMT Institute, ITMO University, 191002 Saint-Petersburg, Russia

## References

1. Bisch, Gaëlle, Sylvie Pagès, John G. McMullen 2nd, S. Patricia Stock, Bernard Duvic, Alain Givaudan, and Sophie Gaudriault. 2015. “Xenorhabdus Bovienii CS03, the Bacterial Symbiont of the Entomopathogenic Nematode Steinernema Weiseri, Is a Non-Virulent Strain against Lepidopteran Insects.” Journal of Invertebrate Pathology 124 (January): 15–22.
2. Herbert, Erin E., and Heidi Goodrich-Blair. 2007. “Friend and Foe: The Two Faces of Xenorhabdus Nematophila.” Nature Reviews. Microbiology 5 (8): 634–46.
3. Ovchinnikova, Olga G., Bin Liu, Dan Guo, Nina A. Kocharova, Magdalena Bialczak-Kokot, Alexander S. Shashkov, Lu Feng, Antoni Rozalski, Lei Wang, and Yuriy A. Knirel. 2012. “Structural, Serological, and Genetic Characterization of the O-Antigen of Providencia Alcalifaciens O40.” FEMS Immunology and Medical Microbiology 66 (3): 382–92.
4. Samuel, Gabrielle, and Peter Reeves. 2003. “Biosynthesis of O-Antigens: Genes and Pathways Involved in Nucleotide Sugar Precursor Synthesis and O-Antigen Assembly.” Carbohydrate Research 338 (23): 2503–19.
5. Wang, Lei, Quan Wang, and Peter R. Reeves. 2010. “The Variation of O Antigens in Gram-Negative Bacteria.” Sub-Cellular Biochemistry 53: 123–52.
6. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler. QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075. doi: 10.1093/bioinformatics/btt086
7. Torsten Seemann, Prokka: rapid prokaryotic genome annotation, Bioinformatics, Volume 30, Issue 14, 15 July 2014, Pages 2068–2069, doi: 10.1093/bioinformatics/btu153
8. Tatusova T, DiCuccio M, Badretdin A, Chetvernin V, Nawrocki EP, Zaslavsky L, Lomsadze A, Pruitt KD, Borodovsky M, Ostell J. NCBI prokaryotic genome annotation pipeline. Nucleic Acids Res. 2016 Aug 19;44(14):6614-24. doi: 10.1093/nar/gkw569
9. Blanca Taboada, Karel Estrada, Ricardo Ciria, Enrique Merino, Operon-mapper: a web server for precise operon identification in bacterial and archaeal genomes, Bioinformatics, Volume 34, Issue 23, 01 December 2018, Pages 4118–4120, doi: 10.1093/bioinformatics/bty496
10. V. Solovyev & A Salamov (2011) Automatic Annotation of Microbial Genomes and Metagenomic Sequences. In Metagenomics and its Applications in Agriculture, Biomedicine and Environmental Studies (Ed. R.W. Li), Nova Science Publishers, p. 61-78
11. Wheeler DL, Church DM, Federhen S, et al. Database resources of the National Center for Biotechnology. Nucleic Acids Res. 2003;31(1):28-33. doi:10.1093/nar/gkg033
12. Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30(14):3059-3066. doi:10.1093/nar/gkf436
13. Sudhir Kumar, Glen Stecher, Michael Li, Christina Knyaz, Koichiro Tamura, MEGA X: Molecular Evolutionary Genetics Analysis across Computing Platforms, Molecular Biology and Evolution, Volume 35, Issue 6, June 2018, Pages 1547–1549, doi: /10.1093/molbev/msy096
14. Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
15. Valentin Zulkower, Susan Rosser, DNA Features Viewer: a sequence annotation formatting and plotting library for Python, Bioinformatics, Volume 36, Issue 15, 1 August 2020, Pages 4350–4352, doi: 10.1093/bioinformatics/btaa213
