If you don't use the [automatic method](https://github.com/rybinaanya/O-antigens/blob/main/search_visualize_operons.py) for visualizing candidate O-antigen operons, you can use scripts from this folder, where the arguments are passed manually.

Prokka and PGAP annotations have some differences in the gb and gff file structure. Therefore, to visualize the operons, we created scripts that use the files of each annotation.

An examples of visualization an operon using a custom script:
* [by gb files from Prokka and PGAP annotations](https://github.com/rybinaanya/O-antigens/blob/main/operon_visualization/draw_operon_from_gb.py) with one contig assemblies 

```python 
from matplotlib import pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO


def draw_operon_one_contig_assembly(gb_file, annotation, start, stop):
    if annotation == "PGAP":
        graphic_translator = BiopythonTranslator()
        graphic_translator.ignored_features_types = ['CDS']
        graphic_record = graphic_translator.translate_record(gb_file)
        graphic_record.crop((start, stop)).plot(figure_width=10,
                                                strand_in_label_threshold=7)
        plt.show()

    if annotation == "Prokka":
        graphic_record = BiopythonTranslator().translate_record(gb_file)
        operon = graphic_record.crop((start, stop))
        operon.plot(figure_width=10, strand_in_label_threshold=7)
        plt.show()


# Prokka
draw_operon_one_contig_assembly("GCA_900635875.1.gbf",
                                'Prokka', 471058, 476657)
# #PGAP
draw_operon_one_contig_assembly("GCA_900635875.1_35661_G01_genomic.gbff",
                                'PGAP', 471058, 476657)
```
Prokka:
![prokka](/operon_visualization/img/Prokka_gb_one_contig.png)

PGAP:
![pgap](/operon_visualization/img/PGAP_gb_one_contig.png)

* [by gb files from Prokka and PGAP annotations](https://github.com/rybinaanya/O-antigens/blob/main/operon_visualization/draw_operon_from_gb.py) with two or more contigs assemblies 

```python
from matplotlib import pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO


def draw_operon_several_contigs_assembly(gb_file, annotation, record_id, start,
                                         stop):
    for record in SeqIO.parse(gb_file, 'genbank'):
        if record.id == record_id:
            if annotation == "PGAP":
                graphic_translator = BiopythonTranslator()
                graphic_translator.ignored_features_types = ['gene']
                graphic_record = graphic_translator.translate_record(record)
                graphic_record.crop((start, stop)).plot(figure_width=10,
                                                        strand_in_label_threshold=7)
                plt.show()

            if annotation == "Prokka":
                graphic_translator = BiopythonTranslator()
                graphic_record = graphic_translator.translate_record(record)
                graphic_record.crop((start, stop)).plot(figure_width=10,
                                                        strand_in_label_threshold=7)
                plt.show()


# Prokka
draw_operon_several_contigs_assembly("GCA_001874625.1.gbf", "Prokka",
                                     "CP017671.1", 4509039, 4521933)
# PGAP
draw_operon_several_contigs_assembly("GCF_010669105.1_ASM1066910v1_genomic.gbff",
                                     "PGAP", "NZ_CP048621.1", 4272652,  4285534)
```
Prokka:
![prokka](/operon_visualization/img/Prokka_gb_two_contigs.png)

PGAP:
![pgap](/operon_visualization/img/PGAP_gb_two_contigs.png)

* [by gff file](https://github.com/rybinaanya/O-antigens/blob/main/operon_visualization/draw_operon_from_gff.py)
   
PGAP:
![pgap_gff](/operon_visualization/img/PGAP_gff.png)

PROKKA:
![prokka_gff](/operon_visualization/img/prokka_gff.png)

* [conserved operon from Prokka annotation](https://github.com/rybinaanya/O-antigens/blob/main/operon_visualization/draw_conserved_operon_prokka.py)

![conserved operon from Prokka](/operon_visualization/img/Prokka_conserved_operon.png)
