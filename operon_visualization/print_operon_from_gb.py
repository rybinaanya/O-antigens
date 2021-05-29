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


# Example
# Prokka
draw_operon_one_contig_assembly("GCA_000259175.1.gbf",
                                'Prokka', 2498562, 2511389)

#PGAP
draw_operon_one_contig_assembly("GCF_013702185.1_ASM1370218v1_genomic.gbff",
                                'PGAP', 3819664, 3832542)


def draw_operon_several_contigs_assembly(gb_file, record_id, start,
                                         stop):

    for record in SeqIO.parse(gb_file, 'genbank'):
        if record.id == record_id:
            graphic_translator = BiopythonTranslator()
            # graphic_translator.ignored_features_types = ['gene']  # for PGAP annotation
            graphic_record = graphic_translator.translate_record(record)
            graphic_record.crop((start, stop)).plot(figure_width=10,
                                                    strand_in_label_threshold=7)
            plt.show()


# Example
# Prokka
draw_operon_several_contigs_assembly("GCA_001874625.1.gbf",
                                     "CP017671.1", 1040675, 1063933)

# PGAP
draw_operon_several_contigs_assembly("GCF_010669105.1_ASM1066910v1_genomic.gbff",
                                     "NZ_CP048621.1", 4490816,  4497410)
