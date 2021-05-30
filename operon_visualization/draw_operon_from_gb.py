from matplotlib import pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO


def draw_operon_one_contig_assembly(gb_file, annotation, start, stop):
    """
    Function produces schematic plot of DNA region specified with coordinates
    by any (Prokka or PGAP) annotation of assembly in one contig;
    :param gb_file: path to GB annotation file
    :param annotation: name annotation: Prokka or PGAP
    :param start: start coordinate of the desired region
    :param stop: end coordinate of the desired region
    :return: None
    """

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


def draw_operon_several_contigs_assembly(gb_file, annotation, record_id, start,
                                         stop):
    """
    Function produces schematic plot of DNA region specified with coordinates
    by any (Prokka or PGAP) annotation of assembly in two or more contigs;
    :param gb_file: path to GB annotation file
    :param annotation: name annotation: Prokka or PGAP
    :param record_id: LOCUS from GB file
    :param start: start coordinate of the desired region
    :param stop: end coordinate of the desired region
    :return: None
    """

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
