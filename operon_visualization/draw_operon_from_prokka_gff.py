from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt


def draw_operon_gff(prokka_gff, gff_without_fasta, start, stop):
    """
    Function create a new GFF file without FASTA record and
    produces schematic plot of DNA region specified with coordinates;
    :param prokka_gff: path to GFF annotation file from Prokka
    :param gff_without_fasta: name of new gff file without FASTA records
    :param start: start coordinate of the desired region
    :param stop: end coordinate of the desired region
    :return: None
    """

    # write new gff file without FASTA record
    with open(prokka_gff, 'r') as gff_with_fasta:
        with open(gff_without_fasta, 'w') as new_gff:
            for line in gff_with_fasta:
                if line.startswith("##FASTA"):
                    break
                new_gff.write(line)

    # draw operon by new gff file
    graphic_record = BiopythonTranslator().translate_record(gff_without_fasta)
    operon = graphic_record.crop((start, stop))
    operon.plot(figure_width=10, elevate_outline_annotations=False)
    plt.show()
