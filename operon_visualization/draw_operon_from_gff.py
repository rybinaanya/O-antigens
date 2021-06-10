from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt


def draw_region_by_coordinates(gff_file, start_, end_, prokka=False, gff_without_fasta=None):
    """
    Function produces schematic plot of DNA region specified with coordinates;
    visualizes only gene feature types
    and labels them by gene or its product.
    :param gff_file: path to GFF annotation file
    :param start_: start coordinate of the desired region
    :param end_: end coordinate of the desired region
    :param prokka: True or False; specifies type of GFF annotation
    :param gff_without_fasta: default None, if Prokka it's used to create new gff
    :return: None
    """
    if prokka:
      with open(gff_file, 'r') as gff_with_fasta:
        with open(gff_without_fasta, 'w') as new_gff:
            for line in gff_with_fasta:
                if line.startswith("##FASTA"):
                    break
                new_gff.write(line)

      # draw operon by new gff file
      graphic_record = BiopythonTranslator().translate_record(gff_without_fasta)
      operon = graphic_record.crop((start_, end_))
      operon.plot(figure_width=10, elevate_outline_annotations=False)
      plt.show()
      
    else:
      biopython_translator = BiopythonTranslator()
      biopython_translator.ignored_features_types = ['CDS']
      biopython_translator.label_fields = ["gene", "product"]
      graphic_record = biopython_translator.translate_record(gff_file)

      graphic_record.crop((start_, end_)).plot(figure_width=10, elevate_outline_annotations=False)
      plt.show()
