from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator


def draw_region_by_coordinates(gff_file, start_, end_):
    """
    Function produces schematic plot of DNA region specified with coordinates;
    visualizes only gene feature types
    and labels them by gene or its product.
    :param gff_file: path to GFF annotation file from PGAP
    :param start_: start coordinate of the desired region
    :param end_: end coordinate of the desired region
    :return: None
    """
    biopython_translator = BiopythonTranslator()
    biopython_translator.ignored_features_types = ['CDS']
    biopython_translator.label_fields = ["gene", "product"]
    graphic_record = biopython_translator.translate_record(gff_file)

    graphic_record.crop((start_, end_)).plot()
