from matplotlib import pyplot as plt
from dna_features_viewer import BiopythonTranslator


wec_operon = ['wecA', 'wzzE', 'wecB', 'wecC', 'rffG', 'rffH', 
              'wecD', 'wecE', 'wzxE', 'wecF', 'wzyE', 'wecG']

gal_operon = ['galE', 'galT', 'galK', 'galM', 'galR']


def print_operon(prokka_gb, start, stop):
    graphic_record = BiopythonTranslator().translate_record(prokka_gb)
    operon = graphic_record.crop((start, stop))
    operon.plot(figure_width=10, strand_in_label_threshold=7)
    plt.show()


def visualize_conserved_operon(prokka_gff, prokka_gb, operon):
    """
    Function produces schematic plot of DNA region with boundaries determined
    by searching for the coordinates of the outermost genes of the operon;
    :param prokka_gff: path to GFF annotation file from Prokka
    :param prokka_gb: path to GB annotation file from Prokka
    :param operon: name of conserved operon: gal or wec
    :return: None
    """
    # find_conserved_operon_boundaries
    conserved_operon_boundaries = []
    with open(prokka_gff) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            if "gene" in line:
                line = line.strip().split()
                info = line[8]
                start = line[3]
                stop = line[4]
                strand = line[6]

                # find_coordinate_of_first_gene_in_operon
                if operon[0] in info:
                    if strand == "+":
                        # if gene in + strand, first coordinate of gene is conserved operon start
                        conserved_operon_start = int(start)
                        conserved_operon_boundaries.append(
                            conserved_operon_start)
                    else:
                        # else second coordinate of gene is conserved operon start
                        conserved_operon_start = int(stop)
                        conserved_operon_boundaries.append(
                            conserved_operon_start)

                # find_coordinate_of_first_gene_in_operon
                if operon[-1] in info:
                    if strand == "+":
                        # if gene in + strand, second coordinate of gene is conserved operon end
                        conserved_operon_stop = int(stop)
                        conserved_operon_boundaries.append(
                            conserved_operon_stop)
                    else:
                        # else second coordinate of gene is conserved operon end
                        conserved_operon_stop = int(start)
                        conserved_operon_boundaries.append(
                            conserved_operon_stop)

        print_operon(prokka_gb, conserved_operon_boundaries[0]-1,
                     conserved_operon_boundaries[1])


# Example
visualize_conserved_operon("GCA_016618195.1.gff",
                           "GCA_016618195.1.gbf",
                           wec_operon)
