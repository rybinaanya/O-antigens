from matplotlib import pyplot as plt
from dna_features_viewer import BiopythonTranslator


wec_operon = ['wecA', 'wecB', 'rffG', 'rffH', 'wecD', 'wecE', 'wzxE', 'wecF',
              'wzyE', 'wzyE', 'wecG']

gal_operon = ['galE', 'galT', 'galK', 'galM', 'galR']


def print_operon(prokka_gb, start, stop):
    graphic_record = BiopythonTranslator().translate_record(prokka_gb)
    operon = graphic_record.crop((start, stop))
    operon.plot(figure_width=10, strand_in_label_threshold=7)
    plt.show()


def visualize_conserved_operon(prokka_gff, prokka_gb, operon):

    # find_conserved_operon_boundaries
    conserved_operon_boundaries = []
    with open(prokka_gff) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            if "gene" in line:
                line = line.strip().split()
                info = line[8]
                strand = line[6]
                if operon[0] in info:
                    if strand == "+":
                        conserved_operon_start = int(line[3])
                        conserved_operon_boundaries.append(
                            conserved_operon_start)
                    else:
                        conserved_operon_start = int(line[4])
                        conserved_operon_boundaries.append(
                            conserved_operon_start)
                if operon[-1] in info:
                    if strand == "+":
                        conserved_operon_start = int(line[4])
                        conserved_operon_boundaries.append(
                            conserved_operon_start)
                    else:
                        conserved_operon_start = int(line[3])
                        conserved_operon_boundaries.append(
                            conserved_operon_start)

        print_operon(prokka_gb, conserved_operon_boundaries[0],
                     conserved_operon_boundaries[1])


# Example
visualize_conserved_operon("GCA_900635875.1.gff",
                           "GCA_900635875.1.gbf",
                           gal_operon)
