from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt


def draw_operon_gff(prokka_gff, gff_without_fasta, start, stop):

    # write new gff file without FASTA records
    with open(prokka_gff, 'r') as gff_with_fasta:
        with open(gff_without_fasta, 'w') as new_gff:
            for line in gff_with_fasta:
                if line.startswith("##FASTA"):
                    break
                new_gff.write(line)

    graphic_record = BiopythonTranslator().translate_record(gff_without_fasta)
    operon = graphic_record.crop((start, stop))
    operon.plot(figure_width=10, elevate_outline_annotations=False)
    plt.show()

# Example
draw_operon_gff("/GCA_001874625.1.gff", "new_gff_file.gff", 4509039, 4521933)
