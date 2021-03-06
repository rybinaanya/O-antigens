import pandas as pd
from matplotlib import pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator


O_ag_gene_list = ['yibK', 'wzc', 'wzb', 'wza', 'galE', 'gne', 'wpaD',
                  'ugd', 'wpaC', 'wpaB', 'wzy', 'wpaA', 'wzx', 'qdtB',
                  'qdtA', 'rmlA', 'qdtf', 'cpxA', 'wec', 'rffG', 'rffH']


def find_o_ag_genes_data_from_gff(gff_file, O_ag_gene_list, prokka=False):
    """
    Function gets data (start and end coordinates, strand (+,-)) on candidate O-antigen operon genes from GFF annotation
    using list of O-antigen operon genes defined from literature
    :param gff_file: path to GFF annotation file (PGAP or PROKKA)
    :param O_ag_gene_list: list of O-antigen operon genes defined from literature
    :param prokka: True or False; specifies type of GFF annotation
    :return: 2-D array with 4 columns: start coordinate, end coordinate, strand, gene name
    """
    if prokka:
        with open(gff_file) as gff_with_fasta:
            i = 0
            for line in gff_with_fasta:
                i += 1
                if line.startswith("##FASTA"):
                    start_fasta = i
            end_fasta = i

        ncbi_gff = pd.read_csv(gff_file,
                               sep='\t', comment="#", skipfooter=end_fasta - start_fasta + 1).dropna()
        ncbi_gff.columns = [i for i in range(1, len(ncbi_gff.columns) + 1)]
        ncbi_gff_gene = ncbi_gff.loc[ncbi_gff[9].str.contains('Name=')]
        ncbi_gff_gene[10] = ncbi_gff_gene[9].apply(lambda x: x.split('Name=')[1].split(';')[0])

    else:
        ncbi_gff = pd.read_csv(gff_file, sep='\t', comment="#").dropna()
        ncbi_gff.columns = [i for i in range(1, len(ncbi_gff.columns) + 1)]
        ncbi_gff_gene = ncbi_gff.loc[ncbi_gff[3] == 'gene']
        ncbi_gff_gene.loc[:, 10] = ncbi_gff_gene[9].apply(lambda x: x.split('Name=')[1].split(';')[0])

    return ncbi_gff_gene[
               ncbi_gff_gene[9].apply(lambda x: any([k in x for k in O_ag_gene_list]))
           ].loc[:, [4, 5, 7, 10]].to_numpy()


def find_o_ag_operon_numbers(operons_df, coord_genes_array):
    """
    Function gets candidate O-ag operons (as operon numbers) from Operon-mapper output
    using coordinates array of O-ag genes obtained from GFF annotation
    :param operons_df: pandas DataFrame of Operon-mapper output with the description of predicted operons
    :param coord_genes_array: 2-D array with the following columns: start coordinate, end coordinate, strand, gene name.
    :return: numpy array: operon numbers of candidate O-antigen operons
    """
    return operons_df.dropna().loc[
        (operons_df.dropna().PosLeft.isin(coord_genes_array[:, 0])) &
        (operons_df.postRight.isin(coord_genes_array[:, 1])) |
        (operons_df.dropna().Function.str.contains('O-antigen'))
    ]['Operon'].unique()


def get_operon_boundaries(operons_df, operon_numbers_list):
    """
    Function gets boundary coordinates of the DNA region specified as the list of operons.
    :param operons_df: pandas DataFrame of Operon-mapper output with the description of predicted operons
    :param operon_numbers_list: list of operon numbers (int) predicted by Operon-mapper
    :return: tuple consisting of 2 integers:
    (<start_coordinate_of_the_first_gene_from_the_first_operon>, <end_coordinate_of_the_last_gene_from_the_last_operon>)
    """
    operon_boundaries = operons_df[
                            operons_df.Operon.isin(operon_numbers_list)
                        ].iloc[[1, -1], :].loc[:, ['PosLeft', 'postRight']].to_numpy()
    return tuple(map(int, [operon_boundaries[0][0], operon_boundaries[1][1]]))


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

if __name__ == '__main__':
  # Example
  coord_genes_array = find_o_ag_genes_data_from_gff('GCF_000973125.1_ASM97312v1_genomic.gff', O_ag_gene_list)

  # Parse output of operon-mapper
  operons_df = pd.read_csv("59353/list_of_operons_59353", sep='\t')
  operons_df["Operon"] = operons_df["Operon"].ffill()

  selected_operons = find_o_ag_operon_numbers(operons_df, coord_genes_array)
  
  # Get boundaries coordinates of a pair of operons from 'selected_operons'
  start, end = get_operon_boundaries(operons_df, selected_operons[:2])
  
  # Visualize the first pair of operons from list 'selected_operons'
  draw_region_by_coordinates(gff_file, start, end)
