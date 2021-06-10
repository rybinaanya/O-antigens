import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Specify arguments
gff_file = sys.argv[1]
output_folder = sys.argv[2]
conserved_operon = sys.argv[3]

conserved_genes = ['wecA', 'wzzE', 'wecB', 'wecC', 'rffG', 'rffH', 
                   'wecD', 'wecE', 'wzxE', 'wecF', 'wzyE', 'wecG', 
                   'galE', 'galT', 'galK', 'galM', 'galR']


dict_conserved_genes = {
    'wecA': 'wecA',
    'wzzE': 'wzzE',
    'wecB': 'wecB',
    'wecC': 'wecC',
    'rffG': 'rffG',
    'rfbA': 'rffH',
    'rffC': 'wecD',
    'rffA': 'wecE',
    'wzxE': 'wzxE',
    'wzyE': 'wzyE',
    'rffM': 'wecG',
    'galR': 'galR',
    'galE': 'galE',
    'galK': 'galK',
    'galM': 'galM',
    'galT': 'galT'
}


def get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, i, gene_names,
                                 operon_strand, dict_conserved_genes):
    """
    Function constructs SeqRecord object of a gene sequence based on given data
    :param fna: path to genome FASTA file
    :param gff: path to GFF annotation file (PGAP)
    :param arr_id_coord_strand: numpy array with 4 columns:
    - id of record sequence in genome FASTA file on which a gene is located;
    - gene start coordinate;
    - gene end coordinate;
    - DNA strand on which a gene is located
    :param gene_names: list of conserved operon genes
    :param operon_strand: '+' or '-' strand on which the operons is located
    (assuming all operon genes are on the same DNA strand)
    :return: SeqRecord object:
    - a gene sequence
    - id consists of the gene name and assembly id;
    - name consists of the gene name and assembly id of organism from which a gene sequence was extracted;
    - description includes name of organism from which a gene sequence was extracted
    """
    fixed_gene_names = [dict_conserved_genes[i]
                        if i in dict_conserved_genes else 'wecF'
                        for i in gene_names]

    for record in SeqIO.parse(fna, 'fasta'):
        if record.id == arr_id_coord_strand[0][0]:
            # for i in range(arr_id_coord_strand.shape[0]):
            if operon_strand == '-':
                return SeqRecord(record.seq[arr_id_coord_strand[i][1] - 1:
                                            arr_id_coord_strand[i][
                                                2]].reverse_complement(),
                                 #id=gene_names[i] + ' gene|' + '_'.join(
                                  #   ['GCF', gff.split('_')[1].split('/')[0]]),
                                 id=fixed_gene_names[i] + ' gene|' + '_'.join(
                                     ['GCF', gff.split('_')[1].split('/')[0]]),
                                 #name=gene_names[i] + ' gene|' + '_'.join(
                                  #   ['GCF', gff.split('_')[1].split('/')[0]]),
                                 name=fixed_gene_names[i] + ' gene|' + '_'.join(
                                     ['GCF', gff.split('_')[1].split('/')[0]]),
                                 description=' '.join(
                                     record.description.split(' ')[1:5]))
            else:
                return SeqRecord(record.seq[arr_id_coord_strand[i][1] - 1:
                                            arr_id_coord_strand[i][2]],
                                 id=fixed_gene_names[i] + ' gene|' + '_'.join(
                                     ['GCF', gff.split('_')[1].split('/')[0]]),
                                 #name=gene_names[i] + ' gene|' + '_'.join(
                                  #   ['GCF', gff.split('_')[1].split('/')[0]]),
                                 name=fixed_gene_names[i] + ' gene|' + '_'.join(
                                     ['GCF', gff.split('_')[1].split('/')[0]]),
                                 description=' '.join(
                                     record.description.split(' ')[1:5]))



def get_conserved_genes_seq(gff, conserved_operon):
    """
    Functions searches for conserved genes from specified GFF annotation file
    and for each gene in operon (if found), gets sequence, id, description and creates a SeqRecord object
    :param gff: GFF PGAP annotation files
    :param conserved_operon: the name of conserved operon whose genes will be recorded in the fasta file
    :return: list of SeqRecord objects constructed for each gene in the operon
    """
    global edge_operon_gene, operon_length_, last

    if conserved_operon == "gal":
        edge_operon_gene = "galR"
        operon_length_ = 5
        last = True

    if conserved_operon == "wec":
        edge_operon_gene = "wecA"
        operon_length_ = 12
        last = False
    lst_new_records = []
    fna = gff.split('_genomic')[0] + '_genomic.fna'
    # load GFF annotation file as pandas DataFrame
    df = pd.read_csv(gff,
                     sep='\t', comment="#").dropna()
    df.columns = [i for i in range(1, len(df.columns) + 1)]
    # extract data from GFF annotation corresponding to gene or pseudogene feature types only
    df_gene = df[(df[3] == 'gene') | (df[3] == 'pseudogene')]
    df_gene.index = list(range(df_gene.shape[0]))
    # get index of the first gene from the conserved operon
    idx = df_gene.index[df_gene[9].str.contains(edge_operon_gene)][0]
    # get strand of the first gene from the conserved operon
    # (genes from this operon are located on the same strand)
    strand = df_gene.loc[idx, 7]
    if last:
        if strand == '-':
            # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
            arr_id_coord_strand = df_gene.loc[idx:(idx + (operon_length_ - 1)),
                                  [1, 4, 5, 7]].to_numpy()
            # get names of conserved operon genes from PGAP GFF annotation
            gene_names = df_gene.loc[idx:(idx + (operon_length_ - 1)), 9].apply(
                lambda x:
                x.split('Name=')[1].split(';')[0]).to_list()

            for i in range(arr_id_coord_strand.shape[0]):
                new_record = get_seqrecord_conserved_gene(fna, gff,
                                                          arr_id_coord_strand,
                                                          i, gene_names, '-', dict_conserved_genes)
                lst_new_records.append(new_record)
        else:
            # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
            arr_id_coord_strand = df_gene.loc[(idx - (operon_length_ - 1)):idx,
                                  [1, 4, 5, 7]].to_numpy()
            # get names of conserved operon genes from PGAP GFF annotation
            gene_names = df_gene.loc[(idx - (operon_length_ - 1)):idx, 9].apply(
                lambda x:
                x.split('Name=')[1].split(';')[0]).to_list()

            for i in range(arr_id_coord_strand.shape[0]):
                new_record = get_seqrecord_conserved_gene(fna, gff,
                                                          arr_id_coord_strand,
                                                          i, gene_names, '+', dict_conserved_genes)
                lst_new_records.append(new_record)
    else:
        if strand == '-':
            # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
            arr_id_coord_strand = df_gene.loc[(idx - (operon_length_ - 1)):idx,
                                  [1, 4, 5, 7]].to_numpy()
            # get names of conserved operon genes from PGAP GFF annotation
            gene_names = df_gene.loc[(idx - (operon_length_ - 1)):idx, 9].apply(
                lambda x: x.split('Name=')[1].split(';')[0]).to_list()

            # new_record = get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, '-')
            for i in range(arr_id_coord_strand.shape[0]):
                new_record = get_seqrecord_conserved_gene(fna, gff,
                                                          arr_id_coord_strand,
                                                          i, gene_names, '-', dict_conserved_genes)
                lst_new_records.append(new_record)
        else:
            # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
            arr_id_coord_strand = df_gene.loc[idx:(idx + (operon_length_ - 1)),
                                  [1, 4, 5, 7]].to_numpy()
            # get names of conserved operon genes from PGAP GFF annotation
            gene_names = df_gene.loc[idx:(idx + (operon_length_ - 1)), 9].apply(
                lambda x:
                x.split('Name=')[1].
                split(';')[0]).to_list()
            for i in range(arr_id_coord_strand.shape[0]):
                new_record = get_seqrecord_conserved_gene(fna, gff,
                                                          arr_id_coord_strand,
                                                          i, gene_names, '+', dict_conserved_genes)
                lst_new_records.append(new_record)
    return lst_new_records


def write_gene_fasta_files(conserved_genes, lst_new_records, output_folder):
    """
    Write multi-FASTA files for each gene from the conserved operon
    :param conserved_genes: list of conserved operon genes
    :param lst_new_records: list of SeqRecord objects of conserved operon genes
    :param output_folder: path to folder where resulting FASTA files would be written to
    :return: None
    """
    for gene in conserved_genes:
        gene_records = []
        with open(f"{output_folder}/{gene}.fasta", 'a') as fout:
            for record in lst_new_records:
                if gene in record.id:
                    gene_records.append(record)
            SeqIO.write(gene_records, fout, "fasta")


if __name__ == '__main__':
    lst_records = get_conserved_genes_seq(gff_file, conserved_operon)
    write_gene_fasta_files(conserved_genes, lst_records, output_folder)
