from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, operon_strand):
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
    for record in SeqIO.parse(fna, 'fasta'):
        if record.id == arr_id_coord_strand[0][0]:
            for i in range(arr_id_coord_strand.shape[0]):
                if operon_strand == '-':
                     return SeqRecord(record.seq[arr_id_coord_strand[i][1] - 1:
                                                 arr_id_coord_strand[i][2]].reverse_complement(),
                                      id=gene_names[i] + ' gene|' + gff.split('/')[1],
                                      name=gene_names[i] + ' gene|' + gff.split('/')[1],
                                      description=' '.join(record.description.split(' ')[1:5]))
                else:
                    return SeqRecord(record.seq[arr_id_coord_strand[i][1] - 1:arr_id_coord_strand[i][2]],
                                     id=gene_names[i] + ' gene|' + gff.split('/')[1],
                                     name=gene_names[i] + ' gene|' + gff.split('/')[1],
                                     description=' '.join(record.description.split(' ')[1:5]))


def get_conserved_genes_seq(lst_gff, edge_operon_gene, operon_length_, last=False):
    """
    Functions searches for conserved genes from the list in specified GFF annotation files
    and for each gene in operon (if found), gets sequence, id, description and creates a SeqRecord object
    :param lst_gff: a list of paths to GFF PGAP annotation files
    :param edge_operon_gene: the first or the last gene of the operon; if the last one then flag 'last' should be True
    :param last: flag; default: False; True if edge_operon_gene is the last gene of the operon
    :param operon_length_: number of genes in the operon
    :return: list of SeqRecord objects constructed for each gene in the operon
    """
    lst_new_records = []
    for gff in lst_gff:
        fna = gff[:-5] + '.fna'
        # load GFF annotation file as pandas DataFrame
        df = pd.read_csv(gff,
                         sep='\t', comment="#").dropna()
        df.columns = [i for i in range(1, len(df.columns) + 1)]
        # extract data from GFF annotation corresponding to gene or pseudogene feature types only
        df_gene = df[(df[3] == 'gene') | (df[3] == 'pseudogene')]
        df_gene.index = range(df_gene.shape[0])
        # get index of the first gene from the conserved operon
        idx = df_gene[df_gene[9].str.contains(edge_operon_gene)].index[0]
        # get strand of the first gene from the conserved operon
        # (genes from this operon are located on the same strand)
        strand = df_gene.loc[idx, 7]
        if last:
            if strand == '-':
                # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
                arr_id_coord_strand = df_gene.loc[idx:(idx + (operon_length_-1)), [1, 4, 5, 7]].to_numpy()
                # get names of conserved operon genes from PGAP GFF annotation
                gene_names = df_gene.loc[idx:(idx + (operon_length_-1)), 9].apply(lambda x:
                                                                 x.split('gene=')[1].split(';')[0]).to_list()
                new_record = get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, '-')
                lst_new_records.append(new_record)
            else:
                # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
                arr_id_coord_strand = df_gene.loc[(idx - (operon_length_-1)):idx, [1, 4, 5, 7]].to_numpy()
                # get names of conserved operon genes from PGAP GFF annotation
                gene_names = df_gene.loc[(idx - (operon_length_-1)):idx, 9].apply(lambda x:
                                                                 x.split('gene=')[1].split(';')[0]).to_list()
                new_record = get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, '+')
                lst_new_records.append(new_record)
        else:
            if strand == '-':
                # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
                arr_id_coord_strand = df_gene.loc[(idx - (operon_length_-1)):idx, [1, 4, 5, 7]].to_numpy()
                # get names of conserved operon genes from PGAP GFF annotation
                gene_names = df_gene.loc[(idx - (operon_length_-1)):idx, 9].apply(lambda x:
                                                                 x.split('gene=')[1].split(';')[0]).to_list()
                new_record = get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, '-')
                lst_new_records.append(new_record)
            else:
                # get id, coordinates, and strand of all conserved operon genes from PGAP GFF annotation
                arr_id_coord_strand = df_gene.loc[idx:(idx + (operon_length_-1)), [1, 4, 5, 7]].to_numpy()
                # get names of conserved operon genes from PGAP GFF annotation
                gene_names = df_gene.loc[idx:(idx + (operon_length_-1)), 9].apply(lambda x:
                                                                                  x.split('gene=')[1].
                                                                                  split(';')[0]).to_list()
                new_record = get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, '+')
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
        with open(f"{output_folder}/{gene}.fasta", 'w') as fout:
            for record in lst_new_records:
                if gene in record.id:
                    gene_records.append(record)
            SeqIO.write(gene_records, fout, "fasta")


# Example of usage (getting sequences of the galETKMR operon genes)
lst_gal_gff = ['./GCF_010669105/GCF_010669105.1_ASM1066910v1_genomic2.gff',
               './GCF_001874625/GCF_001874625.1_ASM187462v1_genomic2.gff',
               './GCF_010669105/GCF_010669105.1_ASM1066910v1_genomic2.gff',
               './GCF_013702025/GCF_013702025.1_ASM1370202v1_genomic2.gff',
               './GCA_016618195/GCA_016618195.1_ASM1661819v1_genomic2.gff',
               './GCF_014652175/GCF_014652175.1_ASM1465217v1_genomic2.gff',
               './GCF_900455075/GCF_900455075.1_28869_H01_genomic2.gff',
               './GCF_900635875/GCF_900635875.1_35661_G01_genomic2.gff']

gal_genes = ['galE', 'galT', 'galK', 'galM', 'galR']
lst_new_records = get_conserved_genes_seq(lst_gal_gff, edge_operon_gene='galR', operon_length_=5, last=True)
write_gene_fasta_files(gal_genes, lst_new_records, 'conserved_gal_genes')
