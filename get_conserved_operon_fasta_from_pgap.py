from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def get_seqrecord_conserved_gene(fna, gff, arr_id_coord_strand, gene_names, operon_strand):
    """

    :param fna:
    :param gff:
    :param arr_id_coord_strand:
    :param gene_names:
    :param operon_strand:
    :return:
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


def get_conserved_genes_seq(lst_gal_gff, edge_operon_gene, operon_length_, last=False):
    """

    :param lst_gal_gff:
    :param edge_operon_gene:
    :param last:
    :param operon_length_:
    :return:
    """
    lst_new_records = []
    for gff in lst_gal_gff:
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


def write_gene_fasta_files(gal_genes, lst_new_records, output_folder):
    """

    :param gal_genes:
    :param lst_new_records:
    :param output_folder:
    :return:
    """
    for gene in gal_genes:
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
