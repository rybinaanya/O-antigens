from Bio import SeqIO

 
def get_sequence_fragment(start_, end_, strand, genome_fasta,
                          up_=0, down_=0):
    """
    Function find a sequence of specified DNA fragment (e.g., a gene).
    By default extracts a sequence strictly from start to end coordinates.
    To get a extended sequence ( e.g. a gene with upstream and/or downstream regions),
    change arguments up_ and/or down_ to required number of nucleotides. 
    Useful for further checking promoter or terminator region on the obtained <upstream_gene_downstream> sequence
    or for validating annotation of a gene.
    :param start_pos: int, start coordinate of desired fragment
    :param end_pos: int, end coordinate of desired fragment
    :param strand: '+' or '-' DNA strand on which fragment is located
    :param genome_fasta: path to genome FASTA file
    :param up_: int specifies how much nucleotide to take additionally upstream of desired fragment
    :param down_: int specifies how much nucleotide to take additionally downstream of desired fragment
    :return: sequence of desired DNA fragment
    """
    for record in SeqIO.parse(genome_fasta, 'fasta'):
        if strand == "+":
            return record.seq[(start_ - up_):(end_ + down_)]
        else:

            return record.seq[(start_ - up_):(end_ + down_)].reverse_complement()
        
