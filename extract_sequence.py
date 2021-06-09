from Bio import SeqIO


def get_gene(path_to_fasta, start, stop, strand):
    """
    Function extract gene sequence from fasta file.
    Useful for BLAST searching.
    :param path_to_fasta: path to genome FASTA file
    :param start: start coordinate of desired fragment
    :param stop: end coordinate of desired fragment
    :param strand: '+' or '-' DNA strand on which fragment is located
    :return: sequence of desired DNA fragment
    """
    for fasta in SeqIO.parse(path_to_fasta, "fasta"):
        sequence = fasta[start: stop]
        if strand == '+':
            return sequence.seq
        else:
            return sequence.reverse_complement().seq

 
def get_fragment_sequence(start_, end_, strand, genome_fasta,
                          up_=300, down_=300):
    """
    Function find a sequence of specified DNA fragment (e.g., a gene).
    Useful for further checking promoter or terminator region on the obtained <upstream_gene_downstream> sequence.
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

