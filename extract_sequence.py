from Bio import SeqIO


def extract_sequence(path_to_fasta, start, stop, strand):
    for fasta in SeqIO.parse(path_to_fasta, "fasta"):

        # 300 nucleotides added to every border to search for promoter and terminator
        sequence = fasta[start - 300: stop + 300]
        if strand == '+':
            return sequence.seq
        else:
            return sequence.reverse_complement().seq
