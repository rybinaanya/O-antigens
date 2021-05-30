from Bio import SeqIO


# list of names of the studied assemblies
assemblies = ["GCA_000259175.1", "GCA_001874625.1", "GCA_010669105.1",
              "GCA_013702025.1", "GCA_013702185.1", "GCA_014652175.1",
              "GCA_016618195.1", "GCA_900455075.1", "GCA_900635875.1"]

# list of genes from conserved operons
genes_from_conserved_operons = ['wecA', 'wzzE', 'wecB', 'wecC', 'rffG', 'rffH', 
                                'wecD', 'wecE', 'wzxE', 'wecF', 'wzyE', 'wecG', 
                                'galE', 'galT', 'galK', 'galM', 'galR']


def record_multifasta(dir_with_assemblies, output_dir):
    """
    Write multi-FASTA files for each gene from the conserved operons
    for each of the assemblies
    :param dir_with_assemblies: path to folder with .fasta and .gff assemblies files
    :param output_dir: path to folder where resulting FASTA files would be written to
    :return: None
    """
    for gene in genes_from_conserved_operons:

        # create fasta file for each gene
        with open(f"{output_dir}/{gene}.fasta", 'w') as gene_to_multifasta:
            for assembly in assemblies:

                # find gene in gff
                with open(f"{dir_with_assemblies}/{assembly}.gff") as gff:
                    for line in gff:
                        if line.startswith("#"):
                            continue
                        if gene in line:
                            line = line.strip().split()

                            # extract sequence from fna by coordinates from gff
                            with open("%s.fna" % assembly) as assembly_fasta:
                                for fasta in SeqIO.parse(
                                        assembly_fasta, "fasta"):
                                    record = fasta[int(line[3])-1:int(line[4])]

                                    # write assembly name and seq to each gene multifasta file
                                    gene_to_multifasta.write(
                                        '>' + assembly + '\n')
                                    if line[6] == "+":
                                        gene_to_multifasta.write(
                                            str(record.seq) + '\n')
                                    else:
                                        gene_to_multifasta.write(
                                            str(record.reverse_complement().seq)
                                            + '\n')

record_multifasta("/O_ag_project/assemblies/",
                  '/O_ag_project/conserved_operons')
