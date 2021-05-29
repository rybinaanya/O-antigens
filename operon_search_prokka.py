O_antigens_genes = ["cpxA", "qdtf", "rmlA", "qdtA", "qdtB", "wzx", "wpaA",
                     "wpaA", "wzyE", "wpaB", "wpaC", "ugd", "wpaD",
                     "gne", "galE", "wza", "wzb", "wzc", "yibK", "wecA",
                     "wecB", "wecD", "wecE", "wecF", "wecG", "rffG", "rffH",
                     "wbpB", "wbpD", "wbpE", 'mshA', 'rfbM', 'rfaQ']


def find_O_antigens_genes(prokka_gff):
    """
    Function print info (start and end coordinates, strand and gene name)
    about candidate O-antigen operon genes from GFF annotation using list
    of O-antigen operon genes defined from literature
    :param prokka_gff: path to GFF annotation file (Prokka)
    """
    with open(prokka_gff) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            if "gene" in line:
                line = line.strip().split()
                info = line[8]
                start = line[3]
                strand = line[6]
                stop = line[4]
                for gene in O_antigens_genes:
                    if gene in info:
                        for gene_info in line[8].strip().split(';'):
                            if gene_info.startswith("gene="):
                                gene_name = gene_info
                                print(start, stop, strand, gene_name)
