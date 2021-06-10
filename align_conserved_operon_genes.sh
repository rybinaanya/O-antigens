#!/bin/bash

echo -e "\n***Creating a folder to which .fasta files would be saved...\n"
mkdir $(pwd)/Providencia_selected/conserved_genes
echo -e "***Folder has been created\n"

echo -e «\n***Unzip all .gff files for each species\n»
for f in $(pwd)/Providencia_selected/*; do 
	gunzip ${f}/*.gff.gz; done

echo -e "***Starting to run python script get_conserved_operon_fasta.py to extract conserved wec operon genes\n"
for species in $(pwd)/Providencia_selected/GCF*; do echo "Looking for wec operon genes in ${species} assembly..."; python3 get_conserved_operon_fasta_2.py "${species}"/*.gff $(pwd)/Providencia_selected/conserved_genes wec; echo -e "Done with ${species} assembly\n"; done
echo -e "***All multi-fasta files with wec operon genes have been written and located in $(pwd)/Providencia_selected/conserved_genes\n"

echo -e "***Starting to run python script get_conserved_operon_fasta.py to extract conserved gal operon genes\n"

for f in $(pwd)/Providencia_selected/GCF*
do
    if grep -Fxq "galR" ${f}/*gff
    then
        echo "Looking for gal operon genes in ${species} assembly..."
        python3 get_conserved_operon_fasta_2.py $(pwd)/Providencia_selected/${species}/*.gff $(pwd)/Providencia_selected/conserved_genes gal
        echo -e "Done with ${species} assembly\n"
    else
        echo "gal operon not found in ${f} assembly"
    fi
done

echo -e "***All multi-fasta files with gal operon genes have been written and located in $(pwd)/Providencia_selected/conserved_genes\n"

echo -e "***Starting to run mafft for all multi-fasta files\n"
for f in $(pwd)/Providencia_selected/conserved_genes/*.fasta; do mafft ${f} > ${f%%.*}.mafft.fa; done
echo -e "\n*** mafft results have been saved at $(pwd)/Providencia_selected/conserved_genes\n"



