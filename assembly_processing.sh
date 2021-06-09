#!/bin/bash

count=53

#----------- Step 1: get summary report on assemblies from RefSeq database -----------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Starting to download complete list of manually reviewed genomes (RefSeq database) to current directory: \n$(pwd)"
echo -e $(yes "#" | head -n${count})'\n'

wget ftp://ftp.ncbi.nlm.nih.gov//genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/395/GCF_000156395.1_ASM15639v1/GCF_000156395.1_ASM15639v1_genomic.fna.gz

echo -e "\n***Download is complete. File is located here: $(pwd)/assembly_summary_refseq.txt\n"


#-----------       Step 2: get FTP links to Providencia assemblies         ------------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Getting FTP download links to Providencia assemblies of full genome representation..."
echo -e $(yes "#" | head -n${count})'\n'

grep -E 'Providencia' assembly_summary_refseq.txt | grep -v 'phage' | awk -F "\t" '$14=="Full" {print $20}'  | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}'  - > download_fna_files.sh

echo -e "\n***FTP download links have been saved to file $(pwd)/download_fna_files.sh\n"


#-----------       Step 3: Download selected .fna genome files             ------------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Download Providencia .fna genome file to the current directory: \n$(pwd)"
echo -e $(yes "#" | head -n${count})'\n'

echo -e "\n***Creating a folder to which .fna genome files would be saved...\n"
mkdir Providencia_genomes
echo -e "***Folder $(pwd)/Providencia_genomes has been created\n"

echo -e "***Entering the folder $(pwd)/Providencia_genomes \n"

cd Providencia_genomes

echo -e "***Starting to run $(pwd)/test_download_fna_files.sh"
bash ../download_fna_files.sh
echo -e "\n***Finished running the script $(pwd)/download_fna_files.sh"

echo -e "***Leaving the directory $(pwd)/Providencia_genomes"
cd ..

echo -e "\n***Downloaded .fna genome files could be found here: $(pwd)/Providencia_genomes\n"


#-----------           Step 4: Assembly quality assessment                ------------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Evaluating quality of Providencia assemblies via QUAST"
echo -e $(yes "#" | head -n${count})'\n'

echo -e "***Starting to run QUAST ...\n"
quast ./Providencia_genomes/*gz -o Providencia_quast_results -t 10
echo -e "\n***QUAST results have been saved at $(pwd)/Providencia_quast_results\n"


#-----------         Step 5: Select best Providencia assemblies            ------------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Selecting best Providencia assemblies"
echo -e $(yes "#" | head -n${count})'\n'

echo -e "***Getting data on Providencia organism names and assembly IDs\n"
#grep -E 'Providencia' assembly_summary_refseq.txt | grep -v 'phage' | cut -f1,8,9,20 |  awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}'  - > Providencia_strains.txt
grep -E 'Providencia' assembly_summary_refseq.txt | grep -v 'phage' | cut -f1,8,9,20 > Providencia_strains.txt
echo -e "***Temporal file $(pwd)/Providencia_strains.txt has been created\n"

echo -e "***Starting to run python script select_assemblies.py\n"
python3 select_assemblies.py Providencia_quast_results/transposed_report.tsv Providencia_strains.txt selected_Providencia_data.csv download_selected_assemblies_gff.sh download_selected_assemblies_fna.sh
echo -e "***Finished running the script select_assemblies.py\n"
echo -e "***Removing temporal file $(pwd)/Providencia_strains.txt"
rm Providencia_strains.txt
echo -e  "\n***Data on selected Providencia assemblies is provided here: $(pwd)/selected_Providencia_data.csv\n"


#-----------Step 6: Get genomic .gff and .fna files for selected Providencia assemblies ------------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Download gff genomic annotation and genome fasta files of selected Providenica assemblies"
echo -e $(yes "#" | head -n${count})'\n'

echo -e "***Creating a directory to store genomic data on selected assemblies...\n"
mkdir Providencia_selected
echo -e "Directory $(pwd)/Providencia_selected has been created\n"

echo -e "Entering $(pwd)/Providencia_selected\n"
cd Providencia_selected

echo -e "***Starting to run $(pwd)/download_selected_assemblies_gff.sh...\n"
bash ../download_selected_assemblies_gff.sh

echo -e "***Starting to run $(pwd)/download_selected_assemblies_gff.sh...\n"
bash ../download_selected_assemblies_fna.sh



echo -e "***GFF and .fna genomic data on selected Providencia assemblies is located here:\n$(pwd)"

#-----------        Step 6: Annotate selected Providencia genomes using PROKKA          ------------#
echo -e '\n'$(yes "#" | head -n${count})
echo -e "Annotate selected Providencia genomes using PROKKA"
echo -e $(yes "#" | head -n${count})'\n'

echo -e "***Entering $(pwd)/Providencia_selected\n"
echo -e "***Starting to run PROKKA...\n"

for f in $(ls); do gunzip ${f}/*.fna.gz; prokka --outdir  ${f}/${f}_prokka --force --genus Providencia --prefix ${f} ${f}/*.fna; done

echo -e "*** Annotation via PROKKA is done. Check folders at $(pwd)"
# download selected annotations
# downnload fasta 
# running prokka
# create data on all assemblies that we downloaded name assembly link
# removing download files sh
### insert to sccipt that we gunzip gff and fna 
