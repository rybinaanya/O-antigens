import sys
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Specify arguments
quast_report = sys.argv[1]
assembly_report = sys.argv[2]
selected_assemblies = sys.argv[3]
download_gff = sys.argv[4]
download_fna = sys.argv[5]

# Load QUAST metadata results
df_quast = pd.read_csv(quast_report, sep='\t')

# Select assemblies which N50 is not higher than 25% quantile
lst_assembly_N50 = df_quast[df_quast['N50'] <= df_quast['N50'].quantile(q=.25)]['Assembly']

# Select assemblies which 'Total length (>= 50000 bp)' < 25% quantile
lst_assembly_TotLen = df_quast[df_quast['Total length (>= 50000 bp)']
                               < df_quast['Total length (>= 50000 bp)'].quantile(q=.25)]['Assembly']

# Unite assemblies of poor quality (by N50 and by 'Total length (>= 50000 bp)')
lst_assemblies_N50_TotLen = list(set(lst_assembly_N50).union(set(lst_assembly_TotLen)))

# Select good assemblies
df_quast_good = df_quast[~df_quast['Assembly'].isin(lst_assemblies_N50_TotLen)]


# Make in QUAST metadata additional column - 'Assembly ID'
df_quast_good['Assembly ID'] = df_quast_good['Assembly'].apply(lambda x: '_'.join(x.split('_')[0:2]))

# Save IDs of good assemblies as a list
lst_assemblies_good_IDs = df_quast_good['Assembly ID']

# NCBI summary data on Providencia
df_ncbi = pd.read_csv(assembly_report, sep='\t', header=None)

# Extract part of ncbi summary that corresponds to good assemblies
df_ncbi_good = df_ncbi.loc[df_ncbi.iloc[:, 0].isin(lst_assemblies_good_IDs)]

# Add a column with the full name of organism
df_ncbi_good['Strain'] = df_ncbi_good.iloc[:, 1].apply(lambda x: ' '.join(x.split(' ')[:2])) + ' ' + \
                         df_ncbi_good.fillna('').iloc[:, 2].apply(lambda x: x.split('strain=')[-1])


# Select candidate assemblies of high quality
# Get Assembly Accessions and Organism Scientific Name from NCBI metadata
# Make a dictionary: key - Assembly ID, value - Strain
dict_Strain_AssemlyID = dict(zip(df_ncbi_good[0], df_ncbi_good["Strain"]))

# Add a column 'Organism Scientific Name' to QUAST metadata based on dictionary
df_quast_good['Organism Scientific Name'] = df_quast_good['Assembly ID'].map(dict_Strain_AssemlyID)

# Add a column 'Species' to QUAST metadata based on 'Organism Scientific Name' column
df_quast_good['Species'] = df_quast_good['Organism Scientific Name'].apply(lambda x: ' '.join(x.split(' ')[0:2]))

# Get array of unique Species
array_uniq_sp = df_quast_good.sort_values('N50', ascending=False)['Species'].unique()

# Create dataframe where  '# contigs' and N50 values are sorted in ascending and descending order, respectively
df_quast_good_sorted = df_quast_good.sort_values('# contigs', ascending=True).sort_values('N50', ascending=False)

# Select 3 best assembly for each species and save respective indices
idx_lst = [i for sp in array_uniq_sp if sp != 'Providencia sp.'
           for i in df_quast_good_sorted.loc[df_quast_good_sorted['Species'] == sp]['Assembly ID'].head(3).index[:]]


# Save summary data on selected assemblies
df_selected = df_quast_good_sorted.loc[idx_lst][['Organism Scientific Name', 'Species', 'Assembly ID',
                                                 '# contigs', 'N50', 'L50', "Total length"]]

# Get FTP download link to GFF annotation files of selected assemblies
#df_ncbi['FTP download GFF link'] = 'wget' + ' ' + df_ncbi.loc[df_ncbi[0].isin(df_selected['Assembly ID'].to_list())].iloc[:, 3] + '_genomic.gff.gz'
#df_ncbi['FTP download GFF link'] = 'wget ' + df_ncbi.loc[df_ncbi[0].isin(df_selected['Assembly ID'].to_list())].iloc[:, 3].apply(lambda x: x + '/' + x.split('/')[-1] + '_genomic.gff.gz')
df_ncbi['FTP download GFF link'] = 'wget' + ' ' + df_ncbi.loc[df_ncbi[0].isin(df_selected['Assembly ID'].to_list())].iloc[:, 3]

# Get FTP download link to .fna genome files of selected assemblies
#df_ncbi['FTP download FNA link'] = 'wget' + ' ' + df_ncbi.loc[df_ncbi[0].isin(df_selected['Assembly ID'].to_list())].iloc[:, 3] + '_genomic.fna.gz'
df_ncbi['FTP download FNA link'] = 'wget' + ' ' + df_ncbi.loc[df_ncbi[0].isin(df_selected['Assembly ID'].to_list())].iloc[:, 3]

# Save FTP GFF download links to file
with open(download_gff, 'w') as f:
    f.write("#!/bin/bash" + "\n")
    for i in df_ncbi.dropna()['FTP download GFF link']:
        echo_line = 'echo -e \"***Creating a directory $(pwd)/' + ''.join(i.split('/')[-1].split('_')[:2]) + '\\n\"'
        f.write(echo_line + '\n')
        mkdir_line="mkdir " + ''.join(i.split('/')[-1].split('_')[:2])
        f.write(mkdir_line + '\n')
        echo_line = 'echo -e \"***Entering $(pwd)/' + ''.join(i.split('/')[-1].split('_')[:2]) + '\\n\"'
        f.write(echo_line + '\n')
        cd_line_in = "cd " + ''.join(i.split('/')[-1].split('_')[:2])
        f.write(cd_line_in + "\n")
        f.write('echo -e \"***Starting to download GFF annotation file\\n\"' + "\n")        
        f.write(i + '/' + i.split('/')[-1] + "_genomic.gff.gz\n")
        cd_line_out = "cd .."
        f.write('echo -e \"***Leaving the directory $(pwd)/' + ''.join(i.split('/')[-1].split('_')[:2]) + '\\n\"' + '\n')
        f.write(cd_line_out + "\n")
    f.write('echo -e \"***Download is complete \\n\"')
        

# Save FTP .fna download links to file
with open(download_fna, 'w') as f:
    f.write("#!/bin/bash" + "\n")
    for i in df_ncbi.dropna()['FTP download FNA link']:
        echo_line = 'echo -e \"***Entering $(pwd)/' + ''.join(i.split('/')[-1].split('_')[:2]) + ' \\n\"'
        f.write(echo_line + '\n')
        cd_line_in = "cd " + ''.join(i.split('/')[-1].split('_')[:2])
        f.write(cd_line_in + "\n")
        f.write('echo -e \"***Getting .fna genomic file\\n\"' + "\n")
        f.write(i + '/' + i.split('/')[-1] + "_genomic.fna.gz\n")
        f.write('echo -e \"***Leaving the directory $(pwd)/' + ''.join(i.split('/')[-1].split('_')[:2]) + '\\n\"' + "\n")
        cd_line_out = "cd .."
        f.write(cd_line_out + "\n")
    f.write('echo -e \"***Download is complete \\n\"')

# Save summary data on selected assemblies to file
df_selected.to_csv(selected_assemblies, index=False)


