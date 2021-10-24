# -*- coding: utf-8 -*-
"""
Created on Mon May 31 18:02:57 2021

@author: Kieran Drake

This is the main script and is used to call functions and perform analysis

This script should be run in sections as and when required to perform analysis

"""

# Set home directory
home_path = input("Enter home directory: ")
# set as current directory
import os # required for changing filenames once in folder
os.chdir(home_path)
import contact_functions # module containing functions for specific tasks

######################################################
# Obtain list of TCR-pMHC complexes from www.IMGT.org
######################################################

#TCR_pMHC_list = contact_functions.TCR_pMHC_IMGT()
# Automation not working so search of website must be done manually and 
# files created manually e.g. IMGT_list_TR_MHC1.txt and IMGT_list_TR_MHC2.txt
PDB_code_list, TCR_pMHC1_list, TCR_pMHC2_list = contact_functions.manual_TCR_pMHC_IMGT()

##################################################################################
# Download PDB files for TCR-pMHC complexes
##################################################################################
PDB_file_dir = home_path + "PDB_files/"
PDB_files_confirm = contact_functions.get_PDB_files(PDB_code_list,home_path,PDB_file_dir)
print(PDB_files_confirm)

##################################################################################
# Extract various data from PDB files for each TCR-pMHC complex
##################################################################################
complex_dict, molecule_name_errors, all_aa_lengths = contact_functions.PDB_file_extract(PDB_code_list,TCR_pMHC1_list, TCR_pMHC2_list,home_path,PDB_file_dir)

# plot frequency of the number of amino acids in the chains
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
all_aa_lengths_array = np.asarray(all_aa_lengths)
plt.hist(all_aa_lengths_array) #, bins=np.arange(all_aa_lengths_array.min(), all_aa_lengths_array.max()+1))

###########################################################
# Categorise chains as ligand/peptide, TCR, MHC or other
###########################################################
ligand_limit = 11 # Anything below this will be treated as a ligand by PDBsum website
peptide_limit = 30 # 6v1a.pdb chains C (13 amino acids) and D work on format http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6v1a&chain1=D&chain2=C
complex_dict, molecule_list, other_molecule_list = contact_functions.categorise_molecule(PDB_code_list,complex_dict,ligand_limit,peptide_limit)

# Code to produce list of chains labelled as 'other' in function above so that 
# can be manually reviewed
chain_count = 0
for i in range(0,len(PDB_code_list),1):
    #print(PDB_code_list[i])
    no_of_chains = len(complex_dict['dict_'+ PDB_code_list[i]]['Chains_info_2'][0])
    for n in range(0,no_of_chains,1):
        #print(str(n+1) + ' of ' + str(no_of_chains) + ' chains')
        chain = complex_dict['dict_'+ PDB_code_list[i]]['Chains_info_2'][0][n][0]
        chain_length = int(complex_dict['dict_'+ PDB_code_list[i]]['Chains_info_2'][0][n][1])
        molecule_name = complex_dict['dict_'+ PDB_code_list[i]]['Chains_info_2'][0][n][2]
        molecule_category = complex_dict['dict_'+ PDB_code_list[i]]['Chains_info_2'][0][n][3]
        if molecule_category == 'other':
            chain_count = chain_count +1
            print(PDB_code_list[i],chain,chain_length,molecule_name,molecule_category)
print(chain_count)
# This printed list was then used in Excel for the manual recategorisation based 
# on visual inspection of the structures on the PDB website
# Once completed the .csv file 'Other_chains_recategorised.csv' was saved
# The format of the file is 
# PDB code,Chain,Chain length,Molecule Name,Molecule category assigned by 'categorise_molecule' function,Manual category assigned 

# Code for changing the chain/molecule category outside of the 'categorise_molecule' function
# This uses the file created manually as described above
import os
os.chdir(home_path)
filename = 'Other_chains_recategorised_tab.txt'
with open(filename,'r') as f:
    other_chains_new_categories = f.read().splitlines() 
for i in range(1,len(other_chains_new_categories),1):
    new_cat_list = other_chains_new_categories[i].split('\t')
    PDB_code = new_cat_list[0]
    chain = new_cat_list[1]
    new_molecule_category = new_cat_list[5]
    for m in range(0,len(complex_dict['dict_'+ PDB_code]['Chains_info_2'][0]),1):
        if complex_dict['dict_'+ PDB_code]['Chains_info_2'][0][m][0] == chain:
            # Need to convert from tuple to list and back to tuple before replacing in dictionary
            t = complex_dict['dict_'+ PDB_code]['Chains_info_2'][0][m]
            lst = list(t)
            lst[3] = new_molecule_category
            t = tuple(lst)
            complex_dict['dict_'+ PDB_code]['Chains_info_2'][0][m] = t

# The code above does not work for PDB code 2e71 due to Excel interpreting 
# as an exponential so can be done manually with the code below

# Code for manually changing the chain/molecule category outside of the 'categorise_molecule' function
# Only need to change PDB_code, chain, and new molecule category
PDB_code = '2e7l'
chain = ['C','D']
new_molecule_category = ['TCR','TCR']
for n in range(0,len(chain),1):
    for m in range(0,len(complex_dict['dict_'+ PDB_code]['Chains_info_2'][0]),1):
        if complex_dict['dict_'+ PDB_code]['Chains_info_2'][0][m][0] == chain[n]:
            # Need to convert from tuple to list and back to tuple before replacing in dictionary
            t = complex_dict['dict_'+ PDB_code]['Chains_info_2'][0][m]
            lst = list(t)
            lst[3] = new_molecule_category[n]
            t = tuple(lst)
            complex_dict['dict_'+ PDB_code]['Chains_info_2'][0][m] = t


############################################################
# Add chain sequences to dictionary
############################################################

# Create dictionary of amino acids
dict_aa_3to1 = {}
dict_aa_1to3 = {}
three_letter_list = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS',
                     'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP',
                     'TYR','VAL','ASX','GLX']
one_letter_list = ['A','R','N','D','C','Q','E','G','H',
                     'I','L','K','M','F','P','S','T','W',
                     'Y','V','B','Z']
for n in range(0,len(three_letter_list),1):
    dict_aa_3to1[three_letter_list[n]] = one_letter_list[n]
    dict_aa_1to3[one_letter_list[n]] = three_letter_list[n]
        
# Search PDB files for sequences and add to dictionary of TCR-pMHC data
complex_dict,non_standard_list = contact_functions.add_sequence(PDB_code_list,home_path,complex_dict,PDB_file_dir,dict_aa_3to1,dict_aa_1to3)


############################################################
# Save dictionary to file
############################################################
import json
# Serialize data into file:
json.dump( complex_dict, open( "complex_dict.json", 'w' ) )
# Read data from file:
data = json.load( open( "complex_dict.json" ) )


###########################################################
# Create list of chain interactions to download from PDB Sum
###########################################################
complex_dict = contact_functions.chain_chain(PDB_code_list, complex_dict)
# Serialize data into file:
json.dump( complex_dict, open( "complex_dict.json", 'w' ) )


################################################
# Work out cut off between peptide and ligand on PDB Sum
################################################

# Produce list of PDB codes, chains and chain lengths sorted by chain length
chain_lengths = []
chains = []
PDB_list = []
chain_types = []
chain_names=[]
chain_sequences = []
MHC_class_list = []
resolution_list = []
for item in PDB_code_list:
    a = complex_dict['dict_'+item]['Chains_info_2'][0]
    MHC_class = complex_dict['dict_'+item]['MHC_class']
    resolution = float(complex_dict['dict_'+item]['Resolution'])
    for item2 in a:
        PDB_list.append(item)
        MHC_class_list.append(MHC_class)
        resolution_list.append(resolution)
        chains.append(item2[0])
        chain_lengths.append(item2[1])#int(item2[0][1]))
        chain_names.append(item2[2])
        chain_types.append(item2[3])
        chain_sequences.append(''.join(item2[5]))
output = list(zip(PDB_list,MHC_class_list, resolution_list,chains,chain_lengths,chain_names,chain_types,chain_sequences))
#from operator import itemgetter
#output.sort(key=itemgetter(3))
#output[0:10]

# write to file
import csv
with open('chain_info.csv','w') as out:
    csv_out=csv.writer(out)
    csv_out.writerow(['PDB code','MHC class','Resolution','Chain','Chain length','Chain name','Chain type','Chain sequence'])
    for row in output:
        csv_out.writerow(row)


################################################
# Download TCR-pMHC contact files from PDB sum
################################################
PDBsum_files_confirm, PDBsum_directory = contact_functions.get_pdbsum_files(PDB_code_list,home_path,complex_dict)
print(PDBsum_files_confirm)


################################################
# Analyse TCR-pMHC contact files from PDB sum
################################################
## If need to manually set a particular PDBsum file directory
#PDBsum_directory = input("Enter PDBsume directory: ")

# Create dataframe of information produced by extraction from PDBsum information and analysis
contact_df = contact_functions.analyse_PDBsum(PDB_code_list,home_path,complex_dict,PDBsum_directory)

# Add peptide starting residue number so can adjust the residue number
contact_df = contact_functions.add_residue_start(home_path,complex_dict,PDB_file_dir,contact_df)

# Write data frame to csv file
contact_df.to_csv('contact_df.csv', sep='\t', encoding='utf-8')

# Add peptide sequence and length to dataframe
contact_df_final,pep_seq_df = contact_functions.add_pep_seq(PDB_code_list,home_path,complex_dict,PDB_file_dir,contact_df)

# Adjust residue numbers for peptides of length 8 from 1-8 to 2-9
import numpy as np

contact_df_final['Residue no. adjusted'] = np.where(contact_df_final['Peptide chain length'] != 8,
                contact_df_final['Residue no. adjusted'],
                contact_df_final['Residue no. adjusted'] + 1)

# Write data frame to csv file
contact_df_final.to_csv('contact_df_final.csv', sep='\t', encoding='utf-8')
pep_seq_df.to_csv('pep_seq_df.csv', sep='\t', encoding='utf-8')

# count number of complexes with peptides of different length
import pandas as pd
# store in temp dataframe
df_temp = contact_df_final
complex_list = pd.unique(df_temp['PDB code'])
count = 0
pep_length_list = []
for item in complex_list: 
    count += 1
    for n in range(0,len(pep_seq_df),1):
        if pep_seq_df['PDB code'][n] in pep_length_list_PDB:
            #do nothing
        else:
            if pep_seq_df['PDB code'][n] == item:
                pep_length_list_PDB.append(item) # create list of PDB codes that have searched pep_seq_df
                pep_length_list.append(pep_seq_df['Peptide chain length'][n]) # create list of PDB codes that have searched pep_seq_df
# dataframe
peptide_length_df = pd.DataFrame({'PDB code': pep_length_list_PDB[:, 0],
            'Peptide length': pep_length_list[:, 1],})

    
################################################
# Count number of complexes with more than one peptide
################################################

# cycles through PDB codes and prints information for complexes with 
# more than 1 peptide in the structure
counter = 0
for PDB_code in PDB_code_list:
    number_of_ligand = len(complex_dict['dict_'+PDB_code]['ligand'])
    number_of_pep = len(complex_dict['dict_'+PDB_code]['peptide'])
    if number_of_pep > 1:
        MHC_class = complex_dict['dict_'+PDB_code]['MHC_class']
        resolution = complex_dict['dict_'+PDB_code]['Resolution']
        counter += 1
        pep_chain = complex_dict['dict_'+PDB_code]['peptide'][0]
        for i in range(0, len(complex_dict['dict_'+PDB_code]['Chains_info_2'][0]),1):
            if pep_chain == complex_dict['dict_'+PDB_code]['Chains_info_2'][0][i][0]:
                pep_chain_length = complex_dict['dict_'+PDB_code]['Chains_info_2'][0][i][1]
        print(counter, PDB_code, MHC_class, resolution, number_of_pep, pep_chain_length)
    if number_of_ligand > 1:
        MHC_class = complex_dict['dict_'+PDB_code]['MHC_class']
        resolution = complex_dict['dict_'+PDB_code]['Resolution']
        counter += 1
        pep_chain = complex_dict['dict_'+PDB_code]['ligand'][0]
        for i in range(0, len(complex_dict['dict_'+PDB_code]['Chains_info_2'][0]),1):
            if pep_chain == complex_dict['dict_'+PDB_code]['Chains_info_2'][0][i][0]:
                pep_chain_length = complex_dict['dict_'+PDB_code]['Chains_info_2'][0][i][1]
        print(counter, PDB_code, MHC_class, resolution, number_of_ligand, pep_chain_length)
        

##############################################################################
# Create dataframes filtered by different parameters
##############################################################################

# store in temp dataframe
df_temp = contact_df_final

# Filter by resolution = 3A or below
resolution = 3
df_res3 = df_temp.loc[df_temp['Resolution']<=resolution]
df_res3.to_csv('df_res3.csv', sep='\t', encoding='utf-8')

# Filter by resolution = 2.5A or below
resolution = 2.5
df_res2_5 = df_temp.loc[df_temp['Resolution']<=resolution]
df_res2_5.to_csv('df_res2_5.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I
MHC_class = 1
df_MHC1 = df_temp.loc[df_temp['MHC class']==MHC_class]
df_MHC1.to_csv('df_MHC1.csv', sep='\t', encoding='utf-8')

# Filter by MHC class II
MHC_class = 2
df_MHC2 = df_temp.loc[df_temp['MHC class']==MHC_class]
df_MHC2.to_csv('df_MHC2.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I and resolution = 3A or below
resolution = 3; MHC_class = 1
df_res3_MHC1 = df_temp.loc[(df_temp['Resolution']<=resolution) & (df_temp['MHC class']==MHC_class)]
df_res3_MHC1.to_csv('df_res3_MHC1.csv', sep='\t', encoding='utf-8')

# Filter by MHC class II and resolution = 3A or below
resolution = 3; MHC_class = 2
df_res3_MHC2 = df_temp.loc[(df_temp['Resolution']<=resolution) & (df_temp['MHC class']==MHC_class)]
df_res3_MHC2.to_csv('df_res3_MHC2.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I and resolution = 2.5A or below
resolution = 2.5; MHC_class = 1
df_res25_MHC1 = df_temp.loc[(df_temp['Resolution']<=resolution) & (df_temp['MHC class']==MHC_class)]
df_res25_MHC1.to_csv('df_res25_MHC1.csv', sep='\t', encoding='utf-8')

# Filter by MHC class II and resolution = 2.5A or below
resolution = 2.5; MHC_class = 2
df_res25_MHC2 = df_temp.loc[(df_temp['Resolution']<=resolution) & (df_temp['MHC class']==MHC_class)]
df_res25_MHC2.to_csv('df_res25_MHC2.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I and 9mer or less peptides
MHC_class = 1; peptide_chain_length = 9
df_MHC1_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class)]
df_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I, 9mer or less peptides and resolution = 3A or lower
MHC_class = 1; peptide_chain_length = 9; resolution = 3
df_res3_MHC1_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Resolution']<=resolution)]
df_res3_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I, 9mer or less peptides and resolution = 2.5A or lower
MHC_class = 1; peptide_chain_length = 9; resolution = 2.5
df_res25_MHC1_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Resolution']<=resolution)]
df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')


# Using dataframe containing TCR-pMHC complexes where the MHC/HLA is HLA-A*0201
# This is created in the section below ***# Finding complexes with MHC HLA 201***
df_temp = df_hla_a_0201
# Filter by MHC class I (HLA-A*0201) and 9mer or less peptides
MHC_class = 1; peptide_chain_length = 9
df_hla_a_0201_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class)]
#df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

# Filter by MHC class I (HLA-A*0201) and 10mer peptides
df_temp = df_hla_a_0201
MHC_class = 1; peptide_chain_length = 10
df_hla_a_0201_p10 = df_temp.loc[(df_temp['Peptide chain length']==peptide_chain_length) & (df_temp['MHC class']==MHC_class)]
#df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

# Class
# remove multi peptide per PDB
df_exclusion_check = df_temp.loc[df_temp['PDB code'] == '6avg']
df_exclusion_check = df_temp.loc[df_temp['PDB code'] == '2ak4']


################################################
# Produce Logo plots
################################################
# https://logomaker.readthedocs.io/en/latest/examples.html#crp-energy-logo

pip install logomaker
import logomaker
import matplotlib.pyplot as plt
#logomaker.demo('fig1b')
# Choose option for y-axis plot
logo_option = 'number'
logo_option = 'probability'

# Include all peptide chains for each PDB code or just one 
chain_option = 'all'
chain_option = 'unique'

# Choose which dataframe to use in logo plot
df = contact_df_final
df = df_res3
df = df_res2_5
df = df_MHC1
df = df_MHC2
df = df_res3_MHC1
df = df_res3_MHC2
df = df_res25_MHC1
df = df_res25_MHC2
df = df_MHC1_p9
df = df_res3_MHC1_p9
df = df_res25_MHC1_p9
df = df_hla_a_0201_p9
df = df_hla_a_0201_p10

# Create logo plot
logo_image = contact_functions.make_logo_plot(df,logo_option,chain_option)


################################################
# Logo plots for section 4.1 in project report
################################################
# https://logomaker.readthedocs.io/en/latest/examples.html#crp-energy-logo

MHC_class = 1; peptide_chain_length = 50; chain_1_type = 'MHC'
df_temp = contact_df_final
df_temp2 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Chain 1 type']==chain_1_type)]

# Complexes with top 2 p-TCR contact positions: 8, 9
PDB_codes_logo_89 = ['6d78']
PDB_codes_logo_78 = ['1fo0','2e7l','2oi9','3e2h','3tpu','4mvb','4n0c','5sws','5swz','5wlg']
PDB_codes_logo_68 = ['2jcc','2uwe','2vlr','5m00','5m01']
PDB_codes_logo_67 = ['1mi5','3sjv','4n5e']
PDB_codes_logo_58 = ['1kj2','4mji','1ao7','1bd2','1qrn','2esv','2gj6','2vlj','2vlk','3d39','3d3v','3gsn','3h9s','3pwp','3qfj','3tf7','4ftv','5euo','5hhm','5hho','5isz','5jhd','5w1v','5w1w','6rpb','6rsy']
PDB_codes_logo_57 = ['1g6r','1jtr','1mwa','1nam','2ckb','2ol3','1qse','5tez','6eqb']
PDB_codes_logo_56 = ['1lp9','2j8u']
PDB_codes_logo_48 = ['3qdj','3tfk','3tjh','4eup','4mnq','4ms8','4mxq','5eu6','5m02','5til','5tje','6g9q','6rp9']
PDB_codes_logo_47 = ['3kpr','3kps','3qeq','6r2l']
PDB_codes_logo_46 = ['3o4l','4qrp','6l9l']
PDB_codes_logo_45 = ['2bnq','2bnr','2f53','2f54','2p5e','2p5w','2pye','5d2n','5nme','5nmf','5nmg','5xot','6bj2','6bj3','6bj8','6q3s','6rpa']
PDB_codes_logo_37 = ['6mtm']
PDB_codes_logo_138 = ['6vmx']
PDB_codes_logo_15 = ['1qsf','5bs0','6eqa']
PDB_codes_logo_14 = ['5brz']


df_logo = df_temp.loc[df_temp['PDB code'].isin(PDB_codes_logo_14)] 

pip install logomaker
import logomaker
import matplotlib.pyplot as plt
#logomaker.demo('fig1b')
#logo_option = 'number'
logo_option = 'probability'

#chain_option = 'all'
chain_option = 'unique'

logo_image = contact_functions.make_logo_plot(df_logo,logo_option,chain_option)

##################################################
# Analysing peptide-MHC contacts to determine which peptide residue positions 
# in groove for peptides with 10 or more amino acids
##################################################

# Re-run chain-chain info compiling in complex_dict as some 
# molecules changed from 'Other' category to a particular molecule type eg MHC
complex_dict = contact_functions.chain_chain(PDB_code_list, complex_dict)
# Serialize data into file:
json.dump( complex_dict, open( "complex_dict.json", 'w' ) )

# Download PDBSum files again so have complete list
PDBsum_files_confirm, PDBsum_directory = contact_functions.get_pdbsum_files(PDB_code_list,home_path,complex_dict)
print(PDBsum_files_confirm)

# Analyse PDBsum files - now with additional MHC molecules and chains labelled
## If need to manually set a particular PDBsum file directory
#PDBsum_directory = input("Enter PDBsume directory: ")
contact_df = contact_functions.analyse_PDBsum(PDB_code_list,home_path,complex_dict,PDBsum_directory)

# Add peptide starting residue number so can adjust the residue number
contact_df = contact_functions.add_residue_start(home_path,complex_dict,PDB_file_dir,contact_df)

# Write data frame to csv file
contact_df.to_csv('contact_df_2021_10_10_inc_MHC.csv', sep='\t', encoding='utf-8')

# Add peptide sequence and length to dataframe
contact_df_final,pep_seq_df = contact_functions.add_pep_seq(PDB_code_list,home_path,complex_dict,PDB_file_dir,contact_df)

# Adjust residue numbers for peptides of length 8 from 1-8 to 2-9
import numpy as np

contact_df_final['Residue no. adjusted'] = np.where(contact_df_final['Peptide chain length'] != 8,
                contact_df_final['Residue no. adjusted'],
                contact_df_final['Residue no. adjusted'] + 1)

# Write data frame to csv file
contact_df_final.to_csv('contact_df_final_2021_10_10_incMHC.csv', sep='\t', encoding='utf-8')
pep_seq_df.to_csv('pep_seq_df.csv', sep='\t', encoding='utf-8')

# count number of complexes with peptides of different length
import pandas as pd
# store in temp dataframe
df_temp = contact_df_final
complex_list = pd.unique(df_temp['PDB code'])
count = 0
pep_length_list = []
pep_length_list_PDB = []
for item in complex_list: 
    count += 1
    for n in range(0,len(pep_seq_df),1):
        if pep_seq_df['PDB code'][n] in pep_length_list_PDB: #pep_length_list_PDB
            #do nothing
            something = 0
        else:
            if pep_seq_df['PDB code'][n] == item:
                pep_length_list_PDB.append(item) # create list of PDB codes that have searched pep_seq_df
                pep_length_list.append(pep_seq_df['Peptide chain length'][n]) # create list of PDB codes that have searched pep_seq_df
# dataframe
#peptide_length_df = pd.DataFrame({'PDB code': pep_length_list_PDB[:, 0],
#            'Peptide length': pep_length_list[:, 1],})
peptide_length_df = pd.DataFrame({'PDB code': pep_length_list_PDB,'Peptide length': pep_length_list})

###########################################################
# Code to filter dataframe
###########################################################
import numpy as np

# store in temp dataframe
df_temp = contact_df_final

MHC_class = 1; peptide_chain_length = 9; #chain_1_type = 'TCR'
df_MHC1_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class)]
df_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'TCR'
df_MHC1_TCR_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Chain 1 type']==chain_1_type)]
df_MHC1_TCR_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')
df_MHC1_TCR_p9_stats = df_MHC1_TCR_p9.describe()

MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'TCR'; resolution = 3
df_MHC1_TCR_p9_res3 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Chain 1 type']==chain_1_type) & (df_temp['Resolution']<=resolution)]
df_MHC1_TCR_p9_res3.to_csv('df_MHC1_p9_res3.csv', sep='\t', encoding='utf-8')
df_MHC1_TCR_p9_res3_stats = df_MHC1_TCR_p9_res3.describe()

MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'TCR'; resolution = 2.5
df_MHC1_TCR_p9_res25 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Chain 1 type']==chain_1_type) & (df_temp['Resolution']<=resolution)]
df_MHC1_TCR_p9_res25.to_csv('df_MHC1_p9_res25.csv', sep='\t', encoding='utf-8')
df_MHC1_TCR_p9_res25_stats = df_MHC1_TCR_p9_res25.describe()

MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'MHC'; resolution = 2.5
df_MHC1_MHC_p9_res25 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Chain 1 type']==chain_1_type) & (df_temp['Resolution']<=resolution)]
df_MHC1_MHC_p9_res25.to_csv('df_MHC1_MHC_p9_res25.csv', sep='\t', encoding='utf-8')
df_MHC1_MHC_p9_res25_stats = df_MHC1_MHC_p9_res25.describe()

MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'MHC'
df_MHC1_MHC_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['MHC class']==MHC_class) & (df_temp['Chain 1 type']==chain_1_type)]
df_MHC1_MHC_p9.to_csv('df_MHC1_MHC_p9.csv', sep='\t', encoding='utf-8')
df_MHC1_MHC_p9_stats = df_MHC1_MHC_p9.describe()

df = df_MHC1_TCR_p9_res25
df = df_MHC1_TCR_p9
df = df_MHC1_MHC_p9_res25
df = df_MHC1_MHC_p9

df_temp = df_hla_a_0201
MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'MHC'
df_hla_a_0201_MHC_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['Chain 1 type']==chain_1_type)]
#df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

df_temp = df_hla_a_0201
MHC_class = 1; peptide_chain_length = 10; chain_1_type = 'MHC'
df_hla_a_0201_MHC_p10 = df_temp.loc[(df_temp['Peptide chain length']==10) & (df_temp['Chain 1 type']==chain_1_type)]
#df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

df_temp = df_hla_a_0201
MHC_class = 1; peptide_chain_length = 9; chain_1_type = 'TCR'
df_hla_a_0201_TCR_p9 = df_temp.loc[(df_temp['Peptide chain length']<=peptide_chain_length) & (df_temp['Chain 1 type']==chain_1_type)]
#df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

df_temp = df_hla_a_0201
MHC_class = 1; peptide_chain_length = 10; chain_1_type = 'TCR'
df_hla_a_0201_TCR_p10 = df_temp.loc[(df_temp['Peptide chain length']==10) & (df_temp['Chain 1 type']==chain_1_type)]
#df_res25_MHC1_p9.to_csv('df_MHC2_p9.csv', sep='\t', encoding='utf-8')

df = df_hla_a_0201_MHC_p9
df = df_hla_a_0201_MHC_p10
df = df_hla_a_0201_TCR_p9
df = df_hla_a_0201_TCR_p10

# TCR chart
summary_contact_df,summary_contact_stats = contact_functions.analyse_contacts(df)
#create box plot
import matplotlib as plt
import matplotlib.pyplot
ax = summary_contact_df.plot(kind='box',ylim=(0,50))
ax.set_xlabel('Peptide residue position')
ax.set_ylabel('Number of peptide-TCR contacts (non-backbone)')
plt.show()

# MHC chart
summary_contact_df,summary_contact_stats = contact_functions.analyse_contacts(df)
# create box plot
import matplotlib as plt
import matplotlib.pyplot
ax = summary_contact_df.plot(kind='box',ylim=(0,50))
ax.set_xlabel('Peptide residue position')
ax.set_ylabel('Number of peptide-MHC contacts (including backbone)')
plt.show()

# create table
from pandas.plotting import table
table(ax, np.round(summary_contact_df.describe(), 2),loc='upper right')

##############################################################################
# Finding complexes with MHC HLA 201
##############################################################################

# List of complexes with HLA-A 0201 as the MHC molecule
list_hla_a_0201_PDBs = ['1ao7','1bd2','1lp9','2f53','2gj6','3d39','3d3v','3h9s','3pwp','3qdg','3qdj','3qdm','3qeq','3qfj','3utt','4eup','4euq','4ftv','4l3e','4zez','5c0b','5e9d','5isz','5jzi','5tez','5yxn']
df_temp = contact_df_final

df_hla_a_0201 = df_temp.loc[df_temp['PDB code'].isin(list_hla_a_0201_PDBs)] 

for PDB_code in PDB_code_list:
    for item in complex_dict['dict_' + PDB_code]['Chains_info_2'][0]:
        print(item)

df = df_hla_a_0201

# TCR chart
summary_contact_df,summary_contact_stats = contact_functions.analyse_contacts(df)
#create box plot
import matplotlib as plt
import matplotlib.pyplot
ax = summary_contact_df.plot(kind='box')
ax.set_xlabel('Peptide residue position')
ax.set_ylabel('Number of peptide-TCR contacts (non-backbone)')
plt.show()

# MHC chart
summary_contact_df,summary_contact_stats = contact_functions.analyse_contacts(df)
#create box plot
import matplotlib as plt
import matplotlib.pyplot
ax = summary_contact_df.plot(kind='box')
ax.set_xlabel('Peptide residue position')
ax.set_ylabel('Number of peptide-MHC contacts (including backbone)')
plt.show()


# Gather TCR chains
for PDB_code in list_hla_a_0201_PDBs:
    num_chains = len(complex_dict['dict_' + PDB_code]['Chains_info_2'][0])
    for n in range (0,num_chains,1):
        chain_type = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][3]
        if chain_type == 'TCR':
            seq_TCR = ''.join(complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][5])
            name = ''.join(complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][2])
            print(PDB_code+' '+chain_type+' '+name+' '+seq_TCR)
            print(PDB_code,chain_type,name,seq_TCR)


###############################################################################
# Look for identical chains in complexes
###############################################################################

# Returns data on chains that occur more than once in dataset, 
# with a file each of peptide, MHC and TCR molecule types
df_seq_peptide, df_seq_TCR, df_seq_MHC = contact_functions.identical_chains(PDB_code_list,complex_dict)

df_check_seq_peptide, df_check_MHC_peptide, df_check_TCR_peptide = contact_functions.identical_chains_matrix(PDB_code_list,complex_dict)


###############################################################################
# Section 4.3
###############################################################################
PDB_list = ['1jtr','1mwa','2ckb']
PDB_list = ['2vlj','2vlk','2vlr','5euo','5isz','5jhd','5tez']
PDB_list = ['5w1v']

# Gather TCR chains for section 4.3
for PDB_code in PDB_list:
    num_chains = len(complex_dict['dict_' + PDB_code]['Chains_info_2'][0])
    for n in range (0,num_chains,1):
        chain_type = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][3]
        #if chain_type == ('TCR' or 'MHC'):
        seq = ''.join(complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][5])
        name = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][2]
        #print(PDB_code+' '+chain_type+' '+name+' '+seq_TCR)
        print(n,"#",PDB_code,"#",chain_type,"#",name,"#",seq)


###############################################################################
# Create plots for individual complexes
###############################################################################
import numpy as np

# store in temp dataframe
df_temp = contact_df_final

#MHC_class = 1; peptide_chain_length = 9; #
PDB_code = '5w1v'; chain_1_type = 'TCR'
df_5w1v = df_temp.loc[(df_temp['PDB code']==PDB_code) & (df_temp['Chain 1 type']==chain_1_type)]
df_5w1v.to_csv('df_5w1v.csv', sep='\t', encoding='utf-8')

df = df_5w1v
# TCR chart
summary_contact_df,summary_contact_stats = contact_functions.analyse_contacts(df)
#create box plot
import matplotlib as plt
import matplotlib.pyplot
ax = summary_contact_df.plot(kind='bar')
ax.set_xlabel('Peptide residue position')
ax.set_ylabel('Number of peptide-TCR contacts (non-backbone)')
plt.show()

# MHC chart
summary_contact_df,summary_contact_stats = contact_functions.analyse_contacts(df)
#create box plot
import matplotlib as plt
import matplotlib.pyplot
ax = summary_contact_df.plot(kind='box')
ax.set_xlabel('Peptide residue position')
ax.set_ylabel('Number of peptide-MHC contacts (including backbone)')
plt.show()