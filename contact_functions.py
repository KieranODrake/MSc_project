#!/usr/bin/python
"""
Created on Sun May 23 11:14:26 2021

@author: Kieran Drake

Downloads list of TCR-pMHC PDB codes from IMGT.org

"""
#pip install wget
import wget # required for pulling in webpage
import os # required for changing filenames once in folder

def manual_TCR_pMHC_IMGT():
    """
    Extracts list of PDB codes from files manually created from search on 
        www.IMGT.org for TCR-pMHC complexes and returns a list of PDB codes
    :TCR_pMHC1_list: List of PDB codes for TCR-pMHC1 complexes
    :TCR_pMHC2_list: List of PDB codes for TCR-pMHC2 complexes
    :PDB_code_list: Combined list of PDB codes
    """
    #Opens file and reads in list of PDB codes
    # List of MHC class I complexes
    with open("IMGT_list_TR_MH1.txt") as f:
        TCR_pMHC1_list = f.readlines()
        TCR_pMHC1_list = [x.strip() for x in TCR_pMHC1_list]
        #print(TCR_pMHC1_list)
    
    # List of MHC class II complexes
    with open("IMGT_list_TR_MH2.txt") as f:
        TCR_pMHC2_list = f.readlines()
        TCR_pMHC2_list = [x.strip() for x in TCR_pMHC2_list]
        #print(TCR_pMHC2_list)
    
    PDB_code_list = TCR_pMHC1_list + TCR_pMHC2_list 
    print(PDB_code_list)
    
    return PDB_code_list, TCR_pMHC1_list, TCR_pMHC2_list
    
def get_pdbsum_files(PDB_code_list,home_path,complex_dict):
    """
    Downloads contact data from PDB Sum for a list of PDB codes and chain pairs
        Depending on the length of the peptides, they are treated as either
        ligands (10-mers and below) or peptides (11mers and above) by PDBsum, with different URLs used for each
        url format for ligand is http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=1ao7&ligtype=01&ligno=01
        url format for peptide is  http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6v1a&chain1=D&chain2=C 
    :PDB_code_list: List of 'PDB code', 'peptide chain', TCR chain'
        e.g. [['1fyt','C','E'],['1fytab','A','B']]
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    """
    #pip install wget
    import wget # required for pulling in webpage
    import os # required for changing filenames once in folder
    # Create new directory where PDBSum files will be saved each time
    from time import gmtime, strftime
    now = strftime("%Y-%m-%d %H-%M-%S", gmtime())
    PDBsum_base_dir = home_path + "PDBsum files/"
    directory = os.path.join(PDBsum_base_dir,now)
    os.makedirs(directory)
    # set as current directory
    os.chdir(directory)

    # Use list of TCR-pMHC complexes (e.g. 1fyt) to build web address on PDBSum
    # eg http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=1fyt&chain1=C&chain2=E

    # NOTE: The PDBSum web address format depends on the length of the peptide,
    # if it is less than 11 amino acids then PDBSum treats it as a ligand rather than a peptide chain
    # The actual webpage format is the same, although the pages for ligands include information on all chains
    # whereas the pages for peptides only have contact data for a pair of chains
        
    # Cycle through list of PDB codes
    count = 0
    for item in PDB_code_list:
        PDB_code = item
        if len(complex_dict['dict_' + PDB_code]['TCR_ligand']) != 0:
            count+=1
            # web address for ligand gives a single page for all contacts
            url = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetLigInt.pl?pdb=" + PDB_code + "&ligtype=01&ligno=01"
            # Download PDB Sum webpage to a file
            wget.download(url)
            # Rename the downloaded PDBSum file with the PDB code
            filename = 'PDBSum_'+ PDB_code +'_ligand.txt'
            os.rename(os.path.join(directory, "GetLigInt.pl"),os.path.join(directory, filename))
        if len(complex_dict['dict_' + PDB_code]['TCR_peptide']) != 0:
            for chain_pair in complex_dict['dict_' + PDB_code]['TCR_peptide']:
                count+=1
                chain1 = chain_pair[0:1]
                chain2 = chain_pair[1:2]
                url = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=" + PDB_code + "&chain1=" + chain1 + "&chain2=" + chain2
                # Download PDB Sum webpage to a file
                # webpage = wget.download(url)
                wget.download(url)
                # Rename the downloaded PDBSum file with the PDB code and TCR
                filename = 'PDBSum_'+ PDB_code + '_TCR_' + chain1 + '_pep_' + chain2 + '.txt'
                os.rename(os.path.join(directory, "GetIface.pl"),os.path.join(directory, filename))
        if len(complex_dict['dict_' + PDB_code]['MHC_peptide']) != 0:
            for chain_pair in complex_dict['dict_' + PDB_code]['MHC_peptide']:
                count+=1
                chain1 = chain_pair[0:1]
                chain2 = chain_pair[1:2]
                url = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=" + PDB_code + "&chain1=" + chain1 + "&chain2=" + chain2
                # Download PDB Sum webpage to a file
                wget.download(url)
                # Rename the downloaded PDBSum file with the PDB code and MHC
                filename = 'PDBSum_'+ PDB_code + '_MHC_' + chain1 + '_pep_' + chain2 + '.txt'
                os.rename(os.path.join(directory, "GetIface.pl"),os.path.join(directory, filename))
        
    # Count number of files downloaded from PDB Sum
    number_files = len(os.listdir(directory))
    PDBsum_files_confirm = str(number_files) + " out of " + str(count) + " PDBsum files downloaded"
    
    # Switch back to home directory    
    os.chdir(home_path)

    return PDBsum_files_confirm, directory

def get_PDB_files(PDB_code_list,home_path,PDB_file_dir):
    """
    Downloads PDB files for TCR-pMHC complexes in list from www.IMGT.org
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    :PDB_file_dir: Directory where PDB files will be saved
    """
    import os
    #### If want to create new directory with the curent time stamped on it
    # from time import gmtime, strftime
    # now = strftime("%Y-%m-%d %H-%M-%S", gmtime())
    # directory = os.path.join(PDB_file_dir, now)
    directory = os.path.join(PDB_file_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)
        # set as current directory
        os.chdir(directory)

    # download PDB files
    import urllib.request
    for n in PDB_code_list:
        PDB_code = n
        print(PDB_code)
        url = "https://files.rcsb.org/download/" + PDB_code + ".pdb"
        print(url)
        urllib.request.urlretrieve (url, directory+ "/" + PDB_code + ".pdb")
    
    # Count number of files downloaded from PDB website
    number_files = len(os.listdir(directory))
    PDB_files_confirm = str(number_files) + " out of " + str(len(PDB_code_list)) + " PDB files downloaded"
    # Switch back to home directory
    os.chdir(home_path)
    
    return PDB_files_confirm
    
def PDB_file_extract(PDB_code_list,TCR_pMHC1_list, TCR_pMHC2_list,home_path,PDB_file_dir):
    """
    Extracts information on chains from PDB files and saves in a dictionary of dictionaries
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :TCR_pMHC1_list: List of PDB codes for TCR-pMHC1 complexes
    :TCR_pMHC2_list: List of PDB codes for TCR-pMHC2 complexes
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    :PDB_file_dir: Directory where PDB files will be saved
    """
    # Extract information from PDB files and put into dictionary
    import re

    # set PDB file directory as current directory
    os.chdir(PDB_file_dir)
    # COULD ADD CODE TO SCAN DIRECTORY AND EXTRACT FROM ALL FILES IN DIR RATHER 
    # THAN JUST FROM THE LIST - IN CASE SOME FILES DON'T DOWNLOAD FOR SOME REASON

    # Create list of dictionaries
    dict_list=[]
    dict_list=list(range(len(PDB_code_list)))
    count = -1
    for n in range(len(PDB_code_list)):
        count += 1
        dict_list[count] = 'dict_' + PDB_code_list[count]

    # create dictionary of dictionaries    
    dictionaries = {}
    for d in range(0,len(dict_list),1):
        dictionaries[dict_list[d]]={}

    # Add information to dictionaries
    
    # Set searches
        
    # Find resolution in PDB file "REMARK   2 RESOLUTION.    2.69 ANGSTROMS."
    res = re.compile(r'REMARK\s+?\d\sRESOLUTION.\s+?(\d.?\d*?)\s*ANGSTROMS.')
    
    # Add list of chains and number of amino acids in each chain
    # p = re.compile(r'SEQRES   1 A  275  GLY SER HIS SER MET ARG TYR PHE PHE THR SER VAL SER          ')
    # Looking to return 'A' and '275' from the above line
    aa_num = re.compile(r'SEQRES\s{3}1\s*(\w)\s*(\d*)\s*')
    
    # Add molecule name for each chain
    # r = re.compile(r'COMPND   3 CHAIN: A;')
    # Assumes no more than 4 chains per molecule name
    chain = re.compile(r'COMPND\s*\d*\sCHAIN:\s(\w*),?\s?(\w*)?,?\s?(\w*)?,?\s?(\w*)?') #the ; is not included in the regex as it is not present in all lines
    # Match the molecule name
    # s = re.compile(r'COMPND   2 MOLECULE: HLA CLASS II HISTOCOMPATIBILITY ANTIGEN, DR ALPHA CHAIN;')
    name = re.compile(r'COMPND\s*\d*\sMOLECULE:\s(.*)') #the ; is not included in the regex as it is not present in all lines
    # Match 2nd line of molecule name if there is one
    # s = re.compile(r'COMPND   2 MOLECULE: HLA CLASS II HISTOCOMPATIBILITY ANTIGEN, DR ALPHA CHAIN;')
    name2 = re.compile(r'COMPND\s*\d*\s(.*)')
    # Match peptide start residue number
    # s = re.compile(r'COMPND   2 MOLECULE: HLA CLASS II HISTOCOMPATIBILITY ANTIGEN, DR ALPHA CHAIN;')
    residue_start_line = re.compile(r'^DBREF')
    
    # (re)initiate lists
    error_pdb = []
    error_chain = []
    molecule_name_errors = []
    all_aa_lengths = [] # used to gather info on aa length so can categorise as protein, MHC or TCR
    residue_start_list = []
    
    # Read in PDB files
    for n in range(0,len(PDB_code_list),1): # Cycle through the list of PDB codes for TCR-pMHC complexes
        # read PDB files
        filename = PDB_code_list[n] + '.pdb'
        m=len(PDB_code_list)
        print(filename + ": file " + str(n+1) + " of " + str(m))
        with open(filename,'r') as f:
            PDB_content = f.read().splitlines()    
    
        # Add PDB code to dictionaries
        dictionaries[dict_list[n]]['PDB_code'] = PDB_code_list[n]
        
        # (Re)initiate lists
        chain_list = []
        aa_length_list = []
        name_list = []
        chain_list_2 = []
        chains_info = []
        chain_length_list = []

        # Hold 4 lines in variables at a time as information sometimes on multiple lines
        for c in range(0,len(PDB_content)-3,1):
            PDB_content_line_1 = PDB_content[c]
            PDB_content_line_2 = PDB_content[c+1]
            PDB_content_line_3 = PDB_content[c+2]
            PDB_content_line_4 = PDB_content[c+3]
            # Find resolution and add to dictionary
            if res.search(PDB_content_line_1):
                match_res = res.search(PDB_content_line_1)
                dictionaries[dict_list[n]]['Resolution'] = match_res.group(1)
            # Find chain length and add to list
            if aa_num.search(PDB_content_line_1):
                match_length = aa_num.search(PDB_content_line_1)
                chain_list.append(match_length.group(1))
                aa_length_list.append(match_length.group(2))        
                all_aa_lengths.append(match_length.group(2))
            # Find molecule name and add to list
            # Some molecules have names that run on to more than one line 
            # e.g. Chain H in 1nam.pdb (2 lines) and Chains D & H in 4p5t.pdb (3 lines).
            # If this is the case then the regex will not find the name and so 
            # details are dumped into an error log/list
            # Need to determine how many lines the molecule name is on 
            if name.search(PDB_content_line_1): 
                if (name.search(PDB_content_line_1) and chain.search(PDB_content_line_2)):
                       match_chain = chain.search(PDB_content_line_2)
                       match_name = name.search(PDB_content_line_1)
                       for i in range(1,5,1): #COUlD ALTER SO USES MAX NUMBER OF GROUPS IN 'match_chain'
                            if match_chain.group(i) != '': # only performs if group is NOT empty
                                chain_list_2.append(match_chain.group(i))
                                name_list.append(match_name.group(1))
                elif (name.search(PDB_content_line_1) and chain.search(PDB_content_line_3)):
                       match_name1 = name.search(PDB_content_line_1)
                       match_name2 = name2.search(PDB_content_line_2)
                       match_chain = chain.search(PDB_content_line_3)
                       for i in range(1,5,1): #COUlD ALTER SO USES MAX NUMBER OF GROUPS IN 'match_chain'
                            if match_chain.group(i) != '': # only performs if group is NOT empty
                                chain_list_2.append(match_chain.group(i))
                                full_name = match_name1.group(1) + " " + match_name2.group(1)
                                name_list.append(full_name)
                elif (name.search(PDB_content_line_1) and chain.search(PDB_content_line_4)):
                       match_name1 = name.search(PDB_content_line_1)
                       match_name2 = name2.search(PDB_content_line_2)
                       match_name3 = name2.search(PDB_content_line_3)
                       match_chain = chain.search(PDB_content_line_4)
                       for i in range(1,5,1): #COUlD ALTER SO USES MAX NUMBER OF GROUPS IN 'match_chain'
                            if match_chain.group(i) != '': # only performs if group is NOT empty
                                chain_list_2.append(match_chain.group(i))
                                full_name = match_name1.group(1) + " " + match_name2.group(1) + " " + match_name3.group(1)
                                name_list.append(full_name)
                else:
                    name_list.append("Error - Check file")
                    error_pdb.append(PDB_code_list[n])
              
        # Add molecule names to list of chains and chain lengths
        chain_length_list = list(zip(chain_list,aa_length_list))
        chains_names = list(zip(chain_list_2,name_list))
        molecule_name_errors = error_pdb
        
        for k in range(0,len(chain_length_list),1):
            # Searches chains_names list for molecule name for chain letter in chain_length_list
            chain_names_entry = [item for item in chains_names if item[0] == chain_length_list[k][0]] 
            # creates three item tuple of chain letter, chain length and molecule name
            combine = list(zip(chain_length_list[k][0],chain_length_list[k][1:],chain_names_entry[0][1:]))
            chains_info.append(combine)
        chains_info.sort()

        dictionaries[dict_list[n]]['Chains_info'] = chains_info
        
    # Add information to dictionary about which type of MHC class (1 or 2) the complex contains
    for PDB_code in PDB_code_list:# cycle through PDB code dictionary entries
        if PDB_code in TCR_pMHC1_list:
            dictionaries['dict_' + PDB_code]['MHC_class'] = 1
        elif PDB_code in TCR_pMHC2_list:
            dictionaries['dict_' + PDB_code]['MHC_class'] = 2

    # Switch back to home directory    
    os.chdir(home_path)
    
    return dictionaries, molecule_name_errors, all_aa_lengths

def categorise_molecule(PDB_code_list,complex_dict,ligand_limit,peptide_limit):
    """
    Categorises molecules using a combination of chain length and molecule names.
        This categorisation is then added to the information in the dictionary
        for each chain for each PDB code
    :complex_dict: Dictionary containing information (chain, chain length, molecule name) for each TCR-pMHC complexes PDB code
    :ligand_limit: upper limit (not including) for classification as a ligand (necessary for PDB Sum website)
    :peptide_limit: upper limit (not including) for classification as a peptide
    """
    # Extract information from PDB files and put into dictionary
    import re

    # Create regexs to check if chain is a TCR
    name_TCR1 = re.compile(r'(?i)TCR') #(?i) indicates case insensitive
    #name_TCR2 = re.compile(r'tcr')
    name_TCR3 = re.compile(r'(?i)T[\s-]*(?i)Cell')
    
    # Create regexs to check if chain is an MHC (numbering of variables is 
    # not sequential as redundant regexes deleted)
    name_MHC1 = re.compile(r'(?i)MHC')
    name_MHC3 = re.compile(r'(?i)HISTOCOMPATIBILITY')
    name_MHC6 = re.compile(r'(?i)HLA')
    name_MHC7 = re.compile(r'(?i)BETA[\s-]2[\s-]MICROGLOBULIN')
    
    # Initiate lists
    chains_info_2 = []
    molecule_list = []
    other_molecule_list = []
    
    # Cycle through PDB codes/complexes
    for i in range(0,len(PDB_code_list),1):
        
        #(re)initiate lists
        chain_type_list = []
        molecule_name_list = []
        chain_length_list = []
        pdb_codes = []
        chain_list = []
        chains_info_2 =[]
        
        #cycle through dictionary for each pdb code
        for n in range(0,len(complex_dict['dict_'+ PDB_code_list[i]]['Chains_info']),1):
            chain = complex_dict['dict_'+ PDB_code_list[i]]['Chains_info'][n][0][0]
            chain_length = int(complex_dict['dict_'+ PDB_code_list[i]]['Chains_info'][n][0][1])
            molecule_name = complex_dict['dict_'+ PDB_code_list[i]]['Chains_info'][n][0][2]
#            residue_start = int(float(complex_dict['dict_'+ PDB_code_list[i]]['Chains_info'][n][0][3]))
            
            # Assign molecule type based on either chain length
            # or whether the name contains certain terms
            if chain_length < ligand_limit:
                chain_type = 'ligand'
            elif chain_length < peptide_limit:
                chain_type = 'peptide'
            elif (name_TCR1.search(molecule_name) or
                  name_TCR3.search(molecule_name)):
                chain_type = 'TCR'
            elif (name_MHC1.search(molecule_name) or
                  name_MHC3.search(molecule_name) or
                  name_MHC6.search(molecule_name) or
                  name_MHC7.search(molecule_name)):
                chain_type = 'MHC'
            else:
                chain_type = 'other'
            
            # Add the information to lists
            #print(PDB_code_list[i],chain,chain_length,molecule_name,chain_type)
            molecule_list.append(PDB_code_list[i]+';'+chain+';'+str(chain_length)+';'+molecule_name+';'+chain_type)
            chain_type_list.append(chain_type)
            molecule_name_list.append(molecule_name)
            chain_length_list.append(chain_length)
            chain_list.append(chain)

        # combine lists and add to dictionary
        combine = list(zip(chain_list,chain_length_list,molecule_name_list,chain_type_list))
        chains_info_2.append(combine)
        complex_dict["dict_"+ PDB_code_list[i]]['Chains_info_2'] = chains_info_2
        #print(complex_dict["dict_"+ PDB_code_list[i]]['Chains_info'][n][1])
        
        # create list of all molecules with information for checking
        
    # create list that just contains molecules categorised as 'other' 
    # so that they can be manually reviewed
    category_other = re.compile(r'other')
    for m in range(0,len(molecule_list),1):
        if category_other.search(molecule_list[m]):
            other_molecule_list.append(molecule_list[m])
    
    return complex_dict, molecule_list, other_molecule_list

def chain_chain(PDB_code_list, complex_dict):
    """
    Creates a list of chain pairs (peptide/ligand : TCR) or (peptide/ligand : MHC) 
        to search for contact information on PDB Sum 
    :complex_dict: Dictionary containing information (PDB code, resolution,
    chain, chain length, molecule name, molecule type) for each TCR-pMHC complex's PDB code
    """
    
    # cycle through dictionary and create list of peptide/ligand, TCR and MHC 
    # chains for each PDB code
    for i in range(0,len(PDB_code_list),1):
        
        #(Re)initiate lists
        TCR_peptide_list = []
        TCR_ligand_list = []
        TCR_chain_list = []
        MHC_peptide_list = []
        MHC_ligand_list = []
        MHC_chain_list = []
        peptide_chain_list = []
        ligand_chain_list = []
        
        for item in complex_dict['dict_'+ PDB_code_list[i]]['Chains_info_2'][0]:
            chain_letter = item[0]
            chain_category = item[3]
            if chain_category == 'TCR':
                TCR_chain_list.append(chain_letter)
            if chain_category == 'peptide':
                peptide_chain_list.append(chain_letter)
            if chain_category == 'ligand':
                ligand_chain_list.append(chain_letter)
            if chain_category == 'MHC':
                MHC_chain_list.append(chain_letter)
        
        complex_dict['dict_'+ PDB_code_list[i]]['TCR'] = TCR_chain_list
        complex_dict['dict_'+ PDB_code_list[i]]['peptide'] = peptide_chain_list
        complex_dict['dict_'+ PDB_code_list[i]]['ligand'] = ligand_chain_list
        complex_dict['dict_'+ PDB_code_list[i]]['MHC'] = MHC_chain_list
        
        # create chain-chain pairs and write to dictionary
        # cycle through TCRs...
        if len(complex_dict['dict_'+ PDB_code_list[i]]['TCR']) != 0:
            for TCR_item in complex_dict['dict_'+ PDB_code_list[i]]['TCR']:
                # ...and then through peptides...
                if len(complex_dict['dict_'+ PDB_code_list[i]]['peptide']) != 0:
                    for pep_item in complex_dict['dict_'+ PDB_code_list[i]]['peptide']:
                        TCR_chain = TCR_item
                        peptide_chain = pep_item
                        TCR_peptide_list.append(TCR_chain + peptide_chain)
                #...and ligands.
                if len(complex_dict['dict_'+ PDB_code_list[i]]['ligand']) != 0:
                    for lig_item in complex_dict['dict_'+ PDB_code_list[i]]['ligand']:
                        TCR_chain = TCR_item
                        ligand_chain = lig_item
                        TCR_ligand_list.append(TCR_chain + ligand_chain)
        
        # write the lists of chain-chain pairs to dictionary
        complex_dict['dict_'+ PDB_code_list[i]]['TCR_peptide'] = TCR_peptide_list
        complex_dict['dict_'+ PDB_code_list[i]]['TCR_ligand'] = TCR_ligand_list
        
        # create chain-chain pairs and write to dictionary
        # cycle through MHCs...
        if len(complex_dict['dict_'+ PDB_code_list[i]]['MHC']) != 0:
            for MHC_item in complex_dict['dict_'+ PDB_code_list[i]]['MHC']:
                # ...and then through peptides...
                if len(complex_dict['dict_'+ PDB_code_list[i]]['peptide']) != 0:
                    for pep_item in complex_dict['dict_'+ PDB_code_list[i]]['peptide']:
                        MHC_chain = MHC_item
                        peptide_chain = pep_item
                        MHC_peptide_list.append(MHC_chain + peptide_chain)
                #...and ligands.
                if len(complex_dict['dict_'+ PDB_code_list[i]]['ligand']) != 0:
                    for lig_item in complex_dict['dict_'+ PDB_code_list[i]]['ligand']:
                        MHC_chain = MHC_item
                        ligand_chain = lig_item
                        MHC_ligand_list.append(MHC_chain + ligand_chain)
        # write the lists of chain-chain pairs to dictionary
        complex_dict['dict_'+ PDB_code_list[i]]['MHC_peptide'] = MHC_peptide_list
        complex_dict['dict_'+ PDB_code_list[i]]['MHC_ligand'] = MHC_ligand_list
        
    return complex_dict

def analyse_PDBsum(PDB_code_list,home_path,complex_dict,PDBsum_directory):
    """
    Extracts information on chains from PDBsum files and saves in a dataframe
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    :complex_dict: Dictionary containing various information on the TCR-pMHC complexes 
        organised with a sub-dictionary for each PDB code/complex
    :PDBsum_directory: Directory containing files downloaded from PDBSum
    :contact_df: Dataframe containing information on peptide residue contacts 
        with TCRs (excludes backbone atoms) and MHCs (includes backbone atoms)
    """
    os.chdir(PDBsum_directory)
    
    #list of atoms in amino acid backbone
    backbone_list = ['CA','NH1','O','C','N','H','HA']
    
    # initiate dictionaries and lists
    import numpy as np
    PDBsum_contact_data = np.array(['PDB code','PDBsum file type','Resolution','Chain 1 type','Chain 1','Chain 2','Atom 2','Residue 2','Residue no.','Residue no. adjusted for start','Contact distance','MHC class','Contact type'])
    new_contact_data = []
    
    # Set searches
    import re
    HB = re.compile(r'Hydrogen bonds')
    HB_number = re.compile(r'Number of hydrogen bonds:\s+(\d+)')
    NBC = re.compile(r'Non-bonded contacts')
    NBC_number = re.compile(r'Number of non-bonded contacts:\s+(\d+)')
    SB = re.compile(r'Salt bridges')
    SB_number = re.compile(r'Number of salt bridges:\s+(\d+)')
    
    # This function takes c.20 minutes to run on c.1100 files so this tries to
    # estimate time remaining. It is not accurate as written curently because
    # the early files are analysed in a much shorter time that the longer files
    # and so the time remaining estimate is always an underestimate
    import time
    t0 = time.time()
    t1 = 0
    file_count = 0
    number_of_files = len(os.listdir(PDBsum_directory))

    for file in os.listdir(PDBsum_directory):
        file_count += 1
        t2 = t1
        t1 = time.time()
        total = t1-t0
        if file_count > 1:
            time_estimate = (total / (file_count -1)) * (number_of_files - file_count)
            print("File ",file_count," of ",number_of_files," in ",round(t1-t2,1)," seconds. Total time taken: ",round(total,1)," seconds. Estimated time to completion: ",round(time_estimate,1)," seconds")
        else:
            print("File ",file_count," of ",number_of_files)
        
        PDB_code = file[7:11]
        
        # Reset counter and line numbers 
        count = -1
        HB_line_number = ''
        NBC_line_number = ''
        SB_line_number = ''
        HB_start_line = ''
        NBC_start_line = ''
        SB_start_line = ''
        HB_end_line = ''
        NBC_end_line = ''
        SB_end_line = ''
    
        if file[12:18] == 'ligand':
            # NEED TO TREAT THESE DIFFERENTLY AS CONTAIN ALL THE CONTACT INFO
            # AND NOT JUST CONTACTS BETWEEN TWO CHAINS
            with open(file,'r') as f:
                PDBsum_content = f.read().splitlines()
                # Define the start and end lines of the contact information 
                # for each type of contact
                for line in PDBsum_content:
                    count += 1
                    if HB.search(line):
                        HB_line_number = count
                        HB_start_line = HB_line_number + 7
                    elif NBC.search(line):
                        NBC_line_number = count
                        NBC_start_line = NBC_line_number + 7
                    elif SB.search(line):
                        SB_line_number = count
                        SB_start_line = SB_line_number + 7
                    elif HB_number.search(line):
                        HB_number_match = HB_number.search(line)
                        HB_end_line = HB_start_line + int(HB_number_match.group(1))
                    elif NBC_number.search(line):
                        NBC_number_match = NBC_number.search(line)
                        NBC_end_line = NBC_start_line + int(NBC_number_match.group(1))
                    elif SB_number.search(line):
                        SB_number_match = SB_number.search(line)
                        SB_end_line = SB_start_line + int(SB_number_match.group(1))
                
                # Obtain information on hydrogen bond contacts
                if isinstance(HB_start_line, int):
                    for n in range(HB_start_line, HB_end_line,1):
                        # Example contact data line from 1d9k hydrogen bond
                        #        Atom Atom Res  Res              Atom Atom Res  Res
                        #         no. name name no.  Chain        no. name name no.  Chain  Distance
                        # '  1.    211  O   THR   28    A   <-->  4833  NH1 ARG  135    P      3.12'
                        #    g1    g2   g3   g4   g5    g6   g7   g8    g9   g10  g11  g12      g13
                        
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        chain1 = g6
                        chain2 = g12
                        # Record contact data if chain1 is a TCR
                        if (chain1+chain2) in complex_dict['dict_' + PDB_code]['TCR_ligand']:
                            atom2 = g9
                            # only save contact data if atom is not in the 
                            # backbone of amino acid in the peptide
                            if atom2 not in backbone_list:
                                PDBsum_filetype = 'ligand'
                                chain1_type = 'TCR'
                                residue2 = g10
                                residue2_number = int(g11)
                                distance = g13
                                residue2_number_adjusted = residue2_number # In PDBsum ligand files the residue number does not appear to need adjusting i.e. they begin at 1
                                MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                                resolution = complex_dict['dict_' + PDB_code]['Resolution']
                                contact_type = 'Hydrogen bond'
                                new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                                PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                        
                        # Record contact data if chain1 is a MHC
                        if (chain1+chain2) in complex_dict['dict_' + PDB_code]['MHC_ligand']:
                            atom2 = g9
                            ### Don't need to exclude peptide backbone atoms for 
                            ### MHC contacts, but could reinstate line below if 
                            ### want to exclude backbone atoms
                            #if atom2 not in backbone_list:
                            PDBsum_filetype = 'ligand'
                            chain1_type = 'MHC'
                            residue2 = g10
                            residue2_number = int(g11)
                            distance = g13
                            residue2_number_adjusted = residue2_number # In PDBsum ligand files the residue number does not appear to need adjusting i.e. they begin at 1
                            MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                            resolution = complex_dict['dict_' + PDB_code]['Resolution']
                            contact_type = 'Hydrogen bond'
                            new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                            PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                # Obtain information on non-bonded contacts
                if isinstance(NBC_start_line, int):
                    for n in range(NBC_start_line, NBC_end_line,1):
                        contact_data = PDBsum_content[n]
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        chain1 = g6
                        chain2 = g12
                        # Record contact data if chain1 is a TCR
                        if (chain1+chain2) in complex_dict['dict_' + PDB_code]['TCR_ligand']:
                            atom2 = g9
                            # only save contact data if atom is not in the 
                            # backbone of amino acid in the peptide
                            if atom2 not in backbone_list:
                                PDBsum_filetype = 'ligand'
                                chain1_type = 'TCR'
                                residue2 = g10
                                residue2_number = int(g11)
                                distance = g13
                                residue2_number_adjusted = residue2_number # In PDBsum ligand files the residue number does not appear to need adjusting i.e. they begin at 1
                                MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                                resolution = complex_dict['dict_' + PDB_code]['Resolution']
                                contact_type = 'Non-bonded contact'
                                new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                                PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                        # Record contact data if chain1 is a MHC
                        if (chain1+chain2) in complex_dict['dict_' + PDB_code]['MHC_ligand']:
                            atom2 = g9
                            ### Don't need to exclude peptide backbone atoms for 
                            ### MHC contacts, but could reinstate line below if 
                            ### want to exclude backbone atoms
                            #if atom2 not in backbone_list:
                            PDBsum_filetype = 'ligand'
                            chain1_type = 'MHC'
                            residue2 = g10
                            residue2_number = int(g11)
                            distance = g13
                            residue2_number_adjusted = residue2_number # In PDBsum ligand files the residue number does not appear to need adjusting i.e. they begin at 1
                            MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                            resolution = complex_dict['dict_' + PDB_code]['Resolution']
                            contact_type = 'Non-bonded contact'
                            new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                            PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                # Obtain information on salt-bridge contacts
                if isinstance(SB_start_line, int):
                    for n in range(SB_start_line, SB_end_line,1):
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        chain1 = g6
                        chain2 = g12
                        # Record contact data if chain1 is a TCR
                        if (chain1+chain2) in complex_dict['dict_' + PDB_code]['TCR_ligand']:
                            atom2 = g9
                            # only save contact data if atom is not in the
                            # backbone of amino acid in the peptide
                            if atom2 not in backbone_list:
                                PDBsum_filetype = 'ligand'
                                chain1_type = 'TCR'
                                residue2 = g10
                                residue2_number = int(g11)
                                distance = g13
                                residue2_number_adjusted = residue2_number # In PDBsum ligand files the residue number does not appear to need adjusting i.e. they begin at 1
                                MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                                resolution = complex_dict['dict_' + PDB_code]['Resolution']
                                contact_type = 'Salt-bridge contact'
                                new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                                PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                        # Record contact data if chain1 is a MHC
                        if (chain1+chain2) in complex_dict['dict_' + PDB_code]['MHC_ligand']:
                            atom2 = g9
                            ### Don't need to exclude peptide backbone atoms for 
                            ### MHC contacts, but could reinstate line below if 
                            ### want to exclude backbone atoms
                            #if atom2 not in backbone_list:
                            PDBsum_filetype = 'ligand'
                            chain1_type = 'MHC'
                            residue2 = g10
                            residue2_number = int(g11)
                            distance = g13
                            residue2_number_adjusted = residue2_number # In PDBsum ligand files the residue number does not appear to need adjusting i.e. they begin at 1
                            MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                            resolution = complex_dict['dict_' + PDB_code]['Resolution']
                            contact_type = 'Salt-bridge contact'
                            new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                            PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
        elif file[12:15] == 'TCR':
            with open(file,'r') as f:
                PDBsum_content = f.read().splitlines()
                # Define the start and end lines of the contact information 
                # for each type of contact
                for line in PDBsum_content:
                    count += 1
                    if HB.search(line):
                        HB_line_number = count
                        HB_start_line = HB_line_number + 7
                    elif NBC.search(line):
                        NBC_line_number = count
                        NBC_start_line = NBC_line_number + 7
                    elif SB.search(line):
                        SB_line_number = count
                        SB_start_line = SB_line_number + 7
                    elif HB_number.search(line):
                        HB_number_match = HB_number.search(line)
                        HB_end_line = HB_start_line + int(HB_number_match.group(1))
                    elif NBC_number.search(line):
                        NBC_number_match = NBC_number.search(line)
                        NBC_end_line = NBC_start_line + int(NBC_number_match.group(1))
                    elif SB_number.search(line):
                        SB_number_match = SB_number.search(line)
                        SB_end_line = SB_start_line + int(SB_number_match.group(1))
                
                # Obtain information on hydrogen bond contacts
                if isinstance(HB_start_line, int):
                    for n in range(HB_start_line, HB_end_line,1):
                        # Example contact data line from 1d9k hydrogen bond
                        #        Atom Atom Res  Res              Atom Atom Res  Res
                        #         no. name name no.  Chain        no. name name no.  Chain  Distance
                        # '  1.    211  O   THR   28    A   <-->  4833  NH1 ARG  135    P      3.12'
                        #    g1    g2   g3   g4   g5    g6   g7   g8    g9   g10  g11  g12      g13
                        
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        atom2 = g9
                        # only save contact data if atom is not in the 
                        # backbone of amino acid in the peptide
                        if atom2 not in backbone_list:
                            PDBsum_filetype = 'TCR-peptide'
                            chain1_type = 'TCR'
                            chain1 = g6
                            residue2 = g10
                            residue2_number = int(g11)
                            chain2 = g12
                            distance = g13
                            residue2_number_adjusted = ''
                            MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                            resolution = complex_dict['dict_' + PDB_code]['Resolution']
                            contact_type = 'Hydrogen bond'
                            new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                            
                            PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                # Obtain information on non-bonded contacts
                if isinstance(NBC_start_line, int):
                    for n in range(NBC_start_line, NBC_end_line,1):
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        atom2 = g9
                        # only save contact data if atom is not in the
                        # backbone of amino acid in the peptide
                        if atom2 not in backbone_list:
                            PDBsum_filetype = 'TCR-peptide'
                            chain1_type = 'TCR'
                            chain1 = g6
                            residue2 = g10
                            residue2_number = int(g11)
                            chain2 = g12
                            distance = g13
                            residue2_number_adjusted = ''
                            MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                            resolution = complex_dict['dict_' + PDB_code]['Resolution']
                            contact_type = 'Non-bonded contact'
                            new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                            PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                # Obtain information on salt-bridge contacts
                if isinstance(SB_start_line, int):
                    for n in range(SB_start_line, SB_end_line,1):
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        atom2 = g9
                        # only save contact data if atom is not in the 
                        # backbone of amino acid in the peptide
                        if atom2 not in backbone_list:
                            PDBsum_filetype = 'TCR-peptide'
                            chain1_type = 'TCR'
                            chain1 = g6
                            residue2 = g10
                            residue2_number = int(g11)
                            chain2 = g12
                            distance = g13
                            residue2_number_adjusted = ''
                            MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                            resolution = complex_dict['dict_' + PDB_code]['Resolution']
                            contact_type = 'Salt-bridge contact'
                            new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                            PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
        
        # Repeat but for MHC contact files
        elif file[12:15] == 'MHC':
            with open(file,'r') as f:
                PDBsum_content = f.read().splitlines()
                # Define the start and end lines of the contact information 
                # for each type of contact
                for line in PDBsum_content:
                    count += 1
                    if HB.search(line):
                        HB_line_number = count
                        HB_start_line = HB_line_number + 7
                    elif NBC.search(line):
                        NBC_line_number = count
                        NBC_start_line = NBC_line_number + 7
                    elif SB.search(line):
                        SB_line_number = count
                        SB_start_line = SB_line_number + 7
                    elif HB_number.search(line):
                        HB_number_match = HB_number.search(line)
                        HB_end_line = HB_start_line + int(HB_number_match.group(1))
                    elif NBC_number.search(line):
                        NBC_number_match = NBC_number.search(line)
                        NBC_end_line = NBC_start_line + int(NBC_number_match.group(1))
                    elif SB_number.search(line):
                        SB_number_match = SB_number.search(line)
                        SB_end_line = SB_start_line + int(SB_number_match.group(1))
                
                # Obtain information on hydrogen bond contacts
                if isinstance(HB_start_line, int):
                    for n in range(HB_start_line, HB_end_line,1):
                        # Example contact data line from 1d9k hydrogen bond
                        #        Atom Atom Res  Res              Atom Atom Res  Res
                        #         no. name name no.  Chain        no. name name no.  Chain  Distance
                        # '  1.    211  O   THR   28    A   <-->  4833  NH1 ARG  135    P      3.12'
                        #    g1    g2   g3   g4   g5    g6   g7   g8    g9   g10  g11  g12      g13
                        
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        atom2 = g9
                        ### Don't need to exclude peptide backbone atoms for 
                        ### MHC contacts, but could reinstate line below if 
                        ### want to exclude backbone atoms
                        #if atom2 not in backbone_list:
                        PDBsum_filetype = 'MHC-peptide'
                        chain1_type = 'MHC'
                        chain1 = g6
                        residue2 = g10
                        residue2_number = int(g11)
                        chain2 = g12
                        distance = g13
                        residue2_number_adjusted = ''
                        MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                        resolution = complex_dict['dict_' + PDB_code]['Resolution']
                        contact_type = 'Hydrogen bond'
                        new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                        PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                # Obtain information on non-bonded contacts
                if isinstance(NBC_start_line, int):
                    for n in range(NBC_start_line, NBC_end_line,1):
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        atom2 = g9
                        ### Don't need to exclude peptide backbone atoms for 
                        ### MHC contacts, but could reinstate line below if 
                        ### want to exclude backbone atoms
                        #if atom2 not in backbone_list:
                        PDBsum_filetype = 'MHC-peptide'
                        chain1_type = 'MHC'
                        chain1 = g6
                        residue2 = g10
                        residue2_number = int(g11)
                        chain2 = g12
                        distance = g13
                        residue2_number_adjusted = ''
                        MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                        resolution = complex_dict['dict_' + PDB_code]['Resolution']
                        contact_type = 'Non-bonded contact'
                        new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                        PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
                # Obtain information on salt-bridge contacts
                if isinstance(SB_start_line, int):
                    for n in range(SB_start_line, SB_end_line,1):
                        contact_data = PDBsum_content[n]
                        #for data in contact_data:
                        (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13) = contact_data.split()
                        atom2 = g9
                        ### Don't need to exclude peptide backbone atoms for 
                        ### MHC contacts, but could reinstate line below if 
                        ### want to exclude backbone atoms
                        #if atom2 not in backbone_list:
                        PDBsum_filetype = 'MHC-peptide'
                        chain1_type = 'MHC'
                        chain1 = g6
                        residue2 = g10
                        residue2_number = int(g11)
                        chain2 = g12
                        distance = g13
                        residue2_number_adjusted = ''
                        MHC_class = complex_dict['dict_' + PDB_code]['MHC_class']
                        resolution = complex_dict['dict_' + PDB_code]['Resolution']
                        contact_type = 'Salt-bridge contact'
                        new_contact_data = [PDB_code,PDBsum_filetype,resolution,chain1_type,chain1,chain2,atom2,residue2,residue2_number,residue2_number_adjusted,distance,MHC_class,contact_type]
                        PDBsum_contact_data = np.append(PDBsum_contact_data,[new_contact_data])
    
    PDBsum_contact_data_reshaped = PDBsum_contact_data.reshape(int(len(PDBsum_contact_data)/13),13)
    
    
    
    #change back to home directory
    os.chdir(home_path)
    
    #for n in range(0,len(PDBsum_contact_data_reshaped),1): 
    #    print(n,PDBsum_contact_data_reshaped[n,10])
    #    if PDBsum_contact_data_reshaped[n,10] == 'Salt-bridge contact': 
    #        print(PDBsum_contact_data_reshaped[n,0],PDBsum_contact_data_reshaped[n,10])
    
    import pandas as pd
    contact_df = pd.DataFrame({'PDB code': PDBsum_contact_data_reshaped[:, 0],
                            'PDBsum file type': PDBsum_contact_data_reshaped[:, 1],
                            'Resolution': PDBsum_contact_data_reshaped[:, 2],
                            'Chain 1 type': PDBsum_contact_data_reshaped[:, 3],
                            'Chain 1 (TCR or MHC)': PDBsum_contact_data_reshaped[:, 4],
                            'Chain 2 (peptide/ligand)': PDBsum_contact_data_reshaped[:, 5],
                            'Atom 2': PDBsum_contact_data_reshaped[:, 6],
                            'Residue 2': PDBsum_contact_data_reshaped[:, 7],
                            'Residue no.': PDBsum_contact_data_reshaped[:, 8],
                            'Residue no. adjusted': PDBsum_contact_data_reshaped[:, 9],
                            'Contact distance': PDBsum_contact_data_reshaped[:, 10],
                            'MHC class': PDBsum_contact_data_reshaped[:, 11],
                            'Contact type': PDBsum_contact_data_reshaped[:, 12],})
    
    #remove initial headings used in array
    contact_df = contact_df.drop(contact_df.index[[0]])
    
    # https://stackoverflow.com/questions/17071871/how-do-i-select-rows-from-a-dataframe-based-on-column-values
                    
    return contact_df

def add_residue_start(home_path,complex_dict,PDB_file_dir,contact_df):
    """
    Extracts peptide starting position within complex structure and adds to 
        dataframe so the residue positions within the peptide can be calulated
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    :complex_dict: Dictionary containing various information on the TCR-pMHC
        complexes organised with a sub-dictionary for each PDB code/complex
    :PDB_file_dir: Directory containing files downloaded from PDB
    :contact_df: Dataframe containing information on peptide residue contacts
        with TCRs (excludes backbone atoms)
    """
    import os
    import pandas as pd
    import numpy as np
    import re
    
    os.chdir(PDB_file_dir)
    
    # Create unique list of PDB codes in the contact dataframe
    # Ignore peptides treated as ligands in PDB sum as the residue number is 
    # already aligned to starting residue
    peptide_only_df = contact_df.loc[contact_df['PDBsum file type'] != 'ligand']
    res_start_PDB_code_list = pd.unique(peptide_only_df['PDB code'])
    
    # Match peptide start residue number
    #s = re.compile(r'COMPND   2 MOLECULE: HLA CLASS II HISTOCOMPATIBILITY ANTIGEN, DR ALPHA CHAIN;')
    residue_start_line = re.compile(r'^DBREF\s+')
    
    # convert column in data frame to integers so can use in calculations later
    contact_df['Residue no.'] = contact_df['Residue no.'].astype(int)
    
    # Cycle through the list of PDB codes for TCR-pMHC complexes
    for n in range(0,len(res_start_PDB_code_list),1): 
        # read PDB files
        PDB_code = res_start_PDB_code_list[n]
        filename = PDB_code + '.pdb'
        m=len(res_start_PDB_code_list)
        #print(filename + ": file " + str(n+1) + " of " + str(m))
        with open(filename,'r') as f:
            PDB_content = f.read().splitlines()    
        for line in PDB_content:
            if residue_start_line.search(line):
                #print(line)
                split_line = line.split()
                chain = split_line[2]
                residue_start_number = int(split_line[3])
                
                # Replace 'Residue no. adjusted' with real values 
                # i.e. 'Residue no.' from PDBsum files minus the starting 
                # residue number in the PDB files                
                contact_df['Residue no. adjusted'] = np.where(((contact_df['PDB code'] == PDB_code) & (contact_df['Chain 2 (peptide/ligand)'] == chain)),contact_df['Residue no.'] - residue_start_number + 1,contact_df['Residue no. adjusted'])
       
    os.chdir(home_path)
    
    # convert columns in data frame to integers so can use in calculations later
    contact_df['Residue no.'] = contact_df['Residue no.'].astype(int)
    contact_df['Residue no. adjusted'] = contact_df['Residue no. adjusted'].astype(int)
    contact_df['Resolution'] = contact_df['Resolution'].astype(float)
    contact_df['Contact distance'] = contact_df['Contact distance'].astype(float)
    
    return contact_df

def add_pep_seq(PDB_code_list,home_path,complex_dict,PDB_file_dir,contact_df):
    """
    Adds peptide sequence and length to the contact_df dataframe
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    :complex_dict: Dictionary containing various information on the TCR-pMHC complexes 
        organised with a sub-dictionary for each PDB code/complex
    :PDB_file_dir: Directory containing files downloaded from PDB
    :contact_df: Dataframe containing information on peptide residue contacts with TCRs (excludes backbone atoms)
    """
    import os
    import pandas as pd
    import numpy as np
    import re
    
    # set PDB file directory as current directory
    os.chdir(PDB_file_dir)
    # POTENTIAL IMPROVEMENT: COULD ADD CODE TO SCAN DIRECTORY AND EXTRACT FROM 
    # ALL FILES IN DIR RATHER THAN JUST FROM THE LIST - IN CASE SOME FILES 
    # DON'T DOWNLOAD FOR SOME REASON
    
    # Find sequence length ('16' in lines below) and amino acid sequence in PDB file e.g.
    # SEQRES   1 P   16  GLY ASN SER HIS ARG GLY ALA ILE GLU TRP GLU GLY ILE          
    # SEQRES   2 P   16  GLU SER GLY                                         
    seqres = re.compile(r'SEQRES\s+?\d+\s+(\w+)\s+(\d+)\s+(.+)')
    rows = []
    # PDB_code_list_test = ['1d9k']
    
    #for n in range(0,len(PDB_code_list_test),1): # Cycle through the list of PDB codes for TCR-pMHC complexes
    for n in range(0,len(PDB_code_list),1): # Cycle through the list of PDB codes for TCR-pMHC complexes
        # read PDB files
        PDB_code = PDB_code_list[n]
        filename = PDB_code + '.pdb'
        m=len(PDB_code_list)
        #print(filename + ": file " + str(n+1) + " of " + str(m))
        with open(filename,'r') as f:
            PDB_content = f.read().splitlines()    
    
        peptide_chain = []
        if len(complex_dict['dict_' + PDB_code]['peptide']) == 0:
            peptide_chain = complex_dict['dict_' + PDB_code]['ligand']
        else:
            peptide_chain = complex_dict['dict_' + PDB_code]['peptide']
        #print(PDB_code,peptide_chain)
        
        two_line_count = 0
        for c in range(0,len(PDB_content)-5,1):
            # Define two lines in the PDB content as peptide sequence can 
            # sometimes be on two lines
            # Also need to skip a line in PDB_content if a peptide sequence 
            # is found to be on two lines
            if two_line_count != 0:
                PDB_content_line_1 = PDB_content[c + two_line_count]
                PDB_content_line_2 = PDB_content[c + two_line_count + 1]
            else:    
                PDB_content_line_1 = PDB_content[c]
                PDB_content_line_2 = PDB_content[c+1]
            # (Re)initiate variables
            aa_length = -1
            pep_seq = ''
            chain = ''
            chain_2 = ''
    
    
            # Find row with sequence data in 
            if seqres.search(PDB_content_line_1):
                #print(c,PDB_content_line_1)
                match_seqres = seqres.search(PDB_content_line_1)
                if match_seqres.group(1) in peptide_chain: # check if chain relates to a peptide
                    chain = match_seqres.group(1)
                    aa_length = match_seqres.group(2)
                    # Determine whether peptide sequence is taken from 1 or 2 lines in PDB file
                    if (int(aa_length) > len(match_seqres.group(3).split())):
                        match_seqres_2 = seqres.search(PDB_content_line_2)
                        chain_2 = match_seqres_2.group(1)
                        if chain_2 == chain: # determine if chains on line 1 and line 2 match
                            two_line_count = 1 + two_line_count # When a sequence is on 2 lines need the next search to begin on the 3rd line
                            pep_seq = match_seqres.group(3).split() + match_seqres_2.group(3).split()
                            rows.append([PDB_code,chain,aa_length,pep_seq])
                    else:
                        pep_seq = match_seqres.group(3).split()
                        rows.append([PDB_code,chain,aa_length,pep_seq])
    # Create dataframe with rows of data and column headings                                        
    pep_seq_df = pd.DataFrame(rows,columns=['PDB code','Chain 2 (peptide/ligand)','Peptide chain length','Peptide sequence'])
    
    # Merge dataframes
    contact_df_final = pd.merge(contact_df,pep_seq_df,on=['PDB code', 'Chain 2 (peptide/ligand)'],how="left")
    ## Check data types in columns
    #contact_df_final.dtypes
    ## Convert columns in data frame to integers so can use in calculations later
    contact_df_final['MHC class'] = contact_df_final['MHC class'].astype(int)
    contact_df_final['Peptide chain length'] = contact_df_final['Peptide chain length'].astype(int)
    
    # Switch back to home directory    
    os.chdir(home_path)
    
    return contact_df_final, pep_seq_df

def add_sequence(PDB_code_list,home_path,complex_dict,PDB_file_dir,dict_aa_3to1,dict_aa_1to3):
    """
    Adds protein sequences to complex_dict in the Chains_info_2 entry
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :home_path: Directory so can switch between main directory and 
        sub-directory containing PDB Sum files
    :complex_dict: Dictionary containing various information on the TCR-pMHC complexes 
        organised with a sub-dictionary for each PDB code/complex
    :PDB_file_dir: Directory containing files downloaded from PDB
    :dict_aa_3to1: Dictionary of amino acids. 3 letter codes as keys and 1 letter codes as values 
    :dict_aa_1to3: Dictionary of amino acids. 1 letter codes as keys and 3 letter codes as values 
    :non_standard_list: List of non-standard amino acids found in chains
    """
    import os
    import re
    
    # set PDB file directory as current directory
    os.chdir(PDB_file_dir)
    # POTENTIAL IMPROVEMENT: COULD ADD CODE TO SCAN DIRECTORY AND EXTRACT FROM 
    # ALL FILES IN DIR RATHER THAN JUST FROM THE LIST - IN CASE SOME FILES 
    # DON'T DOWNLOAD FOR SOME REASON
    
    # Find sequence length ('16' in lines below) and amino acid sequence in PDB file e.g.
    # SEQRES   1 P   16  GLY ASN SER HIS ARG GLY ALA ILE GLU TRP GLU GLY ILE          
    # SEQRES   2 P   16  GLU SER GLY                                         
    seqres = re.compile(r'SEQRES\s+?\d+\s+(\w+)\s+(\d+)\s+(.+)')
    #PDB_code_list_test = ['1ao7']
    
    non_standard_list = []
    #for n in range(0,len(PDB_code_list_test),1): # Cycle through the list of PDB codes for TCR-pMHC complexes
    for n in range(0,len(PDB_code_list),1): # Cycle through the list of PDB codes for TCR-pMHC complexes
        # read PDB files
        PDB_code = PDB_code_list[n]
        filename = PDB_code + '.pdb'
        m=len(PDB_code_list)
        #print(filename + ": file " + str(n+1) + " of " + str(m))
        
        with open(filename,'r') as f:
            PDB_content = f.read().splitlines()    
    
        next_line = 0
        for c in range(0,len(PDB_content),1):
            if c < next_line:
                continue
            else:    
                # (re)initiate sequence
                sequence = ''
                
                # Find row with sequence data in 
                if seqres.search(PDB_content[c]):
                    match_seqres = seqres.search(PDB_content[c])
                    chain_in_PDB = match_seqres.group(1)
                    chain_length = match_seqres.group(2)
                    # Determine number of lines on which the sequence will be
                    # present in PDB file (round up)
                    sequence_lines = int((int(chain_length) / 13) + (int(chain_length) % 13 > 0)) # Second term gives a true (1) or false (0) so adds 1 if there is a remainder and adds 0 if there is no remainder which has the effect of rounding up
                    for n in range(0,sequence_lines,1):
                        match_seqres = seqres.search(PDB_content[c+n])
                        sequence = sequence + match_seqres.group(3) #.split()
                    # Format sequence as list
                    seq_3_letter = sequence.split()
                    
                    # Create sequence in 1 letter format
                    seq_1_letter = []
                    for m in range(0,len(seq_3_letter),1):
                        # There are some non-standard amino acids in 
                        # some PDB files and so this checks
                        if seq_3_letter[m] in dict_aa_3to1:
                            seq_1_letter.append(dict_aa_3to1[seq_3_letter[m]])
                        elif seq_3_letter[m] == 'PFF': 
                            seq_1_letter.append('F') # standard parent of PFF is PHE https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/PFF
                            non_standard_list.append('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                        elif seq_3_letter[m] == 'F2F': 
                            seq_1_letter.append('F') # standard parent of F2F is PHE https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/F2F
                            non_standard_list.append('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                        elif seq_3_letter[m] == 'UNK': 
                            seq_1_letter.append('U')
                            non_standard_list.append('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                        elif seq_3_letter[m] == 'PCA': 
                            seq_1_letter.append('Q') # standard parent of PCA is GLN https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/PCA
                            non_standard_list.append('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                        elif seq_3_letter[m] == 'CIR': # standard parent of CIR is ARG https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/CIR
                            seq_1_letter.append('R')
                            non_standard_list.append('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                        else:
                            seq_1_letter.append('9') # Using 9 as an error code within the one letter sequence
                            print('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                            non_standard_list.append('Error in ' + PDB_code + ' chain ' + chain_in_PDB + ', amino acid number '+ str(m+1) + ' - amino acid ' + seq_3_letter[m] + ' not in dictionary')
                    
                    # Save 3 and 1 letter sequences in dictionary
                    # Search through dictionary entry to find correct chain
                    for i in range(0,len(complex_dict['dict_' + PDB_code]['Chains_info_2'][0]),1):
                        # Find the correct chain
                        if complex_dict['dict_' + PDB_code]['Chains_info_2'][0][i][0] == chain_in_PDB:
                            # Join tuple and list as tuple
                            new_chains_info_2 = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][i] + (seq_3_letter,seq_1_letter)
                            # Replace old tuple with new tuple that includes chain sequence
                            complex_dict['dict_' + PDB_code]['Chains_info_2'][0][i] = new_chains_info_2
                                                
                    # Need to move to line where next chain starts
                    next_line = c + sequence_lines
        
    # Switch back to home directory    
    os.chdir(home_path)

    return complex_dict, non_standard_list

def make_logo_plot(df,logo_option,chain_option):
    """
    Creates a probability matrix for a particular amino acid at a particular
    position in the peptide sequence
    :df: Dataframe to create the logo plot from - contains information on 
        peptide residue contacts with TCRs (excludes backbone atoms)
    :logo_option: Option to determine whether logo plot shows the number 
        ('number') of amino acids or the probability ('probability')
    :chain_option: Option to determine whether all peptide chains ('all') are used in 
        contact data or just one per PDB code/TCR-pMHC complex ('unique')
    """
    import pandas as pd
    import logomaker
    import matplotlib.pyplot as plt
    
    # Building a probability matrix for logomaker
    
    # inputs dataframe
    
    #df = contact_df_final
    #df = df_MHC1_p9
    
    pep_seq_list = []
    PDB_code_chain_list = []
    
    # Build list of unique PDB codes linked to chains and peptide sequences
    if chain_option == 'all':
        for index in df.index: # the indexing in the df changes when filtered (i.e. it will have gaps in the sequence of numbers) so a for n in range() doesn't work
            PDB_code_chain = df.loc[index][0]+df.loc[index][5]
            if PDB_code_chain not in PDB_code_chain_list:
                # Add PDB code and chain to unique list of PDB codes
                PDB_code_chain_list.append(PDB_code_chain)
                # Add peptide sequence to list
                pep_seq_list.append(df.loc[index][14])
    elif chain_option == 'unique':
        for index in df.index: # the indexing in the df changes when filtered (i.e. it will have gaps in the sequence of numbers) so a for n in range() doesn't work
            PDB_code_chain = df.loc[index][0]
            #print(index,PDB_code_chain)
            if PDB_code_chain not in PDB_code_chain_list:
                # add PDB code and chain to unique list of PDB codes
                PDB_code_chain_list.append(PDB_code_chain)
                # add peptide sequence to list
                pep_seq_list.append(df.loc[index][14])
    
    # Create dataframe with peptide residue positions as column headings and
    # peptide sequences in the rows
    max_pep_seq_length = df['Peptide chain length'].max()
    pep_seq_position_df = pd.DataFrame(pep_seq_list,columns=list(range(1,max_pep_seq_length+1,1)))
    # Need to shift 8-mer peptides along by one so that they align with 
    # 9-mer peptides
    for index in pep_seq_position_df.index:
        if pep_seq_position_df.loc[index][9] == None:
            pep_seq_position_df.loc[index]=pep_seq_position_df.loc[index].shift(1)
    
    # Amino acid codes (1-letter and 3-letter code lists are aligned)
    aa_list =['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z']
    aa_3_list = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','ASX','GLX']
    # Initiate lists
    total_aa_position_count_list = []
    number_aa_in_position_list = []
    number_aa_in_position_by_column = []
    
    # Cycle through columns of dataframe...
    for i in range(1,len(pep_seq_position_df.columns)+1,1):
        number_aa_in_position_list = []
        total_aa_position_count = 0
    
        #...and then cycle through amino acid list counting the number for 
        # each type of amino acid in the column
        for j in range(0,len(aa_3_list),1):
            count_aa_in_position = 0
            # Check if amino acid exists in the list of counts for the column in the dataframe
            if aa_3_list[j] in pep_seq_position_df[i].value_counts(dropna=False):
                count_aa_in_position = pep_seq_position_df[i].value_counts(dropna=False)[aa_3_list[j]]
            else:
                count_aa_in_position = 0
            total_aa_position_count += count_aa_in_position
            number_aa_in_position_list.append(count_aa_in_position) # Builds a list of the total number of each amino acid ordered as per the aa_3_list
        number_aa_in_position_by_column.append(number_aa_in_position_list) # Puts the lists into a list ordered by column number i.e. residue position number
        total_aa_position_count_list.append(total_aa_position_count) # List of the total number of amino acids at each residue position number
        
    # Divide totals for individual amino acids in peptide residue positions by the 
    # total number of amino acids in that peptide residue position in the data set 
    # i.e. the total number of peptides with of length equal to or greater than 
    # the peptide position
    aa_position_prob_matrix = []
    # Cycle through the list of total amino acids in each residue position...
    for m in range(0,len(total_aa_position_count_list),1):
        aa_prob_by_position = []
        #...and then cycle through tatals for the individual amino acids at each position...
        for n in range(0,len(aa_3_list),1):
            #...and divide to get the probabilty of having a particular amino acid
            prob = number_aa_in_position_by_column[m][n] / total_aa_position_count_list[m]
            aa_prob_by_position.append(prob)
        aa_position_prob_matrix.append(aa_prob_by_position)
    
    
    # Create dataframe with peptide residue positions as column headings and 
    # peptide sequences in the rows
    #pep_seq_position_df = pd.DataFrame(pep_seq_list,columns=list(range(1,24,1)))
    #prob_matrix_df = pd.DataFrame(list(range(1,23,1)),columns=aa_list)
    #logo_df = pd.DataFrame(aa_position_prob_matrix,columns=aa_list)
    #logo_df = pd.DataFrame(number_aa_in_position_by_column,columns=aa_list)
    
    if logo_option == 'number':
        logo_df = pd.DataFrame(number_aa_in_position_by_column,columns=aa_list)
    elif logo_option == 'probability':
        logo_df = pd.DataFrame(aa_position_prob_matrix,columns=aa_list)
    else:
        logo_df = pd.DataFrame(aa_position_prob_matrix,columns=aa_list)
        
    # create Logo object
    ss_logo = logomaker.Logo(logo_df,
                             width=.8,
                             vpad=.05,
                             fade_probabilities=False,
                             stack_order='small_on_top',
                             color_scheme='weblogo_protein',#dodgerblue',
                             font_name='Rosewood Std')
    
    # Style using Logo methods
    ss_logo.style_spines(spines=['left', 'right'], visible=False)
    
    # Style using Axes methods
    ss_logo.ax.set_xticks(range(len(logo_df)))
    #ss_logo.ax.set_xticklabels('%+d'%x for x in range(1,23,1))#[-3, -2, -1, 1, 2, 3, 4, 5, 6])
    ss_logo.ax.set_xticklabels(range(1,23+1,1))
    ss_logo.ax.set_yticks(range(0,1,0.2)) #[0, .2, 1])
    ss_logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')
    ss_logo.ax.set_ylabel('probability')
    

    return logo_image

def analyse_contacts(df):
    """
    Creates summary of contacts by residue position for each TCR-pMHC in a dataframe
    :df: Dataframe containing information on peptide residue contacts with columns
        Index, PDB code, PDBsum file type, Resolution, Chain 1 type, Chain 1 (TCR or MHC),
        Chain 2 (peptide/ligand), Atom 2, Residue 2, Residue no., Residue no. adjusted,
        Contact distance, MHC class, Contact type, Peptide chain length, Peptide sequence
    """
    import pandas as pd
    import numpy as np
    
    # Check what the maximum peptide chain length is so can adjust calculations and array
    max_peptide_chain_length = np.max(df['Peptide chain length'])
    
    # Initiate column headings for output dataframe
    #summary_contact_data = np.array(['PDB code','1','2','3','4','5','6','7','8','9'])
    summary_contact_data = np.array(['PDB code'])
    for i in range(1,max_peptide_chain_length+1,1):
        summary_contact_data = np.append(summary_contact_data,i)
    
    # Gather unique PDB codes from input dataframe
    PDB_code_list = pd.unique(df['PDB code'])
    
    # Initiate lists
    contact_list = []
        
    # Cycle through PDB codes... 
    for PDB_code in PDB_code_list:
        #...and peptide residue position numbers counting the number of contacts
        # (i.e the number of rows in the input dataframe)
        contact_list = [PDB_code]
        for n in range(1,max_peptide_chain_length+1,1):
            number_of_contacts = len(df.loc[(df['PDB code'] == PDB_code) & (df['Residue no. adjusted'] == n)])
            contact_list.append(number_of_contacts)
        # Add PDB code and list of the number of contacts at each residue position
        summary_contact_data = np.append(summary_contact_data,[contact_list])
    
    # Reshape array
    summary_contact_data_reshaped = summary_contact_data.reshape(int(len(summary_contact_data)/(max_peptide_chain_length+1)),(max_peptide_chain_length+1))
    
    # Convert array to dataframe
    summary_contact_df = pd.DataFrame({'PDB code': summary_contact_data_reshaped[:, 0],
                        '1': summary_contact_data_reshaped[:, 1],
                        '2': summary_contact_data_reshaped[:, 2],
                        '3': summary_contact_data_reshaped[:, 3],
                        '4': summary_contact_data_reshaped[:, 4],
                        '5': summary_contact_data_reshaped[:, 5],
                        '6': summary_contact_data_reshaped[:, 6],
                        '7': summary_contact_data_reshaped[:, 7],
                        '8': summary_contact_data_reshaped[:, 8],
                        '9': summary_contact_data_reshaped[:, 9],})
                        #'10': summary_contact_data_reshaped[:, 10],})
    
    # Convert columns in data frame to integers so can use in calculations later
    summary_contact_df['1'] = summary_contact_df['1'].astype(int)
    summary_contact_df['2'] = summary_contact_df['2'].astype(int)
    summary_contact_df['3'] = summary_contact_df['3'].astype(int)
    summary_contact_df['4'] = summary_contact_df['4'].astype(int)
    summary_contact_df['5'] = summary_contact_df['5'].astype(int)
    summary_contact_df['6'] = summary_contact_df['6'].astype(int)
    summary_contact_df['7'] = summary_contact_df['7'].astype(int)
    summary_contact_df['8'] = summary_contact_df['8'].astype(int)
    summary_contact_df['9'] = summary_contact_df['9'].astype(int)
    #summary_contact_df['10'] = summary_contact_df['10'].astype(int)
        
    # Remove initial headings used in array
    summary_contact_df = summary_contact_df.drop(summary_contact_df.index[[0]])

    summary_contact_stats = summary_contact_df.describe()
    
    return summary_contact_df,summary_contact_stats

def identical_chains(PDB_code_list,complex_dict):
    """
    Filters set of TCR-pMHC complexes in complex_dict and returns chains that 
        occur more than once and split by molecule type (MHC,TCR and peptide)
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :complex_dict: Dictionary containing various information on the TCR-pMHC complexes 
        organised with a sub-dictionary for each PDB code/complex
    """
    
    import pandas as pd
    
    # Build rows of data for each type of molecule
        
    # Initiate lists
    rows_peptide = []
    rows_TCR = []
    rows_MHC = []
    
    # Cycle through all complexes
    for PDB_code in PDB_code_list:
        # Cycle through all chains in each complex
        for n in range(0,len(complex_dict['dict_' + PDB_code]['Chains_info_2'][0]),1):
            molecule_type = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][3]
            molecule_name = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][2]
            chain = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][0]
            seq_length = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][1]
            sequence = ''.join(complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][5])
            
            if molecule_type == 'peptide':
                rows_peptide.append([PDB_code,chain,seq_length,sequence,molecule_name])
            elif molecule_type == 'ligand':
                rows_peptide.append([PDB_code,chain,seq_length,sequence,molecule_name])
            elif molecule_type == 'MHC':
                rows_MHC.append([PDB_code,chain,seq_length,sequence,molecule_name])
            elif molecule_type == 'TCR':
                rows_TCR.append([PDB_code,chain,seq_length,sequence,molecule_name])
            else:
                continue # don't collect data for the 'other' category of molecule
    
    # Create dataframes
    df_seq_peptide = pd.DataFrame(rows_peptide,columns=['PDB code','Chain','Chain length','Sequence','Molecule name'])
    df_seq_TCR = pd.DataFrame(rows_TCR,columns=['PDB code','Chain','Chain length','Sequence','Molecule name'])
    df_seq_MHC = pd.DataFrame(rows_MHC,columns=['PDB code','Chain','Chain length','Sequence','Molecule name'])
    
    # Change formatting so can filter dataframes for chains appearing multiple times
    list_seq_peptide = list(df_seq_peptide['Sequence'])
    list_seq_peptide_strings = []
    list_seq_TCR = list(df_seq_TCR['Sequence'])
    list_seq_TCR_strings = []
    list_seq_MHC = list(df_seq_MHC['Sequence'])
    list_seq_MHC_strings = []
    
    # Convert list of lists to list of strings so can compare easily
    for item in list_seq_peptide:
        list_seq_peptide_strings.append(str(item))
    for item in list_seq_TCR:
        list_seq_TCR_strings.append(str(item))
    for item in list_seq_MHC:
        list_seq_MHC_strings.append(str(item))
    
    # Searches list of sequences as strings and returns the sequences which
    # are present more than once
    set_seq_peptide = set([i for i in list_seq_peptide_strings if list_seq_peptide_strings.count(i)>1])
    set_seq_TCR = set([i for i in list_seq_TCR_strings if list_seq_TCR_strings.count(i)>1])
    set_seq_MHC = set([i for i in list_seq_MHC_strings if list_seq_MHC_strings.count(i)>1])
    
    # Convert from set to list
    set_list_seq_peptide = list(set_seq_peptide)
    set_list_seq_TCR = list(set_seq_TCR)
    set_list_seq_MHC = list(set_seq_MHC)
    
    #Test
    #a = df_seq_TCR['Sequence'][0]
    #list_seq_TCR_strings.count(a)
    #for n in range(0,len(list_seq_TCR_strings),1):
    #    if list_seq_TCR_strings[n] == a:
    #        print(n,df_seq_TCR.iloc[n,0],df_seq_TCR.iloc[n,1],df_seq_TCR.iloc[n,2],df_seq_TCR.iloc[n,3])
    
    # Cycle through set of peptide sequences that occur more than once
    row = []
    identical_chains_peptide = []
    for n in range(0,len(set_list_seq_peptide),1):
        # Cycle through full set of TCR sequences looking for a match    
        for m in range(0,len(df_seq_peptide),1):
            if str(df_seq_peptide.iloc[m,3]) == set_list_seq_peptide[n]:
                row = [m,df_seq_peptide.iloc[m,0],df_seq_peptide.iloc[m,1],df_seq_peptide.iloc[m,2],df_seq_peptide.iloc[m,3],df_seq_peptide.iloc[m,4]]
                #print(m,df_seq_peptide.iloc[m,0],df_seq_peptide.iloc[m,1],df_seq_peptide.iloc[m,2],df_seq_peptide.iloc[m,3],df_seq_peptide.iloc[m,4])
                identical_chains_peptide.append(row)
    
    # Cycle through set of TCR sequences that occur more than once
    identical_chains_TCR = []
    for n in range(0,len(set_list_seq_TCR),1):
        # Cycle through full set of TCR sequences looking for a match    
        for m in range(0,len(df_seq_TCR),1):
            if str(df_seq_TCR.iloc[m,3]) == set_list_seq_TCR[n]:
                row = [m,df_seq_TCR.iloc[m,0],df_seq_TCR.iloc[m,1],df_seq_TCR.iloc[m,2],df_seq_TCR.iloc[m,3],df_seq_TCR.iloc[m,4]]
                #print(m,df_seq_TCR.iloc[m,0],df_seq_TCR.iloc[m,1],df_seq_TCR.iloc[m,2],df_seq_TCR.iloc[m,3],df_seq_TCR.iloc[m,4])
                identical_chains_TCR.append(row)
                    
    # Cycle through set of MHC sequences that occur more than once
    identical_chains_MHC = []
    for n in range(0,len(set_list_seq_MHC),1):
        # Cycle through full set of MHC sequences looking for a match    
        for m in range(0,len(df_seq_MHC),1):
            if str(df_seq_MHC.iloc[m,3]) == set_list_seq_MHC[n]:
                row = [m,df_seq_MHC.iloc[m,0],df_seq_MHC.iloc[m,1],df_seq_MHC.iloc[m,2],df_seq_MHC.iloc[m,3],df_seq_MHC.iloc[m,4]]
                #print(m,df_seq_MHC.iloc[m,0],df_seq_MHC.iloc[m,1],df_seq_MHC.iloc[m,2],df_seq_MHC.iloc[m,3],df_seq_MHC.iloc[m,4])
                identical_chains_MHC.append(row)
    
    # Write lists to files of chains occuring more than once with a new file
    # for each molecule type
    import csv
         
    filenames = ['identical_chains_peptide.csv','identical_chains_TCR.csv','identical_chains_MHC.csv']
    lists_identical_chains = [identical_chains_peptide,identical_chains_TCR,identical_chains_MHC]
    
    for n in range(0,len(filenames),1):
        with open(filenames[n], 'w') as f:
            for item in lists_identical_chains[n]:
                f.write("%s\n" % item)
        print(filenames[n] + ' created')
    
    return df_seq_peptide, df_seq_TCR, df_seq_MHC

def identical_chains_matrix(PDB_code_list,complex_dict):
    """
    Filters set of TCR-pMHC complexes in complex_dict and returns chains that 
        occur more than once and split by molecule type (MHC,TCR and peptide)
    :PDB_code_list: List of PDB codes for the TCR-pMHC complexes
    :complex_dict: Dictionary containing various information on the TCR-pMHC complexes 
        organised with a sub-dictionary for each PDB code/complex
    """
                    
    import pandas as pd
    
    # Build rows of data for each type of molecule
    check_pep_list = []
    check_MHC_list = []
    check_TCR_list = []
    check_pep_rows = []
    check_MHC_rows = []
    check_TCR_rows = []
        
    #PDB_code_list_test = ['2esv']
    # Time function
    import time
    t0 = time.time()
    t1 = 0
    PDB_count = 0
    number_of_PDBs = len(PDB_code_list)
    PDB_2_counter = 0
 
    # Cycle through all complexes
    for PDB_code_1 in PDB_code_list:
        PDB_count += 1
        t2 = t1
        t1 = time.time()
        total = t1-t0
        if PDB_count > 1:
            time_estimate = (total / (PDB_count -1)) * (number_of_PDBs - PDB_count)
            print("PDB code ",PDB_code_1,"number",PDB_count," of ",number_of_PDBs," in ",round(t1-t2,1)," seconds. Total time taken: ",round(total,1)," seconds. Estimated time to completion: ",round(time_estimate,1)," seconds")
        else:
            print("PDB code ",PDB_code_1," number ",PDB_count," of ",number_of_PDBs)
            
        # (Re)initiate lists for next PDB code / complex
        rows_peptide_PDB_1 = []
        rows_TCR_PDB_1 = []
        rows_MHC_PDB_1 = []
        
        check_pep_list = []
        check_TCR_list = []
        check_MHC_list = []
         
        # Cycle through all chains in each complex
        for n in range(0,len(complex_dict['dict_' + PDB_code_1]['Chains_info_2'][0]),1):
            molecule_type_PDB_1 = complex_dict['dict_' + PDB_code_1]['Chains_info_2'][0][n][3]
            #molecule_name = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][2]
            #chain = complex_dict['dict_' + PDB_code]['Chains_info_2'][0][n][0]
            seq_length_PDB_1 = complex_dict['dict_' + PDB_code_1]['Chains_info_2'][0][n][1]
            sequence_PDB_1 = ''.join(complex_dict['dict_' + PDB_code_1]['Chains_info_2'][0][n][5])
            
            if molecule_type_PDB_1 == 'peptide':
                rows_peptide_PDB_1.append([seq_length_PDB_1,sequence_PDB_1])
            elif molecule_type_PDB_1 == 'ligand':
                rows_peptide_PDB_1.append([seq_length_PDB_1,sequence_PDB_1])
            elif molecule_type_PDB_1 == 'MHC':
                if seq_length_PDB_1 == '100': # ignore matches of Beta-2 Microglobulin which has sequence length 100 as there will be many matches
                    continue
                rows_MHC_PDB_1.append([seq_length_PDB_1,sequence_PDB_1])
            elif molecule_type_PDB_1 == 'TCR':
                rows_TCR_PDB_1.append([seq_length_PDB_1,sequence_PDB_1])
            else:
                continue # don't collect data for the 'other' category of molecule
            
        # Cycle through  all complexes and check if any of the chains match
        for PDB_code_2 in PDB_code_list:
            #PDB_2_counter += 1
            #print(PDB_2_counter)
            # Cycle through all chains in each complex
            # Only need to record a single match so create a switch so bypass molecule types once a match already found
            pep_state = 0
            TCR_state = 0
            MHC_state = 0
            
            for m in range(0,len(complex_dict['dict_' + PDB_code_2]['Chains_info_2'][0]),1):
                molecule_type_PDB_2 = complex_dict['dict_' + PDB_code_2]['Chains_info_2'][0][m][3]
                seq_length_PDB_2 = complex_dict['dict_' + PDB_code_2]['Chains_info_2'][0][m][1]
                sequence_PDB_2 = ''.join(complex_dict['dict_' + PDB_code_2]['Chains_info_2'][0][m][5])
                
                # Check molecule in PDB 2 against the list of molecules created for PDB 1
                if molecule_type_PDB_2 == 'peptide':
                    if pep_state == 1: # If already found a matching peptide then don't look at any more
                        continue
                    for peptide_PDB_1 in rows_peptide_PDB_1:
                        if pep_state == 1: 
                            continue
                        if seq_length_PDB_2 == peptide_PDB_1[0]:
                            if sequence_PDB_2 == peptide_PDB_1[1]:
                                pep_state = 1
                                check_pep_list.append(1)
                
                elif molecule_type_PDB_2 == 'ligand':
                    if pep_state == 1: # If already found a matching peptide then don't look at any more
                        continue
                    for peptide_PDB_1 in rows_peptide_PDB_1:
                        if pep_state == 1: 
                            continue
                        if seq_length_PDB_2 == peptide_PDB_1[0]:
                            if sequence_PDB_2 == peptide_PDB_1[1]:
                                pep_state = 1
                                check_pep_list.append(1)
                        
                elif molecule_type_PDB_2 == 'MHC':
                    if seq_length_PDB_2 == '100': # Ignore matches of Beta-2 Microglobulin which has sequence length 100
                        continue
                    if MHC_state == 1: # If already found a matching peptide then don't look at any more
                        continue
                    for MHC_PDB_1 in rows_MHC_PDB_1:
                        if MHC_state == 1: 
                            continue
                        if seq_length_PDB_2 == MHC_PDB_1[0]:
                            if sequence_PDB_2 == MHC_PDB_1[1]:
                                MHC_state = 1
                                check_MHC_list.append(1)
                
                elif molecule_type_PDB_2 == 'TCR':
                    if TCR_state == 1: # If already found a matching peptide then don't look at any more
                        continue
                    for TCR_PDB_1 in rows_TCR_PDB_1:
                        if TCR_state == 1: 
                            continue
                        if seq_length_PDB_2 == TCR_PDB_1[0]:
                            if sequence_PDB_2 == TCR_PDB_1[1]:
                                TCR_state = 1
                                check_TCR_list.append(1)
                    
            # If no matches found with PDB 2 then record a zero
            if pep_state == 0:
                check_pep_list.append(0)
            if MHC_state == 0:
                check_MHC_list.append(0)    
            if TCR_state == 0:
                check_TCR_list.append(0)  
    
        check_pep_rows.append(check_pep_list)
        check_MHC_rows.append(check_MHC_list)
        check_TCR_rows.append(check_TCR_list)
        
    # Create dataframes
    df_check_seq_peptide = pd.DataFrame(check_pep_rows,columns=PDB_code_list,index=PDB_code_list)
    df_check_seq_TCR = pd.DataFrame(check_TCR_rows,columns=PDB_code_list,index=PDB_code_list)
    df_check_seq_MHC = pd.DataFrame(check_MHC_rows,columns=PDB_code_list,index=PDB_code_list)    

    return df_check_seq_peptide, df_check_MHC_peptide, df_check_TCR_peptide