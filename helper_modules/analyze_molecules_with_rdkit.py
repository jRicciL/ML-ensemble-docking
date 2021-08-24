from glob import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import os
from typing import Callable, List, Dict

def get_files_list(path_to_sdfs: str, 
                   actives_name: str = 'ligand',
                   inactives_name: str = 'decoy',
                   sufix: str = '',
                   sort_func: Callable = (
                       lambda x: int(x.split('/')[-1].split('.')[0].split('_')[1])
                       ) 
                  ) -> List:
    '''Returns a list of path sdf files in a given directory,
    a sort function could be parsed depending on the pattern file_names'''
    # Active molecules
    file_list_ACTIVES = glob(path_to_sdfs + F'/*{actives_name}*{sufix}*.sdf')
    file_list_ACTIVES.sort(key = sort_func) 
    # Inactive molecules
    file_list_INACTIVES = glob(path_to_sdfs + F'/*{inactives_name}*{sufix}*.sdf')
    file_list_INACTIVES.sort(key = sort_func)
    # Join both list
    file_list = file_list_ACTIVES + file_list_INACTIVES
    return file_list

def load_molecules_from_dir(list_of_sdf_files: List) -> Dict:
    '''Function to load molecules from sdf files using rdkit'''
    # Load the molecules in a dictionary
    mols_dict = {}
    sanitized = True
    for sdf_file in list_of_sdf_files:
        # Get the molecule name
        mol_name = sdf_file.split('/')[-1].split('.')[0]
        # Try to load the molecule with sanitize = True
        mol_rd = Chem.SDMolSupplier(sdf_file, sanitize = True)[0]
        if mol_rd is None:
        	try:
	            mol_rd = Chem.SDMolSupplier(sdf_file, sanitize = False)[0]
	            mol_rd.UpdatePropertyCache(strict = False)
	            sanitized = False
	        except AttributeError:
	        	print(f'Error with', mol_name)
	        	mol_rd, sanitized = None, None
        mols_dict[mol_name] = [mol_rd, sanitized]
    return mols_dict

def get_mol_dataframe(mol_dictionary: Dict) -> pd.DataFrame:
    '''Turns a dictionary of molecules into a dataframe'''
    # Convert to a dataframe
    df = pd.DataFrame(mol_dictionary).T
    df.columns = ['mol_rdk', 'sanitized']
    # Activity inactivity column
    act_inact = ['active' if i[:6] == 'ligand' else 'inactive' for i in df.index]
    df['Activity'] = act_inact
    # Naming the columns
    df = df[['Activity', 'mol_rdk', 'sanitized']]
    return df


def load_cocrys_molecules_from_dir(list_sdf_files: List) -> pd.DataFrame:
    '''Function to load molecules from sdf files using rdkit.
    This function load molecules obtained from pdb files (cocristalized with proteins).
    The suffix of each file should be "_from_pdb.sdf" and the directory also shoud
    have a corresponding file "_from_mol2.sdf"'''
    # Get a dataframe with pdbis: lig_names
    pdbId_lig_dic = get_protId_ligName_dic(list_sdf_files)
    
    # Load the molecules in a dictionary
    mols_dict = {}
    ligs_validation_dic = {}
    
    sanitized = True
    for file in list_sdf_files:
        #assert sdf_file[:-12] == 'from_pdb.sdf'
        # Get the pdb id
        pdb_id = file.split("/")[-1].split("_")[0]
        # FIRST TRY: Load the sdf from 'from_pdb.sdf' file
        mol_rd = Chem.SDMolSupplier(file, sanitize = True)[0]
        mols_dict[pdb_id] = mol_rd
        ligs_validation_dic[pdb_id] = 'v1'
        # SECOND TRY: Load the molecule from 'from_mol2.sdf' file
        if mols_dict[pdb_id] is None:
            file = file.replace('_pdb.sdf', '_mol2.sdf')
            if not os.path.isfile(file):
                print("File '_from_mol2.sdf' doesn't exist.")
                return
            mol_rd = Chem.SDMolSupplier(file, sanitize = True)[0]
            mols_dict[pdb_id] = mol_rd
            ligs_validation_dic[pdb_id] = 'v2'
            # FINAL TRY: Molecule is read from '_from_mol2.sdf' without sanitize it
            if mols_dict[pdb_id] is None:
                try:
                    mol_rd = Chem.SDMolSupplier(file, sanitize = False)[0]
                    # Force strict to false to let use some descriptors without sanitization
                    mol_rd.UpdatePropertyCache(strict = False)
                    mols_dict[pdb_id] = mol_rd
                    ligs_validation_dic[pdb_id] = 'v3'
                except AttributeError:
                    print(f'Error with', pdb_id)
                    mol_rd, sanitized = None, None
                    mols_dict[pdb_id] = mol_rd
                    ligs_validation_dic[pdb_id] = 'error'
    # Create the final dataframe
    df = pd.DataFrame.from_dict(pdbId_lig_dic, orient='index', columns=["Lig"])
    df["mol_rdk"] = mols_dict.values()
    df["validation"] = ligs_validation_dic.values()
    return df

def get_protId_ligName_dic(list_sdf_files: List) -> Dict:
    '''Returns a dictionary of pdbids as keys and LigNames (LIG) as values
    given a list of sdf files named as "pdid_LIG_xxxxx.sdf"'''
    pdbId_lig_dic = {file.split("/")[-1].split("_")[0]: file.split("/")[-1].split("_")[1] 
               for file in list_sdf_files}
    return pdbId_lig_dic
