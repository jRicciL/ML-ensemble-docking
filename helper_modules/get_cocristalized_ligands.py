'''Module for get cocristalized ligands inside protein pocket/active sites'''
from glob import glob
from prody import *
import nglview as nv
import numpy as np

class PocketResidues:

    def __init__(self, pdb_id, ligand_resname):
        self.pdb_id = pdb_id
        self.ligand_resname = ligand_resname
        # Load the molecule
        self.system = parsePDB(pdb_id)
        self.protein = self.system.select('protein')
        self.ligand  = self.system.select(f'resname {ligand_resname}')
        assert self.ligand is not None
        # List of residues
        self.pocket_res_list = None


    def get_pocket_residues_as_list(self, cutoff=5):
        '''
        Returns a list of residues surrounding
            Parameters:
                cutoff (int): Distance cutoff to select the protein atoms surrounding the ligand
            Returns:
                pocket_residues (list): A list of residue numbers of the protein
        '''

        pocket_sel = self.system.select(
                        f'within 7 of resname {self.ligand_resname} and protein')
        pocket_residues = [str(i) 
                        for i in np.unique(
                            pocket_sel.getResnums())
                        ]
        pocket_list = ' '.join(pocket_residues)
        self.pocket_res_list = pocket_list
        return pocket_list

    def get_pocket_residues_as_str(self, cutoff=5):    
        '''Converts a list of resid selections to a formated string of atom selections'''
        pocket_list = self.get_pocket_residues_as_list(cutoff = cutoff)
        str_of_prody_selections = ''.join(str(pocket_list)).replace('[', '').replace(']', '').replace(',', '')
        return str_of_prody_selections

    def visualize_pocket(self, cutoff=5):
        str_residues = self.get_pocket_residues_as_str(cutoff = cutoff)
        pocket_atoms = self.protein.select('resid ' + str_residues).getIndices()

        # Load ligand as nglview object
        ligand_nv = nv.ProdyStructure(self.ligand)

        view = nv.show_prody(self.protein)
        view.clear_representations()
        view.add_representation('cartoon', selection='protein', color='white')
        view.add_licorice(selection=pocket_atoms, color='red')
        view.add_cartoon(selection=pocket_atoms, color='red')
        view.add_structure(ligand_nv)
        return view



def get_pocket_ligand(
        pdb_id: str,
        pocket_residues: str,
        raw_lig_dir: str,
        prot_chain_dir: str,
        pk_ligs_dir: str,
        cutoff: int = 3,
        min_weight: float = 97.0, # mw sulfate + 1
        write_files: bool = True):
    
    # 1) Load the molecules using prody
    try:
        lig_file = glob(f'{raw_lig_dir}/{pdb_id}*.pdb')[0]
        lig = parsePDB(lig_file)
        
    except IndexError as e:
        print("The protein", pdb_id, "have NO LIGAND.")
        return(None, None)
    
    # 2) Se carga la proteína y se seleccionan los residuos que definen el pocket
    try:
        prot_file = glob(f'{prot_chain_dir}/{pdb_id}*.pdb')[0]
        protein = parsePDB( prot_file )
    except IndexError as e:
        print("Protein file not found:", pdb_id)
        return(None, None)
    
    protein_pocket = protein.select("resnum " + pocket_residues)
    
    # Se seleccionan los ligandos (RESIDUOS no protéicos) que estén a no más de cutoff A
    # de cualquier átomo de los residuos del pocket de la proteína
    lig_sel = lig.select('within ' + str(cutoff) + ' of prot', prot = protein_pocket)
    if lig_sel is None:
        print("The protein", pdb_id, "has no ligand inside the pocket.")
        return(None, None)
    
    # Se obtine la lista de moléculas que cumplen el criterio anterior
    inhibidor_list = np.unique( lig_sel.getResnames() )
    
    # Calcula el centro geométrico del pocket de la proteína
    prot_pocket_center = calcCenter(protein_pocket)
    
    # Puede que haya más de una molécula en el pocket 
    # (debido a presencia de iones, cosolvente o residuos modificados)
    # Se itera entre cada molecula del ligando, se calcula su masa y se mantiene el ligando con mayor masa 
    distance = 10000
    nearest_chain = ''
    nearest_resnum = ''
    nearest_resname = ''
    current_mass = 0
    true_lig = None

    for resname in inhibidor_list:
            mol = lig.select("resname " + resname)
            
            # Itera entre todos las moléculas definidas por un resname, un resnum y una cadena
            mol_chains = np.unique(mol.getChids())
            for chain in mol_chains:
                mol_resnums = np.unique(mol.select('chain ' + chain).getResnums())
                for resnum in mol_resnums:
                    new_mol = mol.select('chain ' + chain + \
                                  ' and resnum ' + str(resnum))
                    new_mass = new_mol.getMasses().sum()
                    
                    #print(new_mass, resname, chain)
                    if new_mass > current_mass:
                    
                    
                        dist_new_mol = calcDistance(prot_pocket_center, calcCenter(new_mol))
#                     #print(chain, resnum, resname, dist_new_mol)
#                     if dist_new_mol < distance:
                        nearest_chain = chain
                        nearest_resnum = str(resnum)
                        nearest_resname = resname
                        distance = dist_new_mol
                        current_mass = new_mass
                        
    true_lig = lig.select('chain ' + nearest_chain + \
                          ' and resnum ' + nearest_resnum + \
                          ' and resname ' + nearest_resname)
    lig_mass = true_lig.getMasses().sum()

            
    # Se guarda el ligando sólo si su masa es igual o mayor a la de un sulfato
    # si el ligando aparece 2 veces cerca del pocket, se elige la molécula más cercana al centro
    if true_lig != "" and lig_mass > min_weight: # Si el ligando existe y su masa fue mayor 96
        # Se extrae el nombre del ligando
        name_lig = np.unique( true_lig.getResnames() )[0]

        # Se combinan proteína y ligando en un solo objeto AtomGroup, de tal manera que
        # sea posible superponer las estructuras de la proteína de acuerdo a los residuos de su pocket
        # y esto afcete la psoición del ligando
        complejo_PL = protein + true_lig.toAtomGroup()
        if write_files:
            # Se guardan ligando y proteína alineados según la estructura de referencia
            print(F'Protein {pdb_id}:\n   Molecules found: {str(inhibidor_list)} -> ligand {name_lig} saved.')
            writePDB( pk_ligs_dir + "/" + pdb_id + "_" + name_lig + "_LIG.pdb", 
                      complejo_PL.select("resname " + name_lig) )
        return(name_lig, lig_mass)
    else: 
        print(F"The model {pdb_id} HAS NO LIGAND inside the pocket.")
        return(None, None)
