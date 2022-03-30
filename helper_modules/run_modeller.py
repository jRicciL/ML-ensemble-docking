import shutil
import os
from prody import parsePDB, writePDB
from Bio import pairwise2, SeqIO
from modeller import *
from modeller.automodel import *
from helper_modules.find_gaps import *
from glob import glob
from typing import List
import numpy as np


def missing_elements_from_resnums(resnums: List[int], 
                                  max_pos: int) -> List[int]:
    start, end = resnums[0], max(resnums[-1], max_pos)
    missed_elements = sorted(
        set(range(start, end + 1)).difference(resnums))
    return missed_elements

def missing_elements_from_aa_sequence(seq: str):
    missed_elements = [i for i, c 
                        in enumerate(seq, 1) if c == '-' ]
    return missed_elements

def run_modeller(pdb_file: str,
                 target_sequence: str, 
                 output_dir: str = './',
                 keep_original_resnum: bool = True,
                 num_res_window: int = 2, 
                 overwrite: bool = True,
                 max_var_iterations: int = 500, 
                 repeat_optimization: int = 2,
                 chid: str = 'A',
                 verbose: bool = True,
                 start_position = None,
                 end_position = None
                ) -> None:
    '''
    Given a pdb structure with missing regions and the full protein sequence, 
    this function models missing loops inside the structure using Modeller
    '''
    
    base_name = pdb_file.split('/')[-1].split('.')[0]
    output_model_file = f'{output_dir}/{base_name}_mod.pdb'

    # Skip modelling if the model file already exists
    if not overwrite:
        if os.path.isfile(output_model_file):
            return print("Model already exists:", base_name)
    try:
        # Read the pdb file, omit non standard residues and residues with negative numbering
        stc_prot = parsePDB(pdb_file)
        if not start_position is None and not end_position is None:
            try:
                stc_prot = stc_prot.select(
               f'resnum {start_position} to {end_position}').toAtomGroup()
            except:
                None
        stc_prot = stc_prot.select('not nonstdaa and resid > 0') 
        # Get the sequence of the given structure
        seq_cry = stc_prot.select('ca').getSequence()
    except FileNotFoundError as e:
        print(e, "Error at opening:", base_name)
        return
    
    # Asks if the sequences have the same length
    # assert (len(seq_cry) == len(target_sequence))
    # assert (seq_cry == target_sequence)
    
    # Force the structure chain to be name as `A`
    stc_prot.setChids("A")
    
    crys_base_name = pdb_file
    model_base_name = base_name + '.modeller'

    # Performs the sequence alignment
    # pairwise returns multiple possible alingments with similar scores
    alignments = pairwise2.align.globalms(
        sequenceA = seq_cry, 
        sequenceB = target_sequence,  
        match    = 3, 
        mismatch = -1, 
        open     = -5, 
        extend   = 0, 
        gap_char = '-')
    # Find the best alingment
    stc_prot_resnums = np.unique(stc_prot.getResnums())
    stc_prot_missed_residues = missing_elements_from_resnums(
                                    stc_prot_resnums, 
                                    len(target_sequence)
                                )
    
    algn1_struc = alignments[0][0]
    for algn in alignments:
        algn_s = algn[0]
        algn_missed_residues = missing_elements_from_aa_sequence(algn_s)
        if np.all(algn_missed_residues == stc_prot_missed_residues):
            algn1_struc = algn_s
            break

    # Alignment sequences
    algn2_seq   = target_sequence

    ''' NEEDED: There should be 10 fields separated by colons ":".
    Please check the file to make sure your sequences end with the '*' character.
    Nomenclaturas de los campos del header: 
    https://salilab.org/modeller/8v2/manual/node176.html'''

    #  Create the alignment headers
    struc_header = "structureX:" + crys_base_name + ":.:" + chid + ":.:" + chid + ":.:.:.:"
    seq_header = "sequence:" + model_base_name + ":.:.:.:.:.:.:.:"
    
    # Create the alignment file
    alg_filename = base_name + ".alg"
    with open(alg_filename, "w") as handle:
        handle.write("\n>P1;%s\n%s\n%s*\n>P1;%s\n%s\n%s*\n" % (
                crys_base_name, 
                struc_header, algn1_struc, 
                model_base_name, seq_header, algn2_seq
            )
        )

    # Identify the gaps inside the structure sequence
    gaps = find_gaps(algn1_struc, r = num_res_window)
    num_gaps = gaps["num_gaps"]
    gap_i    = gaps["gap_window"]
    
    if verbose:
        print("Modelling protein " + base_name)
        print('Target sequence\n', target_sequence)
        print('Input sequence from .pdb file\n', algn1_struc)
        print('*****')
        print(gaps)
    
    # The following string indicates the position of the first gap
    if num_gaps > 0:
        s = f"self.residue_range('{str(gap_i[0][0])}:{chid}', '{str(gap_i[0][1])}:{chid}')"

        # If there are more gaps, these will be added to the initial string
        for i in range(1, num_gaps):
            s = s + f", self.residue_range('{str(gap_i[i][0])}:{chid}', '{str(gap_i[i][1])}:{chid}')" 
    else:
        s = s = f"self.residue_range('1:{chid}', '1:{chid}')"
        
    #*******************#
    ''' RUN MODELLER '''
    #*******************#

    # Create a Modeller environment
    env = Environ()
    env.io.atom_files_directory = ['.', '.']

    # Modify the `MyModel` class. Specifying which residues will be modeled.
     
    MyModel_code = """
class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(""" + s + """)
""" # the identation is important here
    exec(MyModel_code, globals()) 
    # Reads the alignment file
    a = MyModel(env, alnfile = alg_filename, 
                      knowns = pdb_file, 
                      sequence = model_base_name) 
    # Create just one model
    a.starting_model = 1
    a.ending_model   = 1
    # Define the modelling parameters:
    # more info about these parameters https://salilab.org/modeller/9.21/manual/node19.html
    a.library_schedule = autosched.slow 
    a.max_var_iterations = max_var_iterations
    a.md_level = refine.slow 
    a.repeat_optimization = repeat_optimization
    a.make()

    ###########################
    # These steps will clean the unnecessary files and will rename the output model
    model_file = glob('./' + model_base_name + '*.pdb')
    assert len(model_file) == 1 # There should be only one file
    model_file = model_file[0]
    # Rename the model file
    os.rename(model_file, output_model_file) 

    # Keep the original numbering:
    if keep_original_resnum:
        initial_pos = stc_prot.getResnums()[0]
        print(initial_pos)
        # Open the model
        ref_model = parsePDB(output_model_file)
        # Renumber
        ref_model.setResnums(ref_model.getResnums() + initial_pos - 1)
        ref_model.setChids('A')
        writePDB(output_model_file, ref_model)

    # Delete nonuseful files
    os.rename(f'{base_name}.alg', f'{output_dir}/{base_name}.modeller.alg')
    for f in glob(model_base_name + "*"):
        os.remove(f)