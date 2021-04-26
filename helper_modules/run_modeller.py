import shutil
import os
from prody import parsePDB, writePDB
from Bio import pairwise2, SeqIO
from modeller import *
from modeller.automodel import *
from find_gaps import *
from glob import glob


def run_modeller(pdb_file: str, 
                 seq_prot: str, 
                 output_dir: str = './', 
                 keep_original_resnum: bool = True,
                 num_res_window: int = 2, 
                 force_modeling: bool = False,
                 max_var_iterations: int = 500, 
                 repeat_optimization: int = 2,
                 chid: str = 'A',
                 verbose: bool = True
                ) -> None:
    '''
    Given a pdb structure with missing regions and the full protein sequence, 
    this function models missing loops inside the structure using Modeller
    '''
    
    file_name = pdb_file.split('/')[-1].split('.')[0]
    output_file_name = output_dir + file_name + '_mod.pdb'

    # Skip modelling if the model file already exists
    if not force_modeling:
        if os.path.isfile(output_file_name):
            return print("Model already exists:", file_name)
    #
    try:
        # Read the pdb file, omit non standard residues and residues with negative numbering
        stc_prot = parsePDB(pdb_file)
        print(stc_prot.getResnames())
        stc_prot = stc_prot.select('not nonstdaa and resid > 0') 
        print(pdb_file)
        # Get the sequence of the given structure
        seq_cry = stc_prot.select('ca').getSequence()
        print(seq_cry)
    except FileNotFoundError as e:
        print(e, "Error at opening:", file_name)
        return
    
    # Asks if the sequences have the same length
    same_seq = len(seq_cry) == len(seq_prot) and seq_cry == seq_prot
    
    # Force the structure chain to be name as `A`
    stc_prot.setChids("A")
    
    crys_file_name = pdb_file
    model_file_name = file_name + '_mod'

    # Performs the sequence alignment
    alignment = pairwise2.align.globalms(seq_cry, seq_prot,  3, -1, -10, -.1, gap_char='-')[0]
    # Alignment sequences
    algn1_struc = alignment[0]
    algn2_seq = seq_prot

    ''' NEEDED: There should be 10 fields separated by colons ":".
    Please check the file to make sure your sequences end with the '*' character.
    Nomenclaturas de los campos del header: 
    https://salilab.org/modeller/8v2/manual/node176.html'''

    #  Create the alignment headers
    struc_header = "structureX:" + crys_file_name + ":.:" + chid + ":.:" + chid + ":.:.:.:"
    seq_header = "sequence:" + model_file_name + ":.:.:.:.:.:.:.:"
    
    # Create the alignment file
    alg_filename = file_name + ".alg"
    with open(alg_filename, "w") as handle:
        handle.write("\n>P1;%s\n%s\n%s*\n>P1;%s\n%s\n%s*\n" % (crys_file_name, struc_header, algn1_struc, 
                                                               model_file_name, seq_header, algn2_seq))

    # Identify the gaps inside the structure sequence
    gaps = find_gaps(algn1_struc, r = num_res_window)
    num_gaps = gaps["num_gaps"]
    gap_i = gaps["gap_window"]
    
    if verbose:
        print("Modelling protein " + file_name)
        print(seq_cry, seq_prot)
        print('*****')
        print(alignment)
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
                      sequence = model_file_name) 
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
    model_file = glob('./' + file_name + '*.pdb')[0]
    output_model_file = output_dir + f'/{file_name}_mod.pdb'
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
    for f in glob(file_name + "*"):
        os.remove(f)