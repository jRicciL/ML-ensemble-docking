import requests
import json
from bs4 import BeautifulSoup
from prody import parsePDB
from typing import Tuple, List
import pandas as pd
from Bio import pairwise2, SeqIO
import time

def get_structure_sequence(pdb_id: str, 
                           chain: str = 'A') -> Tuple:
    '''
    Given a PDB identifier and a chain name of a protein 
    this function returns que aa sequence of the
    protein and the numbered positions
    '''
    # Load the structure and select the CA atoms, 
    # Ommit non standard residues and negative numbered residues
    ca_selection = f'protein and chain {chain} and ca and not nonstdaa and resid > 0'
    ref_struc = parsePDB(pdb_id).\
                select(ca_selection)
    # Get the protein sequence
    seq_cry = ref_struc.getSequence()
    # Get the initial and final positions
    positions = ref_struc.getResnums() 
    return seq_cry, positions


def get_seq_from_uniprot(uniprot_id: str, 
                         output_dir: str = './') -> str:
    ''' 
    Saves and returns the fasta sequence of a protein given
    its UNIPROT accession number
    '''
    URL= "https://www.uniprot.org/uniprot/"
    url_fasta = requests.get(URL + uniprot_id + ".fasta")
    
    file_name_fasta = output_dir + uniprot_id + '.fasta'
    open(file_name_fasta, 'wb').write(url_fasta.content)

    # Read the protein sequence
    fasta_prot = SeqIO.read(open(file_name_fasta),'fasta')
    seq_prot = str(fasta_prot.seq)
    return seq_prot


def pdb_ids_from_uniprot(uniprot_id: str, 
                         time_sleep: int = 2) -> pd.DataFrame:
    '''
    This function performs a web scrapping on the uniprot 
    entry page of a given protein accession number.
    The pdb ids related to the protein are listed and their 
    sequences are evaluated to determine if, based on
    their length and silimilarity, these structures could be used for futher analysis.
    '''
    r = requests.get('https://www.uniprot.org/uniprot/' + uniprot_id)
    soup = BeautifulSoup(r.content, "html.parser")
    pdb_tags = soup.find_all('a', {'class': 'pdb'})
    
    pdb_chains = []
    for tag in pdb_tags:
        # if it is a model, skip
        method = tag.find_next('td').text
        if method == 'model':
            continue
        pdb_id = tag.text
        resolution = tag.find_next('td').find_next('td').text
        chain = tag.find_next('td').\
                    find_next('td').\
                    find_next('td').text[0]
        seq_tag = tag.find_next('td').\
                      find_next('td').\
                      find_next('td').\
                      find_next('td')
        seq_text = seq_tag.text
        edges = [int(n) for n in seq_text.split('-')]
        range_len = len(range(edges[0], edges[1] + 1))
        pdb_chains.append((pdb_id.lower(), 
                            method, 
                            resolution, 
                            chain, 
                            edges[0], 
                            edges[1], 
                            range_len))
    
    time.sleep(time_sleep)
    # Return a dataframe
    result_df = pd.DataFrame(
        pdb_chains, 
         columns = ['pdb_id', 
                    'method', 
                    'resolution', 
                    'chain',
                    'start', 
                    'end', 
                    'seq_len'
                   ]
    )
    return result_df


def get_useful_pdbids(df_pdb_ids: pd.DataFrame, 
                      positions: List, 
                      thr_tol: Tuple = (3, 3)) -> pd.DataFrame:
    '''
    For each pdb id in the table, 
    select those with the correct range of aa given by the
    `positions` list. Additionaly, a tolerance value, given
    by the `thr_tol` parameter, allows to increase the range
    of positions.
    '''
    df_ids = df_pdb_ids.query(
        f'start <= {positions[0] + thr_tol[0]} and end >= {positions[-1] - thr_tol[1]}'
    )
    return df_ids


def get_bounded_ligands(pdb_id: str, 
                        entity_id: int = 1) -> List:
    '''
    Given a valid `pdb_id` and an entity number (`entity_d`)
    this function returns the names of the cocristalized ligands.
    '''   
    URL = f'https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/binding_sites/{pdb_id}/{entity_id}'
    r = requests.get(URL)
    ligands_names = []  
    try:
        content = json.loads(r.content)[pdb_id]
        ligands_data = content['data']
        for lig in ligands_data:
            lig_name = lig['accession']
            ligands_names.append(lig_name)
    except KeyError as err:
        None
    return ligands_names


def get_pdb_sequence(pdb_id: str, 
                        entity_id: int = 1) -> List:
    '''
    Given a valid `pdb_id` and an entity number (`entity_d`)
    this function returns the sequence of the protein entity
    '''  
    URL = f'https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/uniprot_mapping/{pdb_id}/{entity_id}'
    r = requests.get(URL)
    try:
        content = json.loads(r.content)[pdb_id]
        sequence = content['sequence']
    except KeyError as err:
        None
    return sequence


def get_identity(seq1: str, 
                 seq2: str, 
                 gap_char: str = '-') -> float:
    '''
    Computes the identity value between two protein sequences
    '''
    if not len(seq1) == len(seq2):
        return 0
    full_len = len(seq1)
    min_len = min(len(seq1.strip(gap_char)), len(seq2.strip(gap_char)))
    counter = 0
    for i, j in zip(seq1, seq2):
        if i == gap_char or j == gap_char:
            continue
        else:
            i == j
            counter += 1
    identity = counter / min_len
        
    return identity