import requests
from bs4 import BeautifulSoup
from prody import parsePDB
from typing import Tuple
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


