import re
from prody import parsePDB
from Bio import pairwise2, SeqIO
from typing import Tuple

def get_gaps_and_coverage(pdb_file: str, 
                          full_sequence: str,
                          chain: str = 'A') -> Tuple:
    '''
    Given a pdb file and a protein sequence, performs an alignment 
    between the pdb sequence and the `full_sequence`
    '''
    structure = parsePDB(pdb_file).getHierView()[chain]
    seq_query = structure.getSequence()
    alignment = pairwise2.align.globalxs(seq_query, full_sequence, 
                                            -10, -1, gap_char='-',
                                        one_alignment_only = True)[0]
    seq_alg = alignment[0]

    coverage = len(seq_query) / len(full_sequence) *100
    gaps = find_gaps(alignment[0])
    num_gaps = gaps["num_gaps"]
    gap_lengths = gaps["gap_lengths"]
    gap_list = gaps["gap_list"]

    return (seq_alg, coverage, gaps)



def find_gaps(seq: str, 
              r: int = 1):
    '''
    Given a sequence, this function finds the number and length of sequence gaps
    (defined by the `-` symbol)
    '''
    seq_len = len(seq)
    gaps = list(re.finditer('[-]+', seq))
    num_gaps = len(gaps); gap_lengths = []
    gap_list = []; gap_window = []
    # Get the start and final position of the gap
    for i , gap in enumerate(gaps, 1):
        start = gap.start() + 1 
        end = gap.end()
        gap_lengths.append(end - start + 1)
        gap_list.append([start, end])
        end_right = end if end + r >= seq_len else end + r
        start_right = start if start - r <= 1 else start - r
        gap_window.append([start_right, end_right])
    gaps_dict = {"num_gaps": num_gaps, 
                 "gap_lengths":gap_lengths,
                 "gap_list": gap_list, 
                 "gap_window": gap_window}
    return gaps_dict


def is_a_gap(x, seq_intervals): 
    gap_res = []
    for i in x: 
            for j in range(i[0], i[1] + 1 ): gap_res.append(j)
    return set(gap_res).isdisjoint(seq_intervals)



