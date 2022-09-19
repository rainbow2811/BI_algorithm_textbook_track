# %% ----------------------------------------
import numpy as np
from src import pattern_base, hamming_distance_base
import importlib
import typing as tp
Kmer = tp.AnyStr
# %% ----------------------------------------
importlib.reload(pattern_base)
# %%    Q) 2A
k = 5
d = 2
# %% ----------------------------------------
dnas = """
CTTAACAACGTAGGGAGCCGAGCGT
TGCAGATTAACCAACCAATGACCCA
CACCGTGTGCTATGTCTGTTAGCCT
GCACTAGCCAGAAGAACGGCATAGG
ATAAGAAGAGAACCAGGTCCCGAGG
TTGTGCTCATTAGTCAGCCGTGGAG
AGTTTAGGGGAACCGGGGCCAGCCT
ATTGAGGATGTTACAGCGTCACCCG
CCCCCTCTCTTACACAACCTGACTC
GCACAATCCGCATCTAAAATAGCGC
"""
dnas
# %% ----------------------------------------
_dna_list = dnas.strip().split('\n')
_dna_list
# %% ----------------------------------------
dna_list = list()
for dna in _dna_list:
    text_list = dna.split(' ')
    dna = text_list[-1]
    dna_list.append(dna)
# %% ----------------------------------------
dna_list = [dna.split(' ')[-1] for dna in _dna_list]
dna_list
# %% ----------------------------------------
def preprocess_dna(raw_dna):
    """
    ex)

    raw_dna = '1  ATTTGGC'

    raw_dna.split(' ') --> ['1', '', 'ATTTGGC']

    output_dna = 'ATTTGGC'

    """
    output_dna = raw_dna.split(' ')[-1] 

    return output_dna
# %% ----------------------------------------
# list(map(preprocess_dna, _dna_list))
# %% ----------------------------------------
def get_dna_list_from_raw_data(dna_data):
    _dna_list = dna_data.strip().split('\n')
    dna_list = [dna.split(' ')[-1] for dna in _dna_list]
    return dna_list
# %% ----------------------------------------
def get_kmer_neighbors_set_from_one_dna_string(dna_string, k, d):
    kmer_pattern_matrix = pattern_base.pattern_matrix(dna_string, k)
    total_kmer_neighbor_set = set()
    for kmer in kmer_pattern_matrix:
        kmer_neighbors_list = hamming_distance_base.neighbors_distance(kmer,d)
        total_kmer_neighbor_set.update(kmer_neighbors_list)
    return total_kmer_neighbor_set
# %% ----------------------------------------
def check_all_True_or_not(data: tp.List[tp.Set[Kmer]], kmer_neighbor:Kmer) -> bool:
    check = list()
    for set_ in data:
        # print( f"==== set{n+1} ====")
        check.append(int(kmer_neighbor in set_))
        # print(f"set{n+1}: {int(kmer_neighbor in s)}")
    # print(check)
    result_check_all_True =  np.prod(check)
    return result_check_all_True
# %% ----------------------------------------
def motif_enumeration(dna_data, k, d):
    dna_list = get_dna_list_from_raw_data(dna_data)
    num_of_dna_list = len(dna_list)

    list_of_kmer_neighbor_set_from_each_string = list()
    for i in range(1, num_of_dna_list):  # 첫번째 dna string 제외한 나머지 string의 neighbors_list 구하기
        dna_string = dna_list[i]
        kmer_neighbor_set= get_kmer_neighbors_set_from_one_dna_string(dna_string,k,d)
        list_of_kmer_neighbor_set_from_each_string.append(kmer_neighbor_set)

    kmer_motif_set = set()
    first_dna = dna_list[0]
    kmer_neighbor_set_from_first_dna_string = get_kmer_neighbors_set_from_one_dna_string(first_dna,k,d)
    for neighbor in kmer_neighbor_set_from_first_dna_string:
        result = check_all_True_or_not(list_of_kmer_neighbor_set_from_each_string, neighbor)
        if result:
            kmer_motif_set.add(neighbor)
    return kmer_motif_set
# %% ----------------------------------------
motif_enumeration(dnas,k,d)
# %% ----------------------------------------
def get_kmer_neighbor_matrix(dna_list, k, d):
    kmer_neighbor_matrix = list()
    for dna in dna_list:
        kmer_neighbor_set = get_kmer_neighbors_set_from_one_dna_string(dna, k, d)
        kmer_neighbor_matrix.append(kmer_neighbor_set)
    return kmer_neighbor_matrix
# %% ----------------------------------------
def neighbor_in_all_dna(neighbor, kmer_neighbor_matrix):
    cond_list = list()
    for i in range(1,len(dna_list)):
        cond = neighbor in kmer_neighbor_matrix[i]
        cond_list.append(cond)
    return all(cond_list)
# %% ----------------------------------------
def get_motif_set_from_neighbor_matrix(kmer_neighbor_matrix):
    motif_set = set()
    for neighbor in kmer_neighbor_matrix[0]:
        if neighbor_in_all_dna(neighbor, kmer_neighbor_matrix):
            motif_set.add(neighbor)
    return motif_set
# %% ----------------------------------------
def motif_enumeration_psp(dna_data, k, d):
    dna_list = get_dna_list_from_raw_data(dna_data)
    kmer_neighbor_matrix = get_kmer_neighbor_matrix(dna_list, k, d)
    motif_set = get_motif_set_from_neighbor_matrix(kmer_neighbor_matrix)
    return motif_set
# %% ----------------------------------------
" ".join(motif_enumeration_psp(dnas, k, d))
# %%
