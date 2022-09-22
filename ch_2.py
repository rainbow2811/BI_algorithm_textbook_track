# %% ----------------------------------------
import numpy as np
from src import pattern_base, hamming_distance_base
import importlib
import typing as tp
import itertools as it
Kmer = tp.AnyStr
# %% ----------------------------------------
importlib.reload(pattern_base)
# %%    Q) 2A
k = 6
d = 1
# %% ----------------------------------------
dnas = """
GTAGTACCCACGATACAGTAGTAGCATAAAATTGCAATGTAC
GTAGTTACATTTTCTCCGGTTGTTAACTGTCACGGGATATAT
GGCAGTAACGTGACCTGATTAGTCCAACGAGTAGTTCGTTCG
GTGTCTTAGCGCGGCATGCACGCCGGGGATGTAGTTTTGTGT
GATCAGGGAAGATGATCCTTGAATACTAGCCATCGTGTAGTG
TCGGACTCATTTTCAGGGCATCGTCACAGGGTAGTGCCCTAA
GACCTGCGTAAATCACGTGGGCTGACCGTGATTTTCGTAGTG
ACTCCTACCGATTGTGATACTCTCGTAGTTGTGAGATCATTT
TACCCGTTCCGGCAGGCAGATTCGTGAGTTGTAGTAATCGGT
TCTTGTTAGGGACACAAGTTGTGATAATAGAATATGGTAGTA
"""
dnas
# %% ----------------------------------------
# _dna_list = dnas.strip().split('\n')
# _dna_list
# # %% ----------------------------------------
# dna_list = list()
# for dna in _dna_list:
#     text_list = dna.split(' ')
#     dna = text_list[-1]
#     dna_list.append(dna)
# # # %% ----------------------------------------
# dna_list = [dna.split(' ')[-1] for dna in _dna_list]
# dna_list
# # %% ----------------------------------------
# def preprocess_dna(raw_dna):
#     """
#     ex)

#     raw_dna = '1  ATTTGGC'

#     raw_dna.split(' ') --> ['1', '', 'ATTTGGC']

#     output_dna = 'ATTTGGC'

#     """
#     output_dna = raw_dna.split(' ')[-1] 

#     return output_dna
# %% ----------------------------------------
# list(map(preprocess_dna, _dna_list))
# %% ----------------------------------------
def get_dna_list_from_raw_data(dna_data):
    _dna_list = dna_data.strip().split('\n')
    dna_list = [dna.split(' ')[-1] for dna in _dna_list]
    return dna_list
# %% ----------------------------------------
dna_list = get_dna_list_from_raw_data(dnas)
dna_list
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
# " ".join(motif_enumeration_psp(dnas, k, d))
# %%    Q) 모티프 찾기 문제 (p110)
# 모티프 점수 최소화하는 kmer 집합 찾기 (Brute Forece Motif Search)
def get_all_possible_kmer_matrix(dna_list, k):
    matrix_for_kmer_list = list()
    for dna in dna_list:
        kmer_list = list(pattern_base.pattern_matrix(dna,k))
        matrix_for_kmer_list.append(kmer_list)
    return matrix_for_kmer_list
# %%
def get_all_possible_motif_matrix(matrix_for_kmer_list):
    combination = list(it.product(*matrix_for_kmer_list))
    return combination
# %%
def get_motif_nt_matrix(one_motif_set):
    kmer_motif_matrix = list(map(lambda kmer:list(kmer), one_motif_set ))
    return np.array(kmer_motif_matrix)
# %%
def get_score_from_nt_count_arr(column_nt_matrix):
    _, counts = np.unique(column_nt_matrix, return_counts=True)
    score = np.sum(counts)-np.max(counts)
    return score
# %%
def get_score(nt_motif_matrix):
    column_nt_matrix = list(zip(*nt_motif_matrix))
    return list(map(get_score_from_nt_count_arr, column_nt_matrix))
# %%
def get_score_from_one_motif_set(one_motif_set):
    nt_motif_matrix = get_motif_nt_matrix(one_motif_set)
    return np.sum(get_score(nt_motif_matrix))
# %%
def find_min_score_motif_set(dna_list, k):
    kmer_matrix = get_all_possible_kmer_matrix(dna_list, k)
    comb = get_all_possible_motif_matrix(kmer_matrix)
    motif_score_result = list(map(get_score_from_one_motif_set, comb))
    min_score = np.min(motif_score_result)
    min_score_idx = np.where(motif_score_result==min_score)[0][0]
    return comb[min_score_idx], f"min_score = {min_score}"
# %% ----------------------------------------
%timeit find_min_score_motif_set(dna_list, k)
# %%  Q) 중앙 문자열 해결 (p141, 2H)
def distance_btw_pattern_and_string(pattern, dna_list):
    k = len(pattern)
    distance = 0
    for string in dna_list:
        hamming_distance = len(string) * 2
        pattern_set_from_one_string = pattern_base.pattern_matrix(string, k)
        for _pattern in pattern_set_from_one_string:
            h_distance = hamming_distance_base.hamming_distance(pattern, _pattern)
            if hamming_distance > h_distance:
                hamming_distance = h_distance
        distance += hamming_distance
    return distance
# %%
# distance_btw_pattern_and_string('ATT', dna_list)
# %% ----------------------------------------
# class DNA:
#     # def __init__(self, str1=None):
#     #     self.str1 = str1
#     pass
# # %% ----------------------------------------
# dna = DNA()
# # %% ----------------------------------------
# dna.str1 = 'TTACCTTAAC'
# dna.str2 = 'GATATCTGTC'
# dna.str3 = 'ACGGCGTTCG'
# dna.str4 = 'CCCTAAAGAG'
# dna.str5 = 'CGTCAGAGGT'
# # dna_list = 
# # %% ----------------------------------------
# dna_list = [ getattr(dna, f"str{n}") for n in range(1,6) ]
# # %%
# dna_list
# # %% ----------------------------------------
# k = 3
# # %%
# distance_btw_pattern_and_string('AAA',dna_list)
# %%
def median_string(dna_list: tp.List[Kmer], k:int) -> Kmer:
    distance = len(dna_list[0]) * 2
    kmer_base_list = pattern_base.kmer_base(k)
    for kmer_pattern in kmer_base_list:
        calculated_distance = distance_btw_pattern_and_string(kmer_pattern, dna_list)
        if distance > calculated_distance:
            distance = calculated_distance
            median_pattern = kmer_pattern
    return median_pattern
# %%
%timeit median_string(dna_list, k)
# %% 
median_string(dna_list,k)
# %% Q) most probable kmer 문제(p118)
text = 'AAAGGCGTGTCTCCGGAAATATCAGAGGCTAATCGTGTTTACAATGACAAGAGTTACTCTTGTGCATCACATTGTCAATATTTGGATGTCACCTACAGCTATCTCTAGACCACGGGCGCATCCGGGCCCCATGGGGTTAAACCCACAATCACACTTACAGGTAATTTGGATGTCGTCTGACCATGTCTACTACCACGGGT'
k = 8
profile_str = """
0.16 0.16 0.16 0.28 0.2 0.12 0.2 0.32
0.4 0.4 0.28 0.28 0.36 0.28 0.36 0.16
0.28 0.24 0.36 0.28 0.2 0.36 0.28 0.08
0.16 0.2 0.2 0.16 0.24 0.24 0.16 0.44
"""
# %%
# profile_str
# # %%
# _profile_list = profile_str.strip().split('\n')
# _profile_list
# # %%
# profile_list = list()
# for row_str in _profile_list:
#     x = row_str.split(' ')
#     profile_list.append(x)
# profile_list
# # %%
# profile_matrix = np.array(profile_list)
# profile_matrix.shape
# %% ----------------------------------------
def get_profile_matrix_from_profile_text(raw_profile):
    _profile_list = profile_str.strip().split('\n')
    profile_list = list()
    for row_str in _profile_list:
        x = row_str.split(' ')
        profile_list.append(x)
    return np.array(profile_list)
# %% ----------------------------------------
def get_prob_kmer(kmer, profile_matrix):
    base = dict(A=0,C=1,G=2,T=3)
    prob_nt_list = list()
    for i, nt in enumerate(kmer):
        actg_idx = base[nt]
        prob_nt = float(profile_matrix[ actg_idx, i ])
        prob_nt_list.append(prob_nt)
    prob_kmer = np.prod(prob_nt_list)
    return prob_kmer
# %%
def get_most_probable_kmer(text, k, text_profile):
    profile_matrix = get_profile_matrix_from_profile_text(text_profile)

    kmer_prob = 0
    kmer_num = len(text)-k+1
    probable_kmer = list()
    for i in range(kmer_num):
        kmer = text[i:i+k]
        _prob_kmer = get_prob_kmer(kmer, profile_matrix)
        if kmer_prob < _prob_kmer:
            kmer_prob = _prob_kmer
            probable_kmer = kmer
    return probable_kmer
# %%
get_most_probable_kmer(text, k, profile_str)
# %%
