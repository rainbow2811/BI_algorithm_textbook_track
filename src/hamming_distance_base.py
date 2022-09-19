# %% ----------------------------------------
import numpy as np
# %% ----------------------------------------
def hamming_distance(pattern1, pattern2):
    match_list = list(map(lambda x, y: 0 if x == y else 1, pattern1, pattern2))
    return np.sum(match_list)
# %% ----------------------------------------
def approximate_patt_count(text, pattern, d):
    count = 0
    for i in range(len(text)-len(pattern)):
        compared_pattern = text[i:i+len(pattern)]
        if hamming_distance(pattern, compared_pattern) <= d:
            count += 1
    return count
# %% ----------------------------------------
def hamming_distance_psp(pattern1:str, pattern2:str):
    a = np.frombuffer(pattern1.encode(), dtype=np.uint8)
    b = np.frombuffer(pattern2.encode(), dtype=np.uint8)
    return np.nonzero(a - b)[0].shape[0]
# %% ----------------------------------------
def neighbors(pattern):
    nt_set = set(['A','G','C','T'])
    neighbors_list = list()
    for i in range(len(pattern)):
        nt_remainder_set = nt_set - set(pattern[i])
        for nt in nt_remainder_set:
            neighbors_pattern = pattern[:i] + nt + pattern[i+1:]
            neighbors_list.append(neighbors_pattern)
    return neighbors_list
# %% ----------------------------------------
def neighbors_distance(pattern, d):
    n_list = [pattern] + neighbors(pattern)
    for _ in range(d-1):
        n_list = np.unique(list(map(neighbors, n_list)))
    return n_list
# %% ----------------------------------------
symbol_dict = dict(A=0, C=1, G=2, T=3)

def last_symbol(pattern):
    last_sb = pattern[-1]
    return last_sb

def prefix(pattern):
    prefix = pattern[:-1]
    return prefix

def pattern2num_recur(pattern):
    if pattern == '':
        return 0
    symbol = last_symbol(pattern)
    _prefix = prefix(pattern)
    return 4*pattern2num_recur(_prefix) + symbol_dict[symbol]
# %% ----------------------------------------
def num2symbol(index):
    index_dict = {v:k for k,v in symbol_dict.items()}
    symbol = index_dict[index]
    return symbol

def num2pattern_recur(index, k):
    if k == 1:
        symbol = num2symbol(index)
        return symbol
    prefix_index = index // 4
    remainder = index % 4
    # print(f"k = {k} : prefix_index={prefix_index},remainder={remainder}")
    symbol = num2symbol(remainder)
    prefix_pattern = num2pattern_recur(prefix_index, k-1)
    pattern = prefix_pattern + symbol
    return pattern
# %%
def computing_frequency_w_mismatches(text, k, d):
    freq_arr = np.zeros(4**k)
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        neighbors_list = neighbors_distance(pattern, d)
        for neighbor_pattern in neighbors_list:
            idx = pattern2num_recur(neighbor_pattern)
            freq_arr[idx] += 1
    return freq_arr      
# %% ----------------------------------------
import typing as tp
class MostFrequentPattern(tp.NamedTuple):
    the_most_frequent_patterns:list
    max_frequency:int

def find_the_most_frequent_pattern_w_mismatch(text, k, d):
    freq_arr = computing_frequency_w_mismatches(text, k, d)
    max_freq = np.max(freq_arr)
    idx_arr_max_freq = np.where( freq_arr == max_freq)[0]

    the_most_frequent_patterns = list()
    for idx in idx_arr_max_freq:
        pattern = num2pattern_recur(idx, k)
        the_most_frequent_patterns.append(pattern)

    return MostFrequentPattern(the_most_frequent_patterns,max_freq)
# %% ----------------------------------------
