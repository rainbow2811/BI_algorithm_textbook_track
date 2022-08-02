# %% ----------------------------------------
import numpy as np
# %%
text = 'ACAACTATGCATACTATCGGGAACTATCCT'
pattern = 'ACTAT'
# %%
def pattern_count(text, pattern):
    count = 0
    k = len(pattern)
    seq_len = len(text)
    for i in range(seq_len - k):
        if text[i:i+k] == pattern:
            count += 1
    return count
# %% ----------------------------------------
pattern_count(text, pattern)
# %%
def frequent_words(text, k):
    freq_pattern_set = set()
    count_list = list()
    for i in range(len(text)-k):
        pattern = text[i:i+k]
        count_num = pattern_count(text, pattern)
        count_list.append(count_num)
    max_count = np.max(count_list)
    for i in range(len(text)-k):
        if count_list[i] == max_count:
            pattern = text[i:i+k]
            freq_pattern_set.add(pattern)
    return freq_pattern_set
# %% ----------------------------------------
%timeit frequent_words(text, 3)
# %%
frequent_words(text, 4)
# %%
def frequent_words_dict(text, k):
    count_dict = dict()
    for i in range(len(text)-k):
        pattern = text[i:i+k]
        count_dict[pattern] = count_dict.get(pattern, 0) + 1
        # count_dict[pattern] = count_dict.get(pattern, []) + i
    return sorted(count_dict.items(), key = lambda x:x[1], reverse=True)
    # return sorted(count_dict.items(),key=lambda x:len(x[1]),reverse=True)
# %%
frequent_words_dict(text,3)
# %% ----------------------------------------
def patt_to_num(pattern):
    pass
# %% ----------------------------------------
actg = list('ACTG')
actg
# %%
def two_mer_pattern():
    two_mer_list = list()
    for e1 in actg:
        for e2 in actg:
            two_mer = e1 + e2
            two_mer_list.append(two_mer)
    return two_mer_list
# %%
two_mer_list = two_mer_pattern()
# %%
two_mer_list
# %%
np_actg = np.array(actg)
# %%
np.char.add(np_actg.reshape((4,1)),np_actg.reshape((1,4)))
# %%
# broadcasting
arr_3d = np.arange(24).reshape((3,4,2))
arr_2d = np.arange(8).reshape((4,2))
arr_3d + arr_2d
# %%
base = np.array(list('ACTG'))
base, base.shape
# %% ----------------------------------------
def kmer_base(k):
    shape_list = [[1]*i+[4]+[1]*(k-i-1) for i in range(k)]
    result = base.reshape(shape_list[0])
    for i in range(1,k):
        result = np.char.add(result, base.reshape(shape_list[i]))
    # return dict(zip(result.flatten(),np.arange(4**k)))
    return result.flatten()
# %%
kmer_base(3)
# %%
def pattern2num(pattern):
    k = len(pattern)
    kmer_tot_arr = kmer_base(k)
    idx = np.where(kmer_tot_arr == pattern)
    return idx
# %% ----------------------------------------
idx = pattern2num('ATGCAA')
# %%
idx
# %%
def num2pattern(index, k):
    kmer_tot_list = kmer_base(k)
    pattern = kmer_tot_list[index]
    return pattern
# %%
num2pattern(5437, 7)
# %%
num2pattern(5437, 8)
# %%
num2pattern(720, 6)
# %% ----------------------------------------
def last_symbol(pattern):
    last_sb = pattern[-1]
    return last_sb
# %% ----------------------------------------
def prefix(pattern):
    prefix = pattern[:-1]
    return prefix
# %% ----------------------------------------
def symbol2num(symbol):
    sb_dict = dict(A=0, C=1, G=2, T=3)
    return sb_dict[symbol]
# %%
def pattern2num_recur(pattern):
    if pattern == '':
        return 0
    symbol = last_symbol(pattern)
    prefix = prefix(pattern)
    return pattern2num_recur(prefix) + symbol2num(symbol)