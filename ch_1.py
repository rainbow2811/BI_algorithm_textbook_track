# %% ----------------------------------------
import numpy as np
# %%  Q) ch1_A
def pattern_count(text, pattern):
    count = 0
    k = len(pattern)
    seq_len = len(text)
    for i in range(seq_len-k+1):
        if text[i:i+k] == pattern:
            count += 1
    return count
# %%
def frequent_words(text, k):
    freq_pattern_set = set()
    count_list = list()
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        count_num = pattern_count(text, pattern)
        count_list.append(count_num)
    max_count = np.max(count_list)
    for i in range(len(text)-k):
        if count_list[i] == max_count:
            pattern = text[i:i+k]
            freq_pattern_set.add(pattern)
    return freq_pattern_set
# %%
def frequent_words_dict(text, k):
    count_dict = dict()
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        count_dict[pattern] = count_dict.get(pattern, 0) + 1
        # count_dict[pattern] = count_dict.get(pattern, []) + i
    return sorted(count_dict.items(), key = lambda x:x[1], reverse=True)
    # return sorted(count_dict.items(),key=lambda x:len(x[1]),reverse=True)
# %% ----------------------------------------
actg = list('ACGT')
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
base = np.array(list('ACGT'))
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
def pattern2num(pattern):
    k = len(pattern)
    kmer_tot_arr = kmer_base(k)
    idx = np.where(kmer_tot_arr == pattern)
    return idx
# %%
def num2pattern(index, k):
    kmer_tot_list = kmer_base(k)
    pattern = kmer_tot_list[index]
    return pattern
# %% ----------------------------------------
def last_symbol(pattern):
    last_sb = pattern[-1]
    return last_sb
# %% ----------------------------------------
def prefix(pattern):
    prefix = pattern[:-1]
    return prefix
# %% ----------------------------------------
symbol_dict = dict(A=0, C=1, G=2, T=3)
index_dict = {v:k for k,v in symbol_dict.items()}
# %%   # 1L
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
# %%
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
# %%  Q) 1K  
def computing_frequencies(text, k):
    frequency_arr = np.zeros(4**k)
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        idx = pattern2num_recur(pattern)
        # print(f"[idx] pattern : [{idx}] {pattern}")
        frequency_arr[idx] += 1
    return frequency_arr
# %%
text = 'AAGCAAAGGTGGG'*100
k = 11
# %% ----------------------------------------
def faster_frequent_words(text, k):
    frequent_patterns = list()
    frequency_arr = computing_frequencies(text,k)
    max_count = np.max(frequency_arr)
    for i in range(4**k):
        if frequency_arr[i] == max_count:
            pattern = num2pattern_recur(i, k)
            frequent_patterns.append(pattern)
    return frequent_patterns
# %%
def faster_frequent_words_np(text, k):
    frequent_patterns = list()
    frequency_arr = computing_frequencies(text,k)
    max_count = np.max(frequency_arr)
    np_idx_list = np.where(frequency_arr == max_count)[0]
    for i in np_idx_list:
        pattern = num2pattern_recur(i,k)
        frequent_patterns.append(pattern)
    return frequent_patterns
# %% ----------------------------------------
faster_frequent_words(text,k)
# %%
%timeit faster_frequent_words(text,k)
# %%
%timeit faster_frequent_words_np(text,k)
# %%
def finding_frequent_words_by_sorting(text,k):
    frequent_patterns = list()
    N_kmer = len(text)-k+1
    idx_arr = np.empty(shape=N_kmer)
    count_arr = np.ones(shape=N_kmer)
    for i in range(N_kmer):
        pattern = text[i:i+k]
        idx_arr[i] = pattern2num_recur(pattern)
        sorted_idx_arr = np.sort(idx_arr)
        for i in range(N_kmer):
            if sorted_idx_arr[i] == sorted_idx_arr[i-1]:
                count_arr[i] = count_arr[i-1] + 1
        