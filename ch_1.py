# %% ----------------------------------------
import numpy as np
from collections import Counter
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
[2]*3
# %% ----------------------------------------
def make_shape(k,i):
    # ex) (k=3, i=0) --> [4,1,1] , (k=3, i=1) --> [1,4,1] , (k=3, i=2) --> [1,1,4]
    shape = np.ones(k)
    shape[i] = 4    # base.shape = (4,)이므로 4로 정함. base = array(['A', 'C', 'G', 'T'])
    return shape
    # return [1]*i+[4]+[1]*(k-i-1) 

def kmer_base(k):
    shape_list = [make_shape(k,i) for i in range(k)]  # array(['A', 'C', 'G', 'T'])를 k 숫자 만큼 shape이 다른 배열들(shape_list)을 만듬
    result = base.reshape(shape_list[0])
    for i in range(1,k):   # 1부터 시작하는 이유는 이미 shape_list[0]은 result로 만들었기 때문
        result = np.char.add(result, base.reshape(shape_list[i]))  #  하나씩 broadcasting해서 kmer 자리수 만큼 늘림(2mer -> 3mer -> 4mer ->...)
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
# %%   Q) 1L
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
# %timeit faster_frequent_words(text,k)
# %timeit faster_frequent_words_np(text,k)
# %%
def finding_frequent_words_by_sorting(text,k):
    freq_pattern_list = list()
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

    max_count = np.max(count_arr) # for

    for i in range(N_kmer): # for
        if count_arr[i] == max_count:
            freq_pattern = num2pattern_recur(sorted_idx_arr[i], k)
            freq_pattern_list.append(freq_pattern)
            
    return freq_pattern_list
# %% ----------------------------------------
text = 'AAGCAAAGGTGGG'
k = 2
# %%
finding_frequent_words_by_sorting(text, k)
# %% ----------------------------------------
def finding_frequent_words_by_sorting_nr(text,k):
    freq_pattern_list = list()
    N_kmer = len(text)-k+1
    idx_arr = np.empty(shape=N_kmer)
    count_arr = np.ones(shape=N_kmer)

    for i in range(N_kmer):
        pattern = text[i:i+k]
        idx_arr[i] = pattern2num_recur(pattern)

    sorted_idx_arr = np.sort(idx_arr)
    max_count = count_arr[0]
    for i in range(1, N_kmer):
        if sorted_idx_arr[i] == sorted_idx_arr[i-1]:
            count_arr[i] = count_arr[i-1] + 1
        if max_count < count_arr[i]:
            max_count = count_arr[i]
    freq_idx_list = np.where(count_arr == max_count)[0]
    freq_sorted_idx_list = sorted_idx_arr[freq_idx_list]
    for i in freq_sorted_idx_list:
        pattern = num2pattern_recur(i,k)
        freq_pattern_list.append(pattern)

    return freq_pattern_list
# %%
# %timeit finding_frequent_words_by_sorting(text, k)
# # %% ----------------------------------------
# %timeit finding_frequent_words_by_sorting_nr(text, k)
# %%
def finding_frequent_words_by_sorting_ph(text,k):
    freq_pattern_list = list()
    N_kmer = len(text)-k+1
    idx_arr = np.empty(shape=N_kmer)
    count_arr = np.ones(shape=N_kmer)

    for i in range(N_kmer):
        pattern = text[i:i+k]
        idx_arr[i] = pattern2num_recur(pattern)

    sorted_idx_arr = np.sort(idx_arr)
    max_count = count_arr[0]  # 1
    max_freq_list = list()
    for i in range(1, N_kmer):
        print(f"index = {i}")
        if sorted_idx_arr[i] == sorted_idx_arr[i-1]:
            count_arr[i] = count_arr[i-1] + 1
            print("add count")
        elif max_count == count_arr[i-1]:
            max_freq_list.append(sorted_idx_arr[i-1])
            print("append max_freq_list")
        elif max_count < count_arr[i-1]:
            max_count = count_arr[i-1]
            max_freq_list = [sorted_idx_arr[i-1]]
            print("update new max_count")
            print("remove old max_freq_list")
            print("make new max_freq_list")
        print(f"count_arr = {count_arr}")
        print(f"max_count = {max_count}")
        print(f"max_freq_list = {max_freq_list}")

    for i in max_freq_list:
        pattern = num2pattern_recur(i,k)
        freq_pattern_list.append(pattern)
    return freq_pattern_list
# %%
# %timeit finding_frequent_words_by_sorting_ph(text, k)
# %%
finding_frequent_words_by_sorting_ph(text, k)
# %%
finding_frequent_words_by_sorting_nr(text, k)
# %% 1C
comp_dict = dict(A='T', C='G', T='A', G='C')
comp_dict
# %%
input_seq = 'AGTCGCATAGT'
# %%
import typing as tp    
class Seq(tp.NamedTuple):
    comp_seq:str
    ori_seq:str
    RNA_seq:bool

    # def __init__(self, comp_seq, ori_seq, RNA_seq):
    #     self.comp_seq = comp_seq
    #     self.ori_seq = ori_seq
    #     self.RNA_seq = RNA_seq
    
    # def __str__(self):
    #     return f"comp_seq = {self.comp_seq}, ori_seq = {self.ori_seq}, RNA_seq = {self.RNA_seq}"

    # def __repr__(self):
    #     return f"comp_seq = {self.comp_seq}, ori_seq = {self.ori_seq}, RNA_seq = {self.RNA_seq}"
# %% ----------------------------------------
def get_comp_seq(input_seq, reverse = True, RNA_seq = False):
    if RNA_seq:
        comp_dict = dict(A='U', C='G', T='A', G='C')
    else:
        comp_dict = dict(A='T', C='G', T='A', G='C')

    input_seq_list = list(input_seq)
    comp_seq = list(map(comp_dict.get, input_seq_list))
    # comp_seq = list(map(lambda x: comp_dict[x], input_seq_list))
    
    # if RNA_seq:
    #     comp_seq = [ ('U' if nt == 'T' else nt) for nt in comp_seq ]

    if reverse:
        comp_seq = ''.join(reversed(comp_seq))
    else:
        comp_seq = ''.join(comp_seq)

    # return comp_seq
    return Seq(comp_seq,input_seq,RNA_seq)
# %%
# get_comp_seq(input_seq, reverse=False), get_comp_seq(input_seq, reverse=True)
# # %%
# comp_dict = dict(A='T', C='G', T='A', G='C')
# # %%
# count_dict = dict()
# for n, k in enumerate(input_seq):
#     count_dict[k] = count_dict.get(k,list())+[n]
# count_dict
# # %%
# comp_dict = dict(A='T', C='G', T='A', G='C')
# input_seq_list = list(input_seq)
# comp_seq = list(map(comp_dict.get, input_seq_list))
# # %%
# comp_seq
# %% ----------------------------------------
# [ ('U' if nt == 'T' else nt) for nt in comp_seq ]
# # %%
# get_comp_seq(input_seq, reverse=False, RNA_seq=False)
# # %%
# get_comp_seq(input_seq, reverse=False, RNA_seq=True)
# # %%
# get_comp_seq(input_seq, reverse=True, RNA_seq=False)
# # %%
# get_comp_seq(input_seq, reverse=True, RNA_seq=True)
# # %%
# get_comp_seq(input_seq, reverse=True, RNA_seq=True)
# %% ----------------------------------------
class CompSeq:
    def __init__(self):
        self.dna_comp_dict = dict(A='U', C='G', T='A', G='C')
        self.rna_comp_dict = dict(A='T', C='G', T='A', G='C')
    
    def __call__(self, input_seq, reverse = True, RNA_seq = False):
        if RNA_seq:
            comp_dict = self.rna_comp_dict
        else:
            comp_dict = self.dna_comp_dict
        
        input_seq_list = list(input_seq)
        comp_seq = list(map(comp_dict.get, input_seq_list))

        if reverse:
            comp_seq = ''.join(reversed(comp_seq))
        else:
            comp_seq = ''.join(comp_seq)

        return Seq(comp_seq,input_seq,RNA_seq)
# %%
_get_comp_seq = CompSeq()
# %%
_get_comp_seq(input_seq,True,False)
# %%     Q) 1E
def clump_finding(genome_seq, k, L, t):
    freq_patterns = list()
    clump_arr = np.zeros(shape=(4**k))       # freq_arr와 동일한 형태 (가능한 모든 kmer 경우의 수)
    for i in range(len(genome_seq)-L+1):
        window = genome_seq[i:i+L]
        freq_arr = computing_frequencies(window, k)
        idx_clump = np.where(freq_arr >= t)[0]
        if len(idx_clump) != 0:
            clump_arr[idx_clump] = 1         # clump count 누적 없이 존재 유무만 체크 (1=있음, 0=없음)
    idx_arr = np.where(clump_arr == 1)[0]    # count 누적 안되게 존재 유뮤만 알 수 있는 idx_arr를 만들기 위함
    for i in idx_arr:
        pattern = num2pattern_recur(i, k)
        freq_patterns.append(pattern)
    return freq_patterns
# %% ----------------------------------------
genome_seq = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
k = 5
L = 75
t = 4
# clump_finding(genome_seq, 5, 75, 4)
# %%   Q) 1E (ver.코드 정리)

def find_clump_idxs(genome_seq, k, L, t):
    # clump_arr = np.zeros(shape=(4**k))    
    clump_idxs = np.array([])
    for i in range(len(genome_seq)-L+1):    
        window = genome_seq[i:i+L]
        freq_arr = computing_frequencies(window, k)
        idx_clump = np.where(freq_arr >= t)[0]
        clump_idxs = np.concatenate([clump_idxs, idx_clump])
    return np.unique(clump_idxs)

    #     if len(idx_clump) != 0:
    #         clump_arr[idx_clump] = 1     
    # idx_arr = np.where(clump_arr == 1)[0]
    # return idx_arr

def num_arr2pattern_list(clump_idx_arr, k):
    freq_patterns = list()
    for i in clump_idx_arr:
        pattern = num2pattern_recur(i, k)
        freq_patterns.append(pattern)
    return freq_patterns

def find_clump_patterns(genome_seq, k, L, t):
    clump_idx_arr = find_clump_idxs(genome_seq, k, L, t)
    freq_patterns = num_arr2pattern_list(clump_idx_arr, k)
    return freq_patterns 
# %% ----------------------------------------
def find_clump_patterns_better(genome_seq, k, L, t):
    clump_arr = np.zeros(shape=(4**k)) 
    window = genome_seq[0:L]
    freq_arr = computing_frequencies(window, k)
    idx_arr = np.where(freq_arr >= t)[0]
    clump_arr[idx_arr] = 1
    for i in range(1, len(genome_seq)-L+1):
        first_pattern = genome_seq[i-1:k]
        first_patt_idx = pattern2num_recur(first_pattern)
        freq_arr[first_patt_idx] -= 1

        window = genome_seq[i:i+L]   
        last_pattern = window[-k:]   # window의 끝에서부터 k개 선택 
        last_patt_idx = pattern2num_recur(last_pattern) 
        freq_arr[last_patt_idx] += 1

        idx_arr = np.where(freq_arr >= t)[0]
        clump_arr[idx_arr] = 1

    clump_idx_arr = np.where(clump_arr == 1)[0]
    freq_patterns_list = list()
    for idx in clump_idx_arr:
        pattern = num2pattern_recur(idx, k)
        freq_patterns_list.append(pattern)
    return freq_patterns_list
# %%
%timeit find_clump_patterns_better(genome_seq, k, L, t)
# %%
%timeit find_clump_patterns(genome_seq, k, L, t)
# %% ----------------------------------------
4000/322
# %%   Q) 1G
genome_seq = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
# %%
# np.fromstring(genome_seq,dtype=np.uint8)
# %%
# %timeit aa = np.fromstring(genome_seq,dtype=np.uint8) % 5 % 4
# %%
num_dict = dict(A=0,T=0,C=-1,G=1)
seq_list = list(genome_seq)
# %% ----------------------------------------
# %timeit aa = np.cumsum(list(map(num_dict.get, list(genome_seq))))
# %% ----------------------------------------
nt2num_seq = list(map(num_dict.get, list(genome_seq)))
# %% ----------------------------------------
nt2num_seq = np.insert(nt2num_seq, obj=0, values=0)
# %% ----------------------------------------
arr = np.cumsum(nt2num_seq)
arr
# %% ----------------------------------------
min_idx_list = np.where(arr == np.min(arr))
# %%
from bokeh.plotting import figure, save
from bokeh.io import output_notebook, show
output_notebook()
# %%
x = np.arange(len(arr))
p = figure(title = 'genome sequence asymmetric diagram', x_axis_label = 'location', y_axis_label = 'G count-C count')
p.line(x,arr, line_width = 2)
show(p)
# %%
min_idx_list
# %% ----------------------------------------
class SeqAsymmetry:
    def __init__(self, genome_seq):
        # get cum_arr and min/max idx list
        num_dict = dict(A=0,T=0,C=-1,G=1)
        nt2num_seq = list(map(num_dict.get, list(genome_seq))) # nt를 숫자로 변환
        nt2num_seq_add_0 = np.insert(nt2num_seq, obj=0, values=0) # 첫번째 dix에 0값 추가
        cum_arr = np.cumsum(nt2num_seq_add_0)
        min_idx_list = np.where(cum_arr == np.min(cum_arr))
        max_idx_list = np.where(cum_arr == np.max(cum_arr))

        # diagram
        x = np.arange(len(arr))
        p = figure(title = 'genome sequence asymmetric diagram', x_axis_label = 'location', y_axis_label = 'G count-C count')
        p.line(x,arr, line_width = 2)

        self.diagram = p
        self.ori_position = min_idx_list
        self.ter_position = max_idx_list
    
    def show_fig(self):
        show(self.diagram)
# %%
seq_asymmetry = SeqAsymmetry(genome_seq)
# %%
seq_asymmetry.ori_position
# %%
seq_asymmetry.ter_position
# %%
seq_asymmetry.show_fig()
# %%
