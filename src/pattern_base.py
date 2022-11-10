# %% ----------------------------------------
import numpy as np
# %% ----------------------------------------
def pattern_matrix (text, k):
    """
    text에서 만들 수 있는 모든 kmer 구하기
    
    (예시)
    input = 'ACCGT', 3

    -->
    output = {'ACC', 'CCG', 'CGT'}
    """
    pattern_set = set()
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        pattern_set.add(pattern)
    return pattern_set
# %% ----------------------------------------
def pattern_count(text, pattern):
    """
    text에서 특정 pattern의 개수 구하기

    (예시)
    input = 'ACCGTACCGT', 'CGT'

    -->
    output = 2
    """
    count = 0
    k = len(pattern)
    seq_len = len(text)
    for i in range(seq_len-k+1):
        if text[i:i+k] == pattern:
            count += 1
    return count
# %% ----------------------------------------
def frequent_words(text, k):
    """
    
    """
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
# %% ----------------------------------------
def frequent_words_dict(text, k):
    count_dict = dict()
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        count_dict[pattern] = count_dict.get(pattern, 0) + 1
        # count_dict[pattern] = count_dict.get(pattern, []) + i
    return sorted(count_dict.items(), key = lambda x:x[1], reverse=True)
# %% ----------------------------------------
base = np.array(list('ACGT'))

def make_shape(k,i):
    # ex) (k=3, i=0) --> [4,1,1] , (k=3, i=1) --> [1,4,1] , (k=3, i=2) --> [1,1,4]
    shape = np.ones(k, dtype=np.int8)
    shape[i] = 4    # base.shape = (4,)이므로 4로 정함. base = array(['A', 'C', 'G', 'T'])
    return shape
    # return [1]*i+[4]+[1]*(k-i-1) 

def kmer_base(k):
    """
    모든 가능한 kmer pattern 구하기 

    (예시)
    input k=2 

    -->
    output =
    ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG',
     'GT', 'TA', 'TC', 'TG', 'TT']
    """

    shape_list = [make_shape(k,i) for i in range(k)]  # array(['A', 'C', 'G', 'T'])를 k 숫자 만큼 shape이 다른 배열들(shape_list)을 만듬
    result = base.reshape(shape_list[0])
    for i in range(1,k):   # 1부터 시작하는 이유는 이미 shape_list[0]은 result로 만들었기 때문
        result = np.char.add(result, base.reshape(shape_list[i]))  #  하나씩 broadcasting해서 kmer 자리수 만큼 늘림(2mer -> 3mer -> 4mer ->...)
    # return dict(zip(result.flatten(),np.arange(4**k)))
    return result.flatten()
# %% ----------------------------------------
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
