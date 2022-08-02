# %% ----------------------------------------
import numpy as np
# %% ----------------------------------------
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
# %%
idx = pattern2num('ATGCAA')
# %%
six_mer = kmer_base(6)
# %%
six_mer[idx]
# %%
idx_2 = ([3,6,9,2], [1,2,3,4], [7,9,0,1])
# %% ----------------------------------------
list(zip(*idx_2))
# %% ----------------------------------------
c1, c2, c3, c4 = zip_ab
# %%
list(zip(*zip_ab))
# %%
def func(*arg, **kwarg):
    return dict(arg=arg, kwarg=kwarg)
# %%
func(1,2,3,x=1,y=2)
# %%
def func2(x,y):
    return x+y
# %% ----------------------------------------
def find_max(*arg):
    return max(arg)
# %%
find_max(4,6)
# %%
func2(1,2)
# %%
func2(**{'x':2,'y':3})
# %%
