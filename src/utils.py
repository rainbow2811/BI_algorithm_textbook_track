# %% ----------------------------------------
import typing as tp
import numpy as np
from src import typing as ptp
# %% ----------------------------------------
def check_all_True_or_not(
    data: tp.List[tp.Set[ptp.Kmer]], 
    kmer_neighbor: ptp.Kmer
    ) -> bool:
    check = list()
    for set_ in data:
        check.append(int(kmer_neighbor in set_))
    result_check_all_True =  np.prod(check)
    return result_check_all_True
# %%
