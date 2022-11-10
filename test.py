# %% ----------------------------------------
import src
import numpy as np

import ch_2
# %%
ch_2.__name__
# %% ----------------------------------------
k = 3
d = 1
# %% ----------------------------------------
dnas = """
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG
"""
dnas
# %%
dna_list = ch_2.get_dna_list_from_raw_data(dna_data=dnas)
# %%
dna_list[0]
# %%
ch_2.get_kmer_neighbors_set_from_one_dna_string(dna_string=dna_list[0], k=3, d=1)
# %%
