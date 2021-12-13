#! /usr/bin/env python
# test file for importing created module pybind11_mcs.
import numpy as np
import mcs_interface
import pandas as pd

MCM_Choice = [384, 64, 32, 16, 8, 4, 2, 1]
Basis_Choice = [3, 5, 9, 48, 65, 129, 272, 81, 1]


N = np.uint(0) 
n = 9
#mcs_interface.Original_Basis(n)
mcs_interface.MCS(n, Basis_Choice, MCM_Choice, "INPUT/SCOTUS_n9_N895.dat")
#a = mcs_interface.read_datafile(N, n, "INPUT/SCOTUS_n9_N895.dat")
#Kset = mcs_interface.build_Kset(a, Basis_Choice, n, False)

