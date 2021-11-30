#! /usr/bin/env python
# test file for importing created module pybind11_mcs.
import numpy as np
import mcs_interface

MCM_Choice = [384, 64, 32, 16, 8, 4, 2, 1]
Basis_Choice = [3, 5, 9, 48, 65, 129, 272, 81, 1]


N = np.uint(0) 
n = np.uint32(9)
mcs_interface.MCS(n, Basis_Choice, MCM_Choice, "INPUT/SCOTUS_n9_N895.dat")
a = mcs_interface.read_datafile(N, n, "INPUT/SCOTUS_n9_N895.dat")

types1 = [type(k) for k in a.keys()]
print(types1)
print(type(a[63]))
a[63]=67.8
print(type(a[63]))

