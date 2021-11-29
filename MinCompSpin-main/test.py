#! /usr/bin/env python
# test file for importing created module pybind11_mcs.

import pybind11_mcs

MCM_Choice = [384, 64, 32, 16, 8, 4, 2, 1]
Basis_Choice = [3, 5, 9, 48, 65, 129, 272, 81, 1]
n = 9
pybind11_mcs.MCS(n, Basis_Choice, MCM_Choice, "INPUT/Dataset_Shapes_n9_N1e5.dat")
