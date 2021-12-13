#! /usr/bin/env python
# test file for importing created module pybind11_mcs.
import mcm_interface

MCM_Choice = [384, 64, 32, 16, 8, 4, 2, 1]
Basis_Choice = [3, 5, 9, 48, 65, 129, 272, 81, 1]


n = 9

#mcm_interface.MCM(n, Basis_Choice, MCM_Choice, "INPUT/SCOTUS_n9_N895.dat")
mcm_interface.MCM(n, Basis_Choice, MCM_Choice, "INPUT/Dataset_Shapes_n9_N1e5.dat")


