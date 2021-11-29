#! /usr/bin/env python
# test file for importing created module pybind11_mcs.

import pybind11_mcs

print(pybind11_mcs.Original_Basis(8))

pybind11_mcs.MCS()