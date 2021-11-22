// pybind11_wrapper.cpp
#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include "library.h"

PYBIND11_MODULE(pybind11_mcs, m) {
    m.doc() = "pybind11 pybind11_mcs plugin"; // Optional module docstring
    m.def("Original_Basis", &Original_Basis);
}