/* pybind11_wrapper.cpp 
 * Includes all function that need to be translated from cpp to python.
 * Called from tasks.py to map functions.
 */

#include <pybind11/pybind11.h>
#include "pybind11/stl.h"   // support for standard library.
#include "main.cpp" // cpp source.

PYBIND11_MODULE(pybind11_mcs, m) {
    m.doc() = "pybind11 pybind11_mcs plugin"; // Optional module docstring
    m.def("Original_Basis", &Original_Basis);
    m.def("main", &main);
}