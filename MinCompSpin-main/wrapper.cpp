/* pybind11_wrapper.cpp 
 * Includes all function that need to be translated from cpp to python.
 * Called from tasks.py to map functions.
 */

#include <pybind11/pybind11.h>
#include "pybind11/stl.h"   // support for standard library.
#include "main.cpp" // cpp source.
#include <vector>

PYBIND11_MODULE(pybind11_mcs, m) {
    m.doc() = "pybind11 pybind11_mcs plugin"; // Optional module docstring
    m.def("Original_Basis", &Original_Basis);
    m.def("Read_BasisOp_BinaryRepresentation", &Read_BasisOp_BinaryRepresentation);
    m.def("Read_BasisOp_IntegerRepresentation", &Read_BasisOp_IntegerRepresentation);
    m.def("PrintTerm_Basis", &PrintTerm_Basis);
    m.def("read_datafile", &read_datafile);
    m.def("build_Kset", &build_Kset);
    m.def("LogL_SubCM", &LogL_SubCM);
    m.def("LogE_SubCM", &LogE_SubCM);
    m.def("GeomComplexity_SubCM", &GeomComplexity_SubCM);
    m.def("ParamComplexity_SubCM", &ParamComplexity_SubCM);
    m.def("LogL_CM", &LogL_CM);
    m.def("LogL_MCM", &LogL_MCM);
    m.def("LogE_MCM", &LogE_MCM);
    m.def("Complexity_MCM", &Complexity_MCM);
    m.def("Create_MCM", &Create_MCM);
    m.def("Read_MCMParts_BinaryRepresentation", &Read_MCMParts_BinaryRepresentation);
    m.def("check_partition", &check_partition);
    m.def("PrintTerminal_MCM_Info", &PrintTerminal_MCM_Info);
    m.def("PrintInfo_All_Indep_Models", &PrintInfo_All_Indep_Models);
    m.def("PrintInfo_All_SubComplete_Models", &PrintInfo_All_SubComplete_Models);
    m.def("MCM_GivenRank_r", &MCM_GivenRank_r);
    m.def("MCM_AllRank_SmallerThan_r_Ordered", &MCM_AllRank_SmallerThan_r_Ordered);
    m.def("MCM_AllRank_SmallerThan_r_nonOrdered", &MCM_AllRank_SmallerThan_r_nonOrdered);
    m.def("MCM", &MCM, "spin model function", pybind11::arg("n") = 9, 
                                              pybind11::arg("datafilename") = "INPUT/Dataset_Shapes_n9_N1e5.dat", 
                                              pybind11::arg("basis_IntegerRepresentation_filename") = "INPUT/Dataset_Shapes_n9_Basis_Integer.dat", 
                                              pybind11::arg("basis_BinaryRepresentation_filename") = "INPUT/Dataset_Shapes_n9_Basis_Binary.dat", 
                                              pybind11::arg("OUTPUT_directory") = "OUTPUT/", 
                                              pybind11::arg("Basis_Choice") =  Def_Basis_Choice,
                                              pybind11::arg("MCM_Choice") =  Def_MCM_Choice);
}