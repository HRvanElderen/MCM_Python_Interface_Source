# MCM_Python_Interface_Source

Python package that interfaces with the C++ implementation of the MCM algorithm for community detection.

**The C++ implementation:** https://github.com/clelidm/MinCompSpin

The library uses Python 3.7.

## Requirements

Pybind11: `pip install pybind11`

## Building

**To build:**  `g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3.7 -m pybind11 --includes) wrapper.cpp MinCompSpin-main/*.cpp -o mcm_interface$(python3.7-config --extension-suffix) -undefined dynamic_lookup`


after this the package can be imported into your python file with the name: `mcm_interface`
