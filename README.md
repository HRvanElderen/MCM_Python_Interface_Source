# MCM_Python_Interface_Source

Python package that interfaces with the C++ implementation of the MCM algorithm for community detection.

**The C++ implementation:** https://github.com/clelidm/MinCompSpin

The library uses Python 3.

## Requirements

Pybind11: `pip install pybind11`

## Building

**To build:**   `g++ -O3 -Wall -shared -std=c++11 -fPIC -undefined dynamic_lookup ``(python3 -m pybind11 --includes)`` wrapper.cpp MinCompSpin/*.cpp -o mcm_interface``( python3-config --extension-suffix)`

or run the following command with `invoke` installed:

`python build.py`

after this the package can be imported into your python file with the name: `mcm_interface`. Two tutorials notebooks are added. A long one containing all functions available in the package. And a short one with a bit of data analysis.

