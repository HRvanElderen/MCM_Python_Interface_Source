# tasks.py
# compiles pybind11 module for functions included in wrapper.cpp.
import invoke

invoke.run(
    "g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3.8 -m pybind11 --includes) wrapper.cpp MinCompSpin-main/*.cpp -o mcm_interface$(python3.8-config --extension-suffix) -undefined dynamic_lookup"
)