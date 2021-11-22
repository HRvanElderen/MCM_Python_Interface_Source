# tasks.py
import invoke

invoke.run(
    "g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3.7 -m pybind11 --includes) wrapper.cpp main.cpp Data_Manipulation.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp MCM_info.cpp -o pybind11_mcs$(python3.7-config --extension-suffix) -undefined dynamic_lookup"
)

