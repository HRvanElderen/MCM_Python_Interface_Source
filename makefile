interface: 
	g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) wrapper.cpp MinCompSpin/*.cpp -o mcm_interface$(python3 -config --extension-suffix) -undefined dynamic_lookup
