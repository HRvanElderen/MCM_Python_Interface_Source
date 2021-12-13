Repo containing the code for an interface between C++ and Python for a community detection algorithm.

The library uses Python 3.8.

**To build:** `python build.py`

after this the package can be imported into your python file with the name: `mcm_interface`


# Usage 

## MCM

currently the `MCM` function is the only tested function. This function runs the entire program. The function takes in four arguments:
- `int` `n` the number of spin variables.
- `list` `Basis_Choice` This can be defined manually or with the function `Original_Basis(n)` which takes as input the number of spin variables.
- `list` `MCM_Choice` Specify the MCM. For example, `MCM_Choice = [384, 64, 32, 16, 8, 4, 2, 1]` defines an MCM with `8` independent parts, based on `n=9` spins.
- `string` `inputfile` The datafile to test.
