#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include "support.h"

using namespace std;


/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
map<uint32_t, unsigned int> read_datafile(unsigned int *N, unsigned int n, string filename)    // O(N)  where N = data set size
{
  string line, line2;    
  (*N) = 0;            // N = dataset size
  cout << endl << "--->> Read \"" << filename << "\",\t Build Nset...";

// ***** data are store in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set
  
  ifstream myfile (filename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);
      Nset[stoi(line2, 0, 2)] += 1;
      //cout << nb.to_ulong() << endl;   //cout << nb << " :  " << bitset<n>(nb) << endl;
      (*N)++;
    }
    myfile.close();
  }
  else cout << "Unable to open file"; 

  cout << "\t\t data size N = " << (*N) << endl;

  return Nset;
}


/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a model (defined by the m basis vector) --> return the new m-state (state in the new m-basis)
// Rem: must have m <= n 
uint32_t transform_mu_basis(uint32_t mu, unsigned int n, list<uint32_t> basis)
{
  uint32_t bit_i = 1;
  uint32_t final_mu = 0;

  list<uint32_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    if ( (countSetBits((*phi_i) & mu) % 2) == 1) // odd number of 1, i.e. sig_i = 1
      {
        final_mu += bit_i;
      }
    bit_i = (bit_i << 1);
  }

  return final_mu;
}


/******************************************************************************/
/************************** K_SET *********************************************/
/******************************************************************************/
// Build Kset for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

map<uint32_t, unsigned int> build_Kset(map<uint32_t, unsigned int> Nset, list<uint32_t> Basis, unsigned int n, bool print_bool=false)
// sig_m = sig in the new basis and cut on the m first spins 
// Kset[sig_m] = #of time state mu_m appears in the data set
{
  map<uint32_t, unsigned int>::iterator it;
  map<uint32_t, unsigned int > Kset;

  uint32_t s;        // initial state
  uint32_t sig_m;    // transformed state and to the m first spins

  unsigned int ks=0; // number of time state s appear in the dataset

  cout << endl << "--->> Build Kset..." << endl;

//Build Kset:
  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    s = it->first;       // state s
    ks = it->second;    // # of times s appears in the data set
    sig_m = transform_mu_basis(s, n, Basis);
//    sig_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    if (print_bool)  {  cout << s << ": \t" << int_to_bstring(s, n) << " \t" << sig_m << ": \t" << int_to_bstring(sig_m, n) << endl; }
    Kset[sig_m] += ks;
    //Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
  cout << endl;

  return Kset;
}
