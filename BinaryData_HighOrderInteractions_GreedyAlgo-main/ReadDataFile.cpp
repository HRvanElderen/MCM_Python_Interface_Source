//#include <iostream>
#include <fstream>
#include <map>
#include <bitset> 

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
map<uint32_t, unsigned int> read_datafile(string datafilename, unsigned int *N)   // O(N)  //N = data set size
{
  string line, line2;     uint32_t nb = 0;
  (*N) = 0;            // N = dataset size (global variable)
  cout << endl << "--->> Read \"" << datafilename << "\",\t Build Nset..." << endl;

// ***** data are store in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  ifstream myfile (datafilename.c_str());

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line
      nb = bitset<n>(line2).to_ulong();   //convert string line2 into a binary integer
      Nset[nb] += 1;
      //cout << line << endl;   //cout << nb << " :  " << bitset<n>(nb) << endl;
      (*N)++;
    }
    myfile.close();
  }
  else cout << "Unable to open file"; 

  cout << "\t\t data size, N = " << (*N) << endl;
  cout << "\t\t number of different states, Nset.size() = " << Nset.size() << endl;

  return Nset;
}

/****************    PRINT Nset in file:    ************************/
void read_Nset (map<uint32_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename)
// map.second = nb of time that the state map.first appears in the data set
{
  map<uint32_t, unsigned int>::iterator it;
  int Ncontrol = 0;

  fstream file(OUTPUTfilename.c_str(), ios::out);
  file << "#N = " << N << endl;
  file << "#Total number of accessible states = " << NOp_tot << endl;
  file << "#Number of visited states, Nset.size() = " << Nset.size() << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    file << it->first << ":\t" << bitset<n>(it->first) << " => " << it->second; // << endl;
    file << "  \t  P = " << it->second / (float) N << endl;
    Ncontrol += it->second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }

  file.close();
}

