//#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <bitset> 
#include <set>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/*********************   OPERATOR EVALUATED on a STATE   **********************/
/***********************   depend on the convention    ************************/
/******************************************************************************/

/****************************   !! Convention !!   ****************************/
/****     For Ising:      {1 in data <--> -1 in model}                     ****/
/****                     {0 in data <--> +1 in model}                     ****/
/******************************************************************************/

int Op_Ising(uint32_t Op, uint32_t state)       // Convention {-1;1} 
  {  return ( (bitset<n>(Op & state).count())%2 )?(-1):1;   } 

int Op_ACE(uint32_t Op, uint32_t state)         // Convention {0;1} 
  {  return ( (Op & state) == Op );   } 

// return the contribution to the Energy H(s) associated to this Interaction
double Ising(Interaction I, uint32_t state, int *Op_s)  // Convention {-1;1} 
{
  *Op_s = Op_Ising(I.Op, state);
  return I.g * (*Op_s);
}

double ACE(Interaction I, uint32_t state, int *Op_s)    // Convention {0;1} 
{
  *Op_s = Op_ACE(I.Op, state);
  return I.g * (*Op_s);
}

/******************************************************************************/
/************* PARTITION FUNCTION and PROBA of ALL STATES *********************/
/******************************************************************************/
// return the probability (not normalised) of all the states and compute the partition function

double* Probability_AllStates(double MConvention(Interaction, uint32_t, int*), list<Interaction> list_I, double *Z)  // Convention {-1;1} 
//!! Convention !!:  {1 in data <--> -1 in model}  and  {0 in data <--> 1 in model} 
{
  double H = 0; // -energy of the state
  int Op_s = 1; // value of the operator for the state s ; \in {-1; 1}

  list<Interaction>::iterator I;
  double* all_P = (double*)malloc((NOp_tot+1)*sizeof(double));

  (*Z) = 0;

  for (uint32_t state = 0; state <= NOp_tot; state++)
  {
    H=0;  // here H is (H = -Hamiltonian)
    for (I = list_I.begin(); I != list_I.end(); I++)
    {
      H += MConvention((*I), state, &Op_s);  //H += (*I).g * Op_s;
      //cout << bitset<n>((*I).Op) << "\t" << bitset<n>(state) << "\t" << Op_s;
      //cout << "\t " << (*I).g << "\t " << H << endl;
    }
    all_P[state] = exp(H);
    (*Z) += all_P[state];
    //cout << all_P[state] << "\t" << (*Z) << endl;
  }

  return all_P;
}

/******************************************************************************/
/*************  PRINT PROBA  --  ALL STATES  and  VS DATA  ********************/
/******************************************************************************/
void print_proba_state(double* P, double Z, string OUTPUTfilename)
{
  double Ztest = 0;
  fstream file(OUTPUTfilename.c_str(), ios::out);

  file << "#Total number of accessible states = " << NOp_tot << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state \t #4: state_bin" << endl;

  for (uint32_t state = 0; state <= NOp_tot; state++)
  {
    file << state << "\t" << bitset<n>(state).count() << "\t" << P[state]/Z << "\t" << bitset<n>(state) << endl;
    Ztest += P[state];
  }  
  cout << "Z = " << Z << "\t Ztest = " << Ztest << endl;

  file.close();
}

void Model_VS_data_ProbaSpace(double *P, double Z, map<uint32_t, unsigned int> Nset_test, unsigned int N_test, string outputfilename)
{
  map<uint32_t, unsigned int>::iterator it;
  uint32_t state = 1;

  fstream file(outputfilename.c_str(), ios::out);
  file << "#N_test = " << N_test << endl;
  file << "#Total number of accessible states = " << NOp_tot << endl;
  file << "#Number of visited states, Nset_test.size() = " << Nset_test.size() << endl;
  file << "#" << endl;
  file << "#1: nb of '1's in the state \t #2: Empirical pba state \t #3: Model pba state \t #4: state" << endl;

  for (it = Nset_test.begin(); it!=Nset_test.end(); ++it)
  {
    state = it->first;
    file << bitset<n>(state).count() << "\t" << (it->second) / (float) N_test  << "\t" << P[state]/Z; // << endl;
    file << "  \t " << state << " \t " << bitset<n>(state) << endl;
  }

  file.close();
}
/******************************************************************************/
/****************************  LogLikelihood  *********************************/
/******************************************************************************/
double LogLikelihood(double *P, double Z, map<uint32_t, unsigned int> Nset)
{
  double LogL = 0;
  //double LogL2 = 0;  unsigned int N = 0;
  map<uint32_t, unsigned int>::iterator it;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    LogL += (it->second) * log( P[(it->first)] / Z );
  //  LogL2 += (it->second) * log( P[(it->first)] );
  //  N += (it->second);
  }
  //LogL2 -= N * log(Z);

  //cout << "N = " << N << "\t LogL = " << LogL2 << endl;

  return LogL;
}

// Given the value of Kset --> return the Loglikelihood
double LogLikelihood_CompleteModel(map<uint32_t, unsigned int> Nset, unsigned int N)  //
{
  map<uint32_t, unsigned int>::iterator it;
  double L = 0.;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    Ks = (it->second);
    if (Ks == 0) {cout << "problem Ks = 0 for mu = " << (it->first) << endl; }
    L += (Ks * log( ((double) Ks) / ((double) N) ) );
    Ncontrol += Ks;
  }
  if (Ncontrol != N) { cout << "Error Likelihood function: Ncontrol != N" << endl;  }

  return L;
}
/******************************************************************************/
/*********************************  Averages  *********************************/
/******************************************************************************/
/************************    Model averages all op    *************************/
void Model_averages(double MConvention(Interaction, uint32_t, int*), double *P, double Z, list<Interaction> &list_I, unsigned int N) 
{
  int Op_s = 1; // value of the operator for the state s ; \in {-1; 1}

  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    (*I).av_M = 0.;
    for (uint32_t state = 0; state <= NOp_tot; state++)
    {
      MConvention((*I), state, &Op_s); //Op_s = ( (bitset<n>((*I).Op & state).count())%2 )?(-1):1;
      //cout << bitset<n>((*I).Op) << "\t" << bitset<n>(state) << "\t" << (bitset<n>((*I).Op & state).count())%2 << "\t" << Op_s;

      (*I).av_M += Op_s * P[state];
      //cout << "\t " << (*I).g << "\t " << H << endl;
    }
    //cout << (*I).av_M << endl;
    (*I).av_M = (*I).av_M / Z;
    //cout << (*I).av_M << endl;
  }
}

double* AllOp_av(double *P, double Z, unsigned int N)  // _all_Model_averages in ISING convention
{
  double* m1 = (double *)malloc((NOp_tot+1)*sizeof(double));
  int Op_s = 1;

  for (uint32_t Op = 1; Op <= NOp_tot; Op++)//NOp_tot
  {
    m1[Op] = 0.;
    for (uint32_t state = 0; state <= NOp_tot; state++)
    {
      Op_s = Op_Ising(Op, state);
      m1[Op] += Op_s * P[state];
    }
    m1[Op] = m1[Op] / Z;
  }
  return m1;
}

double* AllOp_m1(double *P, double Z, unsigned int N)  // _all_Model_averages in ISING convention
{
  double* m1 = (double *)malloc((NOp_tot+1)*sizeof(double));

  for (uint32_t Op = 1; Op <= NOp_tot; Op++)//NOp_tot
  {
    m1[Op] = 0.;
    for (uint32_t state = 0; state <= NOp_tot; state++)
    {
      if( (bitset<n>(Op & state).count())%2 ){ m1[Op] += P[state]; }    // convention: {-1;1} <--> {1, 0}
    }
    m1[Op] = m1[Op] / Z;
  }
  return m1;
}
