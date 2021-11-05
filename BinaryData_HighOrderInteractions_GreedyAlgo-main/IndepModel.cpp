//#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <set>
#include <bitset>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/****************************   Create INDEP MODEL   **************************/
/********************************    from Nset   ******************************/
/******************************************************************************/
/*****************************  in INITIAL BASIS  ****************************/
//probability of an operator to be equal to 1 ( = <phi> in the {0,1} representation )
unsigned int k1_op(map<uint32_t, unsigned int> Nset, uint32_t op)  // Complexity = O(|Nset|)
{
  unsigned int k1=0;
  map<uint32_t, unsigned int>::iterator it;  // iterator on Nset

  for (it = Nset.begin(); it!=Nset.end(); ++it)
    {    k1 += (bitset<n>( ((*it).first) & op ).count() % 2)*((*it).second);   }

  return k1;
}
/******************************************************************************/
/**************************     INDEPENDENT MODEL    **************************/
/******************************************************************************/
list<Interaction> IndepModel_fromBasis(vector<pair<Operator, bool> > best_basis)
{
  Interaction I;
  list<Interaction> list_I;

  vector<pair<Operator, bool> >::iterator BasisOp;

  cout << "Fields associated to the independent model: " << endl;
  cout << "# 1:Op\t 2:p1 \t 3:h_{-1,1} \t 4:h_{0,1}" << endl << endl;

  for (BasisOp = best_basis.begin(); BasisOp != best_basis.end(); BasisOp++)
  {
    I.Op = (*BasisOp).first.bin;
    I.k = bitset<n>(I.Op).size();
    I.g_Ising = 0.5*log( (1.-(*BasisOp).first.p1_D) / (*BasisOp).first.p1_D );
    I.g = I.g_Ising;
    I.g_ACE = -2.*I.g_Ising;
    I.av_D = 0; I.av_M = 0;
    list_I.push_back(I); 
    cout << "Op = " << bitset<n>(I.Op) << "\t p1 = " << (*BasisOp).first.p1_D << " \t\t h_{-1,1} = " << I.g_Ising << " \t h_{0,1} = " << I.g_ACE << endl;
  }
  cout << endl;
  return list_I;
}
/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
uint32_t state_new(uint32_t s_old, list<Interaction> basis)
{
  uint32_t s_new = 0;
  uint32_t un_i = un;

  list<Interaction>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    if ( (bitset<n>( (*phi_i).Op & s_old ).count() % 2) == 1) // odd number of 1, i.e. sig_i = 1
      {  s_new += un_i; }
    un_i = (un_i << 1);
  }

  return s_new;
}
/******************************************************************************/
/************    Independent model in the PROBABILITY SPACE:    ***************/
/******************************************************************************/
double proba_fromIndepModel(vector<pair<Operator, bool> > best_basis, uint32_t state)  //for the moment: work only in the (s1, s2, .., sn) basis
{   
  double proba = 1.;
  vector<pair<Operator, bool> >::iterator phi;

  for(phi = best_basis.begin(); phi != best_basis.end(); ++phi)
  {
    if ( (bitset<n>( ((*phi).first.bin) & state ).count() % 2)  == 1 )
      { proba = proba * ((*phi).first.p1_D); }
    else
      { proba = proba * (1-((*phi).first.p1_D)); }
  }

  return proba;
}

void IndepModel_VS_data_ProbaSpace(vector<pair<Operator, bool> > best_basis, map<uint32_t, unsigned int> Nset_test, unsigned int N_test, string OUTPUT_filename)
{
  map<uint32_t, unsigned int>::iterator it;
  uint32_t state = 1;
  double p_data, p_model, sig, r;

  fstream file((OUTPUT_filename + "_ProbaSpace.dat").c_str(), ios::out);
  fstream file_pos((OUTPUT_filename+"_ProbaSpace_sigpos.dat").c_str(), ios::out);
  fstream file_neg((OUTPUT_filename+"_ProbaSpace_signeg.dat").c_str(), ios::out);

  file << "#N_test = " << N_test << endl;
  file << "#Total number of accessible states = " << NOp_tot << endl;
  file << "#Number of visited states, Nset_test.size() = " << Nset_test.size() << endl;
  file << "#" << endl;
  file << "#1: nb of '1's in the state \t #2: Empirical pba state \t #3: Indep Model pba state \t #4: (x-mu)/sig \t #5: state " << endl;

  for (it = Nset_test.begin(); it!=Nset_test.end(); ++it)
  {
    state = it->first;
    p_data = (it->second) / ((double) N_test);
    p_model = proba_fromIndepModel(best_basis, state);
    sig = sqrt( (p_model*(1.-p_model)) / ((double) N_test) );
    r = (p_data-p_model)/sig;

    file << bitset<n>(state).count() << "\t" << p_data  << "\t" << p_model; // << endl;
    file << "\t" << r;
    file << "  \t " << state << " \t " << bitset<n>(state) << endl;
    if(r > 3) 
    {  
      file_pos << bitset<n>(state).count() << "\t" << p_data  << "\t" << p_model; // << endl;
      file_pos << "\t" << r;
      file_pos << "  \t " << state << " \t " << bitset<n>(state) << endl;
    }
    if(r < -3) 
    {  
      file_neg << bitset<n>(state).count() << "\t" << p_data  << "\t" << p_model; // << endl;
      file_neg << "\t" << r;
      file_neg << "  \t " << state << " \t " << bitset<n>(state) << endl;
    }
  }

  file.close();
  file_pos.close();
  file_neg.close();
}
/******************************************************************************/
/************    Independent model in the PROBABILITY SPACE:    ***************/
/******************************************************************************/
double* Probability_AllStates_Indep(vector<pair<Operator, bool> > best_basis)
{
  double* all_P = (double*)malloc((NOp_tot+1)*sizeof(double));

  for (uint32_t state = 0; state <= NOp_tot; state++)
  {
    all_P[state] = proba_fromIndepModel(best_basis, state);
  }
  return all_P;
}

/******************************************************************************/
/************    Independent model in the FOURIER SPACE:    *******************/
/******************************************************************************/
set<Operator> Rank_m1_Indepmodel(set<Operator> allOp, double *P_indep)
{
  set<Operator>::iterator it;  // Set of all the operators ordered by bias from the data

  Operator Op;
  set<Operator> allOp_buffer;        // New Set of all the operators re-ordered by bias from the data towards the model

  double DKL_diff=0;

  for (it = allOp.begin(); it != allOp.end(); it++)
  {
    Op = (*it);
    Op.p1_M = 0.;
    for (uint32_t state = 0; state <= NOp_tot; state++)
    {
          if( (bitset<n>(Op.bin & state).count())%2 ){ Op.p1_M += P_indep[state]; }    // convention: {-1;1} <--> {1, 0}
    }
    DKL_diff = (Op.p1_D==0 || Op.p1_D==1)? 0: ( (Op.p1_D*log(Op.p1_M) + (1.-Op.p1_D)*log(1.-Op.p1_M)) );
    Op.DKL = - Op.S - DKL_diff;

    allOp_buffer.insert(Op);
  }

  return allOp_buffer;
}


