//#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <list>
#include <bitset> 
#include <set>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

unsigned int k1_op(map<uint32_t, unsigned int> Nset, uint32_t Op); //k1 += (bitset<n>( ((*it).first) & op ).count() % 2)*((*it).second); 

/******************************************************************************/
/**************************** Boltzmann Learning ******************************/
/******************************************************************************/
double BoltzmannLearning_Ising(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);
double BoltzmannLearning_ACE(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);

/******************************************************************************/
/**************************    ADD NEXT OPERATOR:    **************************/
/******************************************************************************/
Interaction Next_Model(set<Operator> &m1_ranked, list<Interaction> &list_I, map<uint32_t, unsigned int> Nset, unsigned int N, double *L)
{
  // select the new relevant interaction:
  Interaction I;
  Operator Op_new = *(m1_ranked.begin());   // take the first element
  m1_ranked.erase(m1_ranked.begin());  // remove the 1rst element

  I.Op = Op_new.bin;
  I.k = bitset<n>(I.Op).count();
  I.g = 0;  I.g_Ising = 0;  I.g_ACE = 0;
  I.av_D = 0; I.av_M = 0;
  list_I.push_back(I);

  // fit the parameters:
  double K = list_I.size();
  cout << "Op[K" << K << "] = " << Op_new.bin << " = " << bitset<n>(Op_new.bin) << endl;
  cout << "--->> Finding best parameters (gradient descent Ising).." << endl;
  cout << "Best Model: K = " << K << endl;
  *L = BoltzmannLearning_Ising (Nset, list_I, N);
  cout << "Max Log-Likelihood: L = " << (*L) << endl;
  cout << "BIC: " << (*L) - K * log( ((double) N) /(2.*M_PI)) / 2. <<  endl << endl;

  return I;
}

/******************************************************************************/
/*********************************  Averages  *********************************/
/******************************************************************************/
set<Operator> Rank_m1_Model(set<Operator> allOp, double *P, double Z)
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
          if( (bitset<n>(Op.bin & state).count())%2 ){ Op.p1_M += P[state]; }    // convention: {-1;1} <--> {1, 0}
    }
    Op.p1_M = Op.p1_M / Z;

    DKL_diff = (Op.p1_D==0 || Op.p1_D==1)? 0: ( (Op.p1_D*log(Op.p1_M) + (1.-Op.p1_D)*log(1.-Op.p1_M)) );
    Op.DKL = - Op.S - DKL_diff;

    allOp_buffer.insert(Op);
  }
  return allOp_buffer;
}

void Fill_m1_model(set<Operator> &allOp, double *P, double Z)  // _all_Model_averages in ISING convention
{
  set<Operator>::iterator it_Op; // Set of all the operators ordered by bias from the data

  for (it_Op = allOp.begin(); it_Op != allOp.end(); it_Op++)//NOp_tot
  {
    (*it_Op).p1_M = 0.;
    for (uint32_t state = 0; state <= NOp_tot; state++)
    {
          if( (bitset<n>((*it_Op).bin & state).count())%2 ){ (*it_Op).p1_M += P[state]; }    // convention: {-1;1} <--> {1, 0}
    }
    (*it_Op).p1_M = (*it_Op).p1_M / Z;
  }
}


/******************************************************************************/
/********************** PRINT FOURIER SPACE ***********************************/
/******************************************************************************/
void print_m1_Data_VS_Model_DKLranked(set<Operator> Op_m1_Model, unsigned int N, string outputfilename)
{
  double Nd = (double) N;
  set<Operator>::iterator Op;

//  cout << "--->> Create outputfile " << outputfilename << " :" << endl;
  fstream fichier(outputfilename.c_str(), ios::out);
  fichier << "## 1:rank 2:p1_data \t 3:p1_Model \t4:S(data)/N \t5:DKL/N \t6:Op \t 7:Op_bin \t 8:order \t 9:layer" << endl;

  double sig_p1 = 0.5/sqrt(Nd);
  double p1 = 0.5 + 3*sig_p1;
  double Lmin_delta = Nd*p1*log(p1) + Nd*(1.-p1)*log(1-p1);
  fichier << "## For non-informative operators, p = 0.5 +/- " << 3*sig_p1 << " (3*sig_p1) ";
  fichier << ",\t Lmin = [" << -Nd*log(2.) << ", " << Lmin_delta << "]"<< endl;
  fichier << "## 99.73\% of the non-informative operators should end within this interval." << endl;
  fichier << "## " << endl;

  int rank = 1;
  for (Op = Op_m1_Model.begin(); Op != Op_m1_Model.end(); Op++)
  {
    fichier << rank; rank++;
    fichier << "\t" << (*Op).p1_D;
    fichier << "\t" << (*Op).p1_M;
    fichier << "\t" << (*Op).S;
    fichier << "\t" << (*Op).DKL;

    fichier << "\t" << (*Op).bin << "\t" << bitset<n>((*Op).bin);
    fichier << "\t" << bitset<n>((*Op).bin).count();
    fichier << "\t" << (*Op).layer;
    fichier << endl;
  }

  fichier.close();
}

/******************************************************************************/
/*************  FOURIER SPACE  --  ALL STATES  and  VS DATA  ******************/
/******************************************************************************/
void print_m1_Data_VS_Model(map<uint32_t, unsigned int> Nset, double* Op_m1_Model, unsigned int N, string outputfilename)
{
  double Nd = (double) N;

//  cout << "--->> Create outputfile " << outputfilename << " :" << endl;
  fstream fichier(outputfilename.c_str(), ios::out);
  fichier << "## 1:p1_data \t 2:p1_Model \t3:Op \t 4:Op_bin \t 5:Op_count" << endl;

  double sig_p1 = 0.5/sqrt(Nd);
  double p1 = 0.5 + 3*sig_p1;
  double Lmin_delta = Nd*p1*log(p1) + Nd*(1.-p1)*log(1-p1);
  fichier << "## For non-informative operators, p = 0.5 +/- " << 3*sig_p1 << " (3*sig_p1) ";
  fichier << ",\t Lmin = [" << -Nd*log(2.) << ", " << Lmin_delta << "]"<< endl;
  fichier << "## 99.73\% of the non-informative operators should end within this interval." << endl;
  fichier << "## " << endl;

  for (uint32_t Op = 1; Op <= NOp_tot; Op++)
  {
    p1 = k1_op(Nset, Op)/Nd;
    fichier << "\t" << p1 ;
    fichier << "\t" << Op_m1_Model[Op] ;

    fichier << "\t" << Op << "\t" << bitset<n>(Op) << "\t" << bitset<n>(Op).count();
    fichier << endl;
  }

  fichier.close();
}
