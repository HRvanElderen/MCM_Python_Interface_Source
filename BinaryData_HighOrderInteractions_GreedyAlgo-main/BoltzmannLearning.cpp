//#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <bitset> 

using namespace std;

//For learning:
//#include <iostream>
#include <Eigen/Core>
#include <LBFGS.h>

using namespace Eigen;  //using Eigen::VectorXd;
using namespace LBFGSpp;


/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"
void PTerm_ListInteraction (list<Interaction> list_I);

/******************************************************************************/
/**************************   MODEL TYPE   *********************************/
/******************************************************************************/
double Ising(Interaction I, uint32_t state, int *Op_s);   // s_i in {-1;1}   !! Convention !! {0,1} in data <-> {1,-1} in model
double ACE(Interaction I, uint32_t state, int *Op_s);    // s_i in {0;1}    !! Convention !! {0,1} in data <-> {0,1} in model; But mapping {0,1} <-> {-1,1} in Ising

/******************************************************************************/
/************* PARTITION FUNCTION and PROBA of ALL STATES *********************/
/******************************************************************************/
// return the probability (not normalised) of all the states and compute the partition function
double* Probability_AllStates(double MConvention(Interaction, uint32_t, int*), list<Interaction> list_I, double *Z);

/******************************************************************************/
/****************************   MODEL STATISTICS   ****************************/
/******************************************************************************/
void Model_averages(double MConvention(Interaction, uint32_t, int*), double *P, double Z, list<Interaction> &list_I, unsigned int N);

/******************************************************************************/
/**************************     EMPIRICAL AVERAGES    *************************/
/******************************************************************************/

/************************     ISING MODEL {+1; -1}    *************************/
unsigned int k1_op(map<uint32_t, unsigned int> Nset, uint32_t op); //k1 += (bitset<n>( ((*it).first) & op ).count() % 2)*((*it).second); 

double op_av_Ising(map<uint32_t, unsigned int> Nset, uint32_t op, unsigned int N)
{
  return ( ((double) N) - 2.*k1_op(Nset, op) ) / ((double) N); // [ [-1 * k1] + [+1 * (N-k1)] ] / N
}

/************************     ACE MODEL {0; 1}    *************************/
unsigned int k1_op_ACE(map<uint32_t, unsigned int> Nset, uint32_t op)  // Complexity = O(|Nset|)
{
  unsigned int k1=0;
  map<uint32_t, unsigned int>::iterator it;  // iterator on Nset

  for (it = Nset.begin(); it!=Nset.end(); ++it)
    {    k1 += ( ( ((*it).first) & op ) == op ) *((*it).second);   }

  return k1;
}

double op_av_ACE(map<uint32_t, unsigned int> Nset, uint32_t op, unsigned int N)
{
  return ( k1_op_ACE(Nset, op) ) / ((double) N); // [ [1 * k1] + [0 * (N-k1)] ] / N
}

double op_av_ACE_0(map<uint32_t, unsigned int> Nset, uint32_t op, unsigned int N)
{
  unsigned int k0=0;
  map<uint32_t, unsigned int>::iterator it;  // iterator on Nset

  for (it = Nset.begin(); it!=Nset.end(); ++it)
    {    k0 += ( ( ((*it).first) & op ) == 0 ) *((*it).second);   }

  return ((double) k0 ) / ((double) N); // [ [1 * k1] + [0 * (N-k1)] ] / N
}

/************************    Empirical averages all op    *************************/
void empirical_averages_Ising(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N) 
{
  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    (*I).av_D = op_av_Ising(Nset, (*I).Op, N);
  }
}

void empirical_averages_ACE(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N) 
{
  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    (*I).av_D = op_av_ACE(Nset, (*I).Op, N);
  }
}

/******************************************************************************/
/********************     EMPIRICAL CORRELATION FUNCTION    *******************/
/******************************************************************************/
void Print_Correlations(map<uint32_t, unsigned int> Nset, unsigned int N)
{
  double **Corr_Ising = (double **)malloc(n*sizeof(double*));  //Correlations <si sj>-<si><sj>
  for (int i1 = 0; i1 < n; i1++) {  Corr_Ising[i1] = (double*)malloc(n*sizeof(double));  }

  double **Cij =  (double **)malloc(n*sizeof(double*));   //Cij = P[si=1 and sj=1]
  for (int i1 = 0; i1 < n; i1++) {  Cij[i1] = (double*)malloc(n*sizeof(double));  }

  double **Cij_fifj =  (double **)malloc(n*sizeof(double*));   //Cij - fi fj = P[si=1 and sj=1] - P[si=1] * P[sj=1] 
  for (int i1 = 0; i1 < n; i1++) {  Cij_fifj[i1] = (double*)malloc(n*sizeof(double));  }

  double **Cij_0 =  (double **)malloc(n*sizeof(double*));   //Cij_0 = P[si=0 and sj=0]
  for (int i1 = 0; i1 < n; i1++) {  Cij_0[i1] = (double*)malloc(n*sizeof(double));  }

  double **Cij_fifj_0 =  (double **)malloc(n*sizeof(double*));   //= P[si=0 and sj=0] - P[si=0] * P[sj=0] 
  for (int i1 = 0; i1 < n; i1++) {  Cij_fifj_0[i1] = (double*)malloc(n*sizeof(double));  }

  unsigned int i1, i2;
  uint32_t Op1, Op2;
  double Op1_av, Op2_av, Op3_av;

  for (i1 = 0; i1 < n; i1++)
  {
    Op1 = 1 << (n-1-i1);
    Op1_av = op_av_Ising(Nset, Op1, N);

    // All the diagonals are set to zero...
    Corr_Ising[i1][i1] = 0; //1-Op1_av*Op1_av;
    Cij[i1][i1] = 0;
    Cij_fifj[i1][i1] = 0;
    Cij_0[i1][i1] = 0;
    Cij_fifj_0[i1][i1] = 0;

    i2 = i1 + 1;
    while(i2 < n)
    {
      Op2 = 1 << (n-1-i2);
      Op2_av = op_av_Ising(Nset, Op2, N);
      Op3_av = op_av_Ising(Nset, Op1+Op2, N);

      //symmetric matrices:
      Corr_Ising[i1][i2] = Op3_av - Op1_av * Op2_av;      Corr_Ising[i2][i1] = Corr_Ising[i1][i2];

      Cij[i1][i2] = op_av_ACE(Nset, Op1+Op2, N);          Cij[i2][i1] = Cij[i1][i2];
      Cij_fifj[i1][i2] = Cij[i1][i2] - op_av_ACE(Nset, Op1, N)*op_av_ACE(Nset, Op2, N);   Cij_fifj[i2][i1] = Cij_fifj[i1][i2];

      Cij_0[i1][i2] = op_av_ACE_0(Nset, Op1+Op2, N);          Cij_0[i2][i1] = Cij_0[i1][i2];
      Cij_fifj_0[i1][i2] = Cij_0[i1][i2] - op_av_ACE_0(Nset, Op1, N)*op_av_ACE_0(Nset, Op2, N);   Cij_fifj_0[i2][i1] = Cij_fifj_0[i1][i2];

      i2++;
    }
  }

  //Cij = P[si=1 and sj=1]
  fstream fichier_Cij((OUTPUT_directory+"Cij.dat").c_str(), ios::out);
  fstream fichier_Cij_fifj((OUTPUT_directory+"Cij_fifj.dat").c_str(), ios::out);

  fstream fichier_Cij_0((OUTPUT_directory+"Cij_0.dat").c_str(), ios::out);
  fstream fichier_Cij_fifj_0((OUTPUT_directory+"Cij_fifj_0.dat").c_str(), ios::out);

  for (i1 = 0; i1 < n; i1++)
  {
    for(i2 = 0; i2 < n; i2++)
    {
      fichier_Cij << Cij[i1][i2] << " ";
      fichier_Cij_fifj << Cij_fifj[i1][i2] << " ";

      fichier_Cij_0 << Cij_0[i1][i2] << " ";
      fichier_Cij_fifj_0 << Cij_fifj_0[i1][i2] << " ";
    }
    fichier_Cij << endl;
    fichier_Cij_fifj << endl;

    fichier_Cij_0 << endl;
    fichier_Cij_fifj_0 << endl;
  }
  fichier_Cij.close();
  fichier_Cij_fifj.close();
  fichier_Cij_0.close();
  fichier_Cij_fifj_0.close();

  //Correlations :
  fstream fichier((OUTPUT_directory+"Correlations_IsingC1.dat").c_str(), ios::out);
  fstream fichier_log((OUTPUT_directory+"CorrelationsLog_IsingC1.dat").c_str(), ios::out);
  for (i1 = 0; i1 < n; i1++)
  {
    for(i2 = 0; i2 < n; i2++)
    {
      fichier << Corr_Ising[i1][i2] << " ";
      if (i1 == i2) { fichier_log << 0 << " ";  }
      else {  fichier_log << log(abs(Corr_Ising[i1][i2])) << " "; }
      cout << Corr_Ising[i1][i2] << " " << log(abs(Corr_Ising[i1][i2]))<< endl;
    }
    fichier << endl;
    fichier_log << endl;
  }
  fichier.close();
  fichier_log.close();
}


/******************************************************************************/
/****************************  LogLikelihood  *********************************/
/******************************************************************************/
double LogL(double MConvention(Interaction, uint32_t, int*), list<Interaction> list_I, unsigned int N)
{
  double Z = 0; Probability_AllStates(MConvention, list_I, &Z);
  double LogLi = 0;
//cout << "Zbis = " << Z << endl;
  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    LogLi += (*I).g * (*I).av_D;
  } 
  LogLi -= log(Z);
  //cout << "LogLbis = " << ((double) N)*LogLi << endl;

  return N*LogLi;
}

double LogL_bis(list<Interaction> list_I, double Z, unsigned int N)
{
  double LogLi = 0;
  list<Interaction>::iterator I;

  for (I = list_I.begin(); I != list_I.end(); I++)
  {
    LogLi += (*I).g * (*I).av_D;
  } 
  LogLi -= log(Z);

  return ((double) N)*LogLi;
}
/******************************************************************************/
/****************************  Boltzmann Learning  ****************************/
/******************************************************************************/

/************************     ISING MODEL {+1; -1}    *************************/
class LogL_class_Ising
{
private:
  list<Interaction> li_I;
  unsigned int N;
public:
 LogL_class_Ising(list<Interaction> list_I, unsigned int N_) : li_I(list_I), N(N_) {}
 double operator()(const VectorXd& x, VectorXd& grad)
 {
  list<Interaction>::iterator I;

  int i=0;
  for (I = li_I.begin(); I != li_I.end(); I++)
  {    (*I).g = x[i];   i++;  }

  double Z = 0; double *P = Probability_AllStates(Ising, li_I, &Z);  //cout << "Z = " << Z << "\t ";
  Model_averages(Ising, P, Z, li_I, N);
  //PTerm_ListInteraction (li_I);

  double LogLi = LogL_bis(li_I, Z, N);      //cout << "LogLi_1 = " << LogLi << endl;
  //LogLi = LogL(li_I, N);   cout << "LogLi_2 = " << LogLi << endl;
  //double LogLi = 0;

//cout << "grad: \t" << endl;
  i=0;
  for (I = li_I.begin(); I != li_I.end(); I++)
  {
    //LogLi += x[i] * (*I).av_D;
    //cout << i << "\t g[i] = " << x[i] << "\t <Op_i> = " << (*I).av_D << "\t LogLi = " << LogLi << endl;
    grad[i] = -((double) N) * ( (*I).av_D - (*I).av_M); 
    //cout << " " << N << " " << (*I).av_D << " " << (*I).av_M << " " << -N * ( (*I).av_D - (*I).av_M) << " " << grad[i] << endl;
    i++;
  } 
  //LogLi -= log(Z);  cout << "LogLi = " << N*LogLi << endl;
  //cout << endl;
  free(P);
  return -LogLi;
 }
};

double BoltzmannLearning_Ising(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N)
{
  unsigned int K = list_I.size();
  empirical_averages_Ising(Nset, list_I, N);

  // Set up parameters:
  LBFGSParam<double> param;
  param.epsilon = 1e-6;  //1e-6;
  param.max_iterations = 100;  //100
 
  // Create solver and function object:
  LBFGSSolver<double> solver(param);
  LogL_class_Ising LogLi(list_I, N);

  // Initial guess:
  VectorXd g = VectorXd::Zero(K);

  // g will be overwritten to be the best point found:
  double maxLogLi;
  int niter = solver.minimize(LogLi, g, maxLogLi);

  //updating the parameters:
  int i=0;  list<Interaction>::iterator I;
  for (I = list_I.begin(); I != list_I.end(); I++)
  {    
    (*I).g_Ising = g[i];   (*I).g = g[i];    i++;
  }

  //updating the model averages:
  double Z = 0; double *P = Probability_AllStates(Ising, list_I, &Z);
  //cout << "Z = " << Z << "\t ";
  Model_averages(Ising, P, Z, list_I, N);

/*
  double Z = 0; double *P = Probability_AllStates(model_convention, list_I, &Z);
  cout << "Z = " << Z << "\t ";
  Model_averages(Ising, P, Z, list_I, N);
  PTerm_ListInteraction (list_I); cout << endl;
  cout << "N = " << N << endl;

  cout << niter << " iterations" << std::endl;
  cout << "g = \n" << g.transpose() << std::endl;*/
  //cout << "maxLogLi = " << -maxLogLi << std::endl;
  free(P);
  return -maxLogLi;
}

/************************     ACE MODEL {0; 1}    *************************/
class LogL_class_ACE
{
private:
  list<Interaction> li_I;
  unsigned int N;
public:
 LogL_class_ACE(list<Interaction> list_I, unsigned int N_) : li_I(list_I), N(N_) {}
 double operator()(const VectorXd& x, VectorXd& grad)
 {
  list<Interaction>::iterator I;

  int i=0;
  for (I = li_I.begin(); I != li_I.end(); I++)
  {    
    (*I).g = x[i];   i++;
  }

  double Z = 0; double *P = Probability_AllStates(ACE, li_I, &Z);
  //cout << "Z = " << Z << "\t ";
  Model_averages(ACE, P, Z, li_I, N);
  //PTerm_ListInteraction (li_I);

  double LogLi = LogL_bis(li_I, Z, N);      //cout << "LogLi_1 = " << LogLi << endl;
  //LogLi = LogL(li_I, N);   cout << "LogLi_2 = " << LogLi << endl;
  //double LogLi = 0;

//cout << "grad: \t" << endl;
  i=0;
  for (I = li_I.begin(); I != li_I.end(); I++)
  {
    //LogLi += x[i] * (*I).av_D;
    //cout << i << "\t g[i] = " << x[i] << "\t <Op_i> = " << (*I).av_D << "\t LogLi = " << LogLi << endl;
    grad[i] = -((double) N) * ( (*I).av_D - (*I).av_M);
    //cout << " " << N << " " << (*I).av_D << " " << (*I).av_M << " " << grad[i] << endl;
    i++;
  } 
  //LogLi -= log(Z);  cout << "LogLi = " << N*LogLi << endl;
  free(P);
  return -LogLi;
 }
};

double BoltzmannLearning_ACE(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N)
{
  unsigned int K = list_I.size();

  empirical_averages_ACE(Nset, list_I, N); 

  // Set up parameters
  LBFGSParam<double> param;
  param.epsilon = 1e-6; // 1e-6
  param.max_iterations = 100; //100
 
  // Create solver and function object
  LBFGSSolver<double> solver(param);
  LogL_class_ACE LogLi(list_I, N);

  // Initial guess
  VectorXd g = VectorXd::Zero(K);

  // x will be overwritten to be the best point found
  double maxLogLi;
  int niter = solver.minimize(LogLi, g, maxLogLi);

  //updating the parameters:
  int i=0;  list<Interaction>::iterator I;
  for (I = list_I.begin(); I != list_I.end(); I++)
  {    
    (*I).g_ACE = g[i];   (*I).g = g[i]; i++;
  }

  //updating the model averages:
  double Z = 0; double *P = Probability_AllStates(ACE, list_I, &Z);
  //cout << "Z = " << Z << "\t ";
  Model_averages(ACE, P, Z, list_I, N);

/*
  PTerm_ListInteraction (list_I); cout << endl;
*/
  //cout << niter << " iterations" << std::endl;
  //cout << "g = \n" << g.transpose() << std::endl;
  free(P);
  return -maxLogLi;
}

/******************************************************************************/
/****************************     INTERACTIONS     ****************************/
/*********************     {0,1} OR {-1,1} CONVENTION     *********************/
/******************************************************************************/
void gIsing_from_gACE(list<Interaction> &li_I)
{
  map<uint32_t, double> h_add;
  list<Interaction>::iterator I;

  uint32_t Op_i = 1, Op_j = 1;

//J_ij
  for (I = li_I.begin(); I != li_I.end(); I++)
  {    
    if((*I).k == 2)
      {   (*I).g_ACE = 4*((*I).g);  h_add[((*I).Op)] = (*I).g; }
  }

//h_i
  for (I = li_I.begin(); I != li_I.end(); I++)
  {    
    if((*I).k == 1)
    {   (*I).g_ACE = 2*((*I).g);  
      Op_i = (*I).Op;
      Op_j = 1;
      for (int i=0; i<n; i++) 
      { 
        if (Op_j != Op_i)
          { (*I).g_ACE -= 2 * h_add[Op_i + Op_j]; }
        Op_j = Op_j << 1;
      }
    }
  }
}



