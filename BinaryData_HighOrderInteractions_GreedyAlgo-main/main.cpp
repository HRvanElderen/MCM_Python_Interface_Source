//*** g++ -I/path/to/eigen -I/path/to/lbfgspp/include -O2 example.cpp ***/
//
//g++ -std=c++11 -I /usr/local/include/eigen3/ -I ./LBFGSpp-master/include -O2 main.cpp ReadDataFile.cpp IndepModel.cpp Models.cpp ModelStatistics.cpp BoltzmannLearning.cpp HeuristicAlgo.cpp BestBasis.cpp
//time ./a.out
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <map>
#include <random>  // for Mersen Twister  // also include "<cmath>" for function "log" (for instance)
#include <set> 
#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono
#include <bitset>

using namespace std;


/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
map<uint32_t, unsigned int> read_datafile(string datafilename, unsigned int *N);   // O(N)  //N = data set size  //READ DATA and STORE them in Nset
void read_Nset (map<uint32_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename); // PRINT Nset in a file  

/******************************************************************************/
/**************************   MODEL TYPE   *********************************/
/******************************************************************************/
double Ising(Interaction I, uint32_t state, int *Op_s);   // s_i in {-1;1}   !! Convention !! {0,1} in data <-> {1,-1} in model
double ACE(Interaction I, uint32_t state, int *Op_s);    // s_i in {0;1}    !! Convention !! {0,1} in data <-> {0,1} in model; But mapping {0,1} <-> {-1,1} in Ising

/******************************************************************************/
/************************** INDEPENDENT MODEL *********************************/
/******************************************************************************/
list<Interaction> IndepModel_fromBasis(vector<pair<Operator, bool> > best_basis);

/************************** Probability space *********************************/
void IndepModel_VS_data_ProbaSpace(vector<pair<Operator, bool> > best_basis, map<uint32_t, unsigned int> Nset_test, unsigned int N_test, string outputfilename);
double* Probability_AllStates_Indep(vector<pair<Operator, bool> > best_basis);

/************************** Fourier space *************************************/
set<Operator> Rank_m1_Indepmodel(set<Operator> allOp, double *P_indep);
set<Operator> Rank_m1_Model(set<Operator> allOp, double *P, double Z);

/******************************** Noise Model *********************************/
list<Interaction> Noise_Model(vector< pair<uint32_t,double> > IndepModel);

/******************************************************************************/
/***************************  ADD NEXT INTERACTION  ***************************/
/******************************************************************************/
Interaction Next_Model(set<Operator> &m1_ranked, list<Interaction> &list_I, map<uint32_t, unsigned int> Nset, unsigned int N, double *L);

/******************************************************************************/
/*************************  MODEL: LIST of INTERACTIONS ***********************/
/******************************************************************************/

/**************************** OTHER MODELS ************************************/
list<Interaction> FullyConnectedPairwise();
list<Interaction> Complete_Model();

void all_int_kbits(int k);

/******************************* ANY MODEL ************************************/
list<Interaction> Build_Model(uint32_t Op_M[], int K);

/************************* PRINT MODEL in TERM ********************************/
void PTerm_ListInteraction (list<Interaction> list_I);

/************************* PRINT MATRIX PAIRWISE INTERACTION ********************************/
void Print_matrix_J_pairwise(string F_name, list<Interaction> list_I);

/******************************************************************************/
/************* PARTITION FUNCTION and PROBA of ALL STATES *********************/
/******************************************************************************/
// return the probability (not normalised) of all the states and compute the partition function
double* Probability_AllStates(double MConvention(Interaction, uint32_t, int*), list<Interaction> list_I, double *Z);

/************************** PRINT in FILE PROBA *******************************/
// model probability of all the states:
void print_proba_state(double* P, double Z, string OUTPUTfilename);
// model probability VS data frequencies of the observed states:
void Model_VS_data_ProbaSpace(double *P, double Z, map<uint32_t, unsigned int> Nset, unsigned int N_test, string outputfilename);

/************************ PRINT in FILE FOURIER SPACE *************************/
void print_m1_Data_VS_Model(map<uint32_t, unsigned int> Nset, double* Op_m1_Model, unsigned int N, string outputfilename);

/******************************************************************************/
/****************************    LOG-LIKELIHOOD    ****************************/
/******************************************************************************/
double LogLikelihood(double *P, double Z, map<uint32_t, unsigned int> Nset);
double LogL_bis(list<Interaction> list_I, double Z, unsigned int N);
double LogL(double MConvention(Interaction, uint32_t, int*), list<Interaction> list_I, unsigned int N);

double LogLikelihood_CompleteModel(map<uint32_t, unsigned int> Nset, unsigned int N);
/******************************************************************************/
/**************************** Operator Averages *******************************/
/******************************************************************************/
void Model_averages(double MConvention(Interaction, uint32_t, int*), double *P, double Z, list<Interaction> &list_I, unsigned int N);

void empirical_averages_Ising(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);
void empirical_averages_ACE(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);

double* AllOp_av(double *P, double Z, unsigned int N);  // _all_Model_averages in ISING convention
double* AllOp_m1(double *P, double Z, unsigned int N);

void Print_Correlations(map<uint32_t, unsigned int> Nset, unsigned int N);

/******************************************************************************/
/******************************   Ranked m1   *********************************/
/******************************************************************************/
set<Operator> AllOp_m1_data_rank(map<uint32_t, unsigned int> Nset, unsigned int N);

void Fill_m1_model(set<Operator> &allOp, double *P, double Z);

//set<Operator> AllOp_m1_Model_ranked(double *P, double Z, unsigned int N);

//void print_m1_Data_VS_Model_Mranked(map<uint32_t, unsigned int> Nset, set<Operator> Op_m1_Model, unsigned int N, string outputfilename);
void print_m1_Data_VS_Model_DKLranked(set<Operator> Op_m1_Model, unsigned int N, string outputfilename);

/******************************************************************************/
/******************************   Pairwise model   ****************************/
/******************************************************************************/
set<Operator> PairOp_m1_data_rank(map<uint32_t, unsigned int> Nset, unsigned int N);

/******************************************************************************/
/****************************  Best Basis  ************************************/
/******************************************************************************/
vector<pair<Operator, bool> > select_best_basis(set<Operator> &allOp_data, double Nd, double *LogL);
void printfile_BestBasis(vector<pair<Operator, bool> > Op, double Nd, string name);

/******************************************************************************/
/**************************** Boltzmann Learning ******************************/
/******************************************************************************/
double BoltzmannLearning_Ising(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);
double BoltzmannLearning_ACE(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);

void gIsing_from_gACE(list<Interaction> &list_I);

/******************************************************************************/
/****************************  Analysed Data  *********************************/
/******************************************************************************/
void USSC_BestModel_all_K45(map<uint32_t, unsigned int>  Nset, unsigned int N);

/******************************************************************************/
/****************************  BEST INDEP MODEL *******************************/
/******************************************************************************/
set<Operator> best_Indep_Model(map<uint32_t, unsigned int> Nset, unsigned int N, list<Interaction> &list_I, double *L, string OUTPUT_filename)
{
  double Nd = (double) N;

  cout << endl << "--->> Search Best Basis.. \t";
  set<Operator> m1_ranked = AllOp_m1_data_rank(Nset, N);
  cout << "Total number of Operators = " << m1_ranked.size() << endl << endl;
  print_m1_Data_VS_Model_DKLranked(m1_ranked, N, OUTPUT_filename + "_FourierSpace_dataVSNoise_DKLranked.dat");

/// This function is updating the layer values: Rem, still an issue with the layers values, 
/// as we remove the basis op from the "m1_Indep_ranked" while searching for the best basis:
  vector<pair<Operator, bool> > best_basis = select_best_basis(m1_ranked, ((double) N), L); 
  printfile_BestBasis(best_basis, (double) N, "");

  //Print BIC:
  fstream BIC_file( (OUTPUT_filename + "_BestIndep_toK"+ to_string(n) +".dat").c_str(), ios::out);
  BIC_file << "##1:K \t 2:maxlogL \t 3: BIC \t 4:worst_MDL \t 5:Op" << endl << "##" << endl;

  int K=0; double L_buff = -Nd * n * log(2.);
  for (vector<pair<Operator, bool> >::iterator Op = best_basis.begin(); Op != best_basis.end(); Op++)
  {
    K++;       BIC_file << K << "\t ";
    L_buff += ( - Nd * (*Op).first.S + Nd * log(2.) );      
    BIC_file << L_buff << "\t ";
    BIC_file << L_buff - K * log( Nd /(2.*M_PI)) / 2. << "\t";   //BIC
    BIC_file << L_buff - K * log( (Nd*M_PI) /2.) / 2. << "\t";   //MDL worst
    BIC_file << (*Op).first.bin << "\t" << bitset<n>((*Op).first.bin) << endl;
  }
  BIC_file.close();

  list_I = IndepModel_fromBasis(best_basis);
  PTerm_ListInteraction (list_I);

  K = list_I.size(); 
  cout << "Number of parameters, K =  " << K << endl;
  cout << "Max Log-Likelihood, L = " << (*L) << endl;
  cout << "BIC: " << (*L) - K * log( (Nd) /(2.*M_PI)) / 2. <<  endl << endl;  

  // Print all probabilities:
  double* P_indep = Probability_AllStates_Indep(best_basis); // Proba each state (normalised)
  IndepModel_VS_data_ProbaSpace(best_basis, Nset, N, OUTPUT_filename);

  // Indep model all m1:
  cout << "Model Rank all DKL.." << endl;
  set<Operator> allOp_buffer = Rank_m1_Indepmodel(m1_ranked, P_indep);
  m1_ranked.clear();
  //m1_ranked.swap(allOp_buffer);

  print_m1_Data_VS_Model_DKLranked(m1_ranked, N, OUTPUT_filename + "_FourierSpace_dataVSIndep_DKLranked.dat");

  return allOp_buffer;
}

set<Operator> best_Indep_PairwiseModel(map<uint32_t, unsigned int> Nset, unsigned int N, list<Interaction> &list_I, double *L, string OUTPUT_filename)
{
  double Nd = (double) N;

  cout << endl << "--->> Search Best Basis.. \t";
  set<Operator> m1_ranked = PairOp_m1_data_rank(Nset, N);
  cout << "Total number of Operators = " << m1_ranked.size() << endl << endl;
  print_m1_Data_VS_Model_DKLranked(m1_ranked, N, OUTPUT_filename + "_Pair_FourierSpace_dataVSNoise_DKLranked.dat");

/// This function is updating the layer values: Rem, still an issue with the layers values, 
/// as we remove the basis op from the "m1_Indep_ranked" while searching for the best basis:
  vector<pair<Operator, bool> > best_basis = select_best_basis(m1_ranked, ((double) N), L); 
  printfile_BestBasis(best_basis, Nd, "Pair_");

  //Print BIC:
  fstream BIC_file( (OUTPUT_filename + "_Pair_BestIndep_toK"+ to_string(n) +".dat").c_str(), ios::out);
  BIC_file << "##1:K \t 2:maxlogL \t 3: BIC \t 4:worst_MDL \t 5:Op" << endl << "##" << endl;

  int K=0; double L_buff = -Nd * n * log(2.);
  for (vector<pair<Operator, bool> >::iterator Op = best_basis.begin(); Op != best_basis.end(); Op++)
  {
    K++;       BIC_file << K << "\t ";
    L_buff += ( - Nd * (*Op).first.S + Nd * log(2.) );      
    BIC_file << L_buff << "\t ";
    BIC_file << L_buff - K * log( Nd /(2.*M_PI)) / 2. << "\t";   //BIC
    BIC_file << L_buff - K * log( (Nd*M_PI) /2.) / 2. << "\t";   //MDL worst
    BIC_file << (*Op).first.bin << "\t" << bitset<n>((*Op).first.bin) << endl;
  }
  BIC_file.close();

  list_I = IndepModel_fromBasis(best_basis);
  PTerm_ListInteraction (list_I);

  K = list_I.size(); 
  cout << "Number of parameters, K =  " << K << endl;
  cout << "Max Log-Likelihood, L = " << (*L) << endl;
  cout << "BIC: " << (*L) - K * log( ((double) N) /(2.*M_PI)) / 2. <<  endl << endl;  

  // Print all probabilities:
  double* P_indep = Probability_AllStates_Indep(best_basis); // Proba each state (normalised)
  IndepModel_VS_data_ProbaSpace(best_basis, Nset, N, OUTPUT_filename + "_Pair_");

  // Indep model all m1:
  cout << "Model Rank all DKL.." << endl;
  set<Operator> allOp_buffer = Rank_m1_Indepmodel(m1_ranked, P_indep);
  m1_ranked.clear();
  //m1_ranked.swap(allOp_buffer);

  print_m1_Data_VS_Model_DKLranked(m1_ranked, N, OUTPUT_filename + "_Pair_FourierSpace_dataVSIndep_DKLranked.dat");

  return allOp_buffer;
}


/******************************************************************************/
/*************************  HEURISTIC ALGORITHM *******************************/
/******************************************************************************/
void Heuristic_Model_K(map<uint32_t, unsigned int> Nset, unsigned int N, list<Interaction> &list_I, set<Operator> &m1_ranked, string OUTPUT_filename, int Kmax = n)
{
  cout << endl << "--->> Search next most relevant Operators.. \t";
  cout << "Number of Operators left = " << m1_ranked.size() << endl << endl;

  //variables:
  int K = list_I.size();
  set<Operator> allOp_buffer;
  double L=0;  // max-logLikelihood
  double Z=0;  // partition function
  double *P;
  double *m1;
  double Nd = (double) N;

  //Print BIC:
  fstream BIC_file( (OUTPUT_filename + "_Best_toK"+ to_string(Kmax) +".dat").c_str(), ios::out);
  BIC_file << "##1:K \t 2:maxlogL \t 3: BIC \t 4:worst_MDL \t 5:Op" << endl << "##" << endl;

  while(K < Kmax)
  {
    K += 1;
    Interaction I_K = Next_Model(m1_ranked, list_I, Nset, N, &L);
    PTerm_ListInteraction (list_I);

    BIC_file << K << "\t ";
    BIC_file << L << "\t ";
    BIC_file << L - K * log( Nd /(2.*M_PI)) / 2. << "\t";
    BIC_file << L - K * log( (Nd*M_PI) /2.) / 2. << "\t";
    BIC_file << I_K.Op << "\t" << bitset<n>(I_K.Op) << endl;

    // UPDATE p1 and DKL:
    P = Probability_AllStates(Ising, list_I, &Z); // Proba each state, non normalised    
    //P = Probability_AllStates(ACE, list_I, &Z); // Proba each state, non normalised // ACE version
    allOp_buffer = Rank_m1_Model(m1_ranked, P, Z);
    m1_ranked.clear();
    m1_ranked.swap(allOp_buffer);

    //PRINT OUT:
    //cout << "--->> Print file: Best model, Fourier Space, m1 of all the states.." << endl;
    //print_proba_state(P, Z, OUTPUT_filename + "_AllProba_K"+to_string(K)+".dat");
    ////Model_VS_data_ProbaSpace(P, Z, Nset, N, OUTPUT_filename + "_ProbaSpace_K"+to_string(K)+".dat");// model probability VS data frequencies of the observed states
    print_m1_Data_VS_Model_DKLranked(m1_ranked, N, OUTPUT_filename + "_FourierSpace_dataVSK" + to_string(K) + "_DKLranked.dat");
    //m1 = AllOp_m1(P, Z, N);

    //print_m1_Data_VS_Model(Nset, m1, N, OUTPUT_filename + "_FourierSpace_dataVSK"+to_string(K)+".dat");
  
    cout << "Number of Operators left = " << m1_ranked.size() << endl << endl;
  }
  BIC_file.close();
}

/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/
//Rem : ---> 2^30 ~ 10^9
int main(int argc, char **argv)
{
//*******************************************
//********** READ DATA FILE:  ***************     -->  data are put in Nset:
//*******************************************    // Nset[mu] = # of time state mu appears in the data set

  cout << "--->> Create OUTPUT Folder: (if needed) ";
  system( ("mkdir " + OUTPUT_directory).c_str() );
  cout << endl;

  string filename;
  if (argc > 1) {    filename = argv[1];  }  // output file name
  else {    cout << "The execution is missing the input filename as an argument." << endl;  }

  string INPUT_filename = INPUT_folder + filename + ".dat";
  cout << "Read the input file: " << INPUT_filename << endl << endl;

  string OUTPUT_filename = OUTPUT_directory + filename;
  //OUTPUT_filename_Data_VS_Model = OUTPUT_directory + filename + "_dataVSmodel.dat"; 

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************  TRAINING = TEST on the WHOLE dataset:  **********************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "************************************  Read dataset:  **************************************" << endl;

  unsigned int N = 0;  // set in 'read_datafile()'
  map<uint32_t, unsigned int> Nset = read_datafile(INPUT_filename, &N);  //  double Nd = (double) N;
  read_Nset(Nset, N, OUTPUT_filename + "_Nset.dat");   // print Nset in a file
  double Nd = (double) N;

//  cout << "--->> Print file: Correlations in Ising C1 convention.." << endl;
//  Print_Correlations(Nset, N);

  int K=0; double L=0;
  
  cout << endl << "***********************************  Complete model:  *************************************" << endl;

  K = NOp_tot;
  L = LogLikelihood_CompleteModel(Nset, N);
  
  cout << "Number of parameters, K =  " << K << endl;
  cout << "Max Log-Likelihood, L = " << L << endl;
  cout << "BIC: " << L - K * log( Nd /(2.*M_PI)) / 2. <<  endl << endl;  

  list<Interaction> list_I;

//  cout << endl << "*****************************  Best Independent model:  **********************************" << endl;
  set<Operator> m1_ranked = best_Indep_Model(Nset, N, list_I, &L, OUTPUT_filename);

  cout << endl << "*****************************  Best model with K interactions:  **********************************" << endl;
  int Kmax = 20;
  //set<Operator> m1_ranked = AllOp_m1_data_rank(Nset, N);
  Heuristic_Model_K(Nset, N, list_I, m1_ranked, OUTPUT_filename, Kmax);

//  cout << endl << "*****************************  Best independent Model, among Pairwise models:  **********************************" << endl;
//  set<Operator> m1_ranked = best_Indep_PairwiseModel(Nset, N, list_I, &L, OUTPUT_filename);

//  cout << endl << "*****************************  Best pairwise model with K interactions:  **********************************" << endl;
//  int Kmax = 40;
//  set<Operator> m1_ranked = PairOp_m1_data_rank(Nset, N);
//  Heuristic_Model_K(Nset, N, list_I, m1_ranked, OUTPUT_filename, Kmax);


  return 0;
}









