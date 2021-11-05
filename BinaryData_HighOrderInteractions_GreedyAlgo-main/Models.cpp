//#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <bitset> 

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/**************************   PRINT MODEL IN TERMINAL   ***********************/
/*****************************    ALL INTERACTIONS    *************************/
/******************************************************************************/
void PTerm_ListInteraction (list<Interaction> list_I)
{
  list<Interaction>::iterator it;

  cout << "# 1:Op \t 2:g_Ising \t 3:g_ACE \t 4:Emp_av \t 5:Model_av" << endl;
  for (it=list_I.begin(); it!=list_I.end(); it++)
  {
    cout << "Op = " << (*it).Op << " = " << bitset<n>((*it).Op) << "\t g_Ising = " << (*it).g_Ising << "\t\t" << "\t g_ACE = " << (*it).g_ACE;
    cout << "\t\t Emp_av = " << (*it).av_D << "\t Model_av = " << (*it).av_M << endl;
  }
  cout << endl;
}

/******************************************************************************/
/******************************      NOISE MODEL     **************************/
/******************************************************************************/
list<Interaction> Noise_Model(vector< pair<uint32_t,double> > IndepModel)
{
  Interaction I;
  list<Interaction> list_I;

  vector< pair<uint32_t,double> >::iterator BasisOp;

  for (BasisOp = IndepModel.begin(); BasisOp != IndepModel.end(); BasisOp++)
  {
    I.Op = (*BasisOp).first;
    I.k = 1;
    I.g_Ising = 0;
    I.g = I.g_Ising;
    I.g_ACE = -2.*I.g_Ising;
    I.av_D = 0; I.av_M = 0;
    list_I.push_back(I); 
  }

  return list_I;
}

/******************************************************************************/
/*********************     FULLY CONNECTED PAIRWISE    ************************/
/******************************************************************************/
list<Interaction> FullyConnectedPairwise()
{
  Interaction I;
  list<Interaction> list_I;

  uint32_t Op1 = 1, Op2 = 1;

  for(int i=1; i<=n; i++) // n fields
  {
    I.Op = Op1;    I.k = 1;
    I.g = 0;  I.g_Ising = 0;  I.g_ACE = 0;
    I.av_D = 0; I.av_M = 0;
    list_I.push_back(I);

    Op1 = Op1 << 1;
  }

  Op1 = 1;

  for(int i=1; i<=n; i++) // n(n-1)/2 pairwise interactions
  {    
    Op2 = Op1 << 1;
    for (int j=i+1; j<=n; j++)
    {
      I.Op = Op1 + Op2;     I.k = 2;
      I.g = 0;  I.g_Ising = 0;  I.g_ACE = 0;
      I.av_D = 0; I.av_M = 0;
      list_I.push_back(I);
      
      Op2 = Op2 << 1; 
    }
    Op1 = Op1 << 1;      
  }
  return list_I;
}

/******************************************************************************/
/***************************     PRINT PARAMETERS    **************************/
/******************************************************************************/
void find_ij(uint32_t Op_pair, unsigned int *i, unsigned int *j)
{
  unsigned int k = 0; uint32_t Op_test = 1;
  *i = n; *j = n;

  while( (k < n) && ((Op_test & Op_pair)!=Op_test) )
  {    Op_test = Op_test << 1; k++;   }
  *i = k;
  Op_test = Op_test << 1; k++;

  while( (k < n) && ((Op_test & Op_pair)!=Op_test) )
  {  Op_test = Op_test << 1; k++;     } 
  *j = k;
}

void Print_matrix_J_pairwise(string F_name, list<Interaction> list_I)
{
  list<Interaction>::iterator it;

  double **J = (double **)malloc(n*sizeof(double*));

  unsigned int i=0, j=0;
  for (i = 0; i < n; i++) 
  {  
    J[i] = (double*)malloc(n*sizeof(double));  
    for (j = 0; j < n; j++) 
      { J[i][j] = 0;  }
  }

  for (it = list_I.begin(); it != list_I.end(); it++)
  {
    find_ij((*it).Op, &i, &j);
    if (i==n) {   cout << "Error in \'Print_matrix_J_pairwise\'." << endl;  }
    else if (j==n) 
      {   
        J[(n-1-i)][(n-1-i)] = (*it).g; //0 
        cout << "field: " << (n-1-i) << " h = " << J[(n-1-i)][(n-1-i)] << endl;  
      }
    else {   J[(n-1-i)][(n-1-j)] = (*it).g;  J[(n-1-j)][(n-1-i)] = (*it).g;  cout << "pairwise: " << (n-1-i) << " " << (n-1-j) << " J = " << (*it).g << endl;  }
  }

  fstream fichier(F_name.c_str(), ios::out);
  for (i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++)
    {
      fichier << J[i][j] << " ";
    }
    fichier << endl;
  }

//  for (i = n-1; i >= 0; i--)
//  {
//    for(j = n-1; j >= 0; j--)
//    {   fichier << J[i][j] << " "; }
//    fichier << endl;
//  }

  fichier.close();
}

