#include <map>
#include <fstream>
#include <iostream>
#include <list>
#include "support.h"

using namespace std;

/******************************************************************************/
/**************** Log-likelihood (LogL), Geometric Complexity *****************/
/*************************  and Log-evidence (LogE) ***************************/
/******************************************************************************/
double LogL_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false);
double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, bool print_bool = false);
double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, unsigned int n, double *C_param, double *C_geom);

double LogE_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false);
double LogL_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, unsigned int n, bool print_bool = false);
double GeomComplexity_SubCM(unsigned int m);
double ParamComplexity_SubCM(unsigned int m, unsigned int N);


/******************************************************************************/
/***************************    Define an MCM   *******************************/
/******************************************************************************/
map<uint32_t, uint32_t> Create_MCM(list<uint32_t> MCM_table)
{
  map<uint32_t, uint32_t> MCM_partition;
  uint32_t integer = 0;

  for (auto const& i : MCM_table) {
    MCM_partition[integer]=i;
    integer++;
  }

  return MCM_partition;
}

/******************************************************************************/
/*** VERSION a) Operators are written as the binary          ******************/
/****           representation of the interactions           ******************/
/******************************************************************************/
map<uint32_t, uint32_t> Read_MCMParts_BinaryRepresentation(string MCM_binary_filename, unsigned int n)
{
  map<uint32_t, uint32_t> MCM_partition;
  uint32_t integer = 0;

  ifstream myfile (MCM_binary_filename.c_str());
  string line, line2;     

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line
      MCM_partition[integer]=stoi(line2, 0, 2);   //convert string line2 into a binary integer
      integer++;
    }
    myfile.close();
  }
  return MCM_partition;
}

/********************************************************************/
/*************    CHECK if "Partition" IS A PARTITION   *************/
/********************************************************************/
//check if *Partition* is an actual partition of r basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.

bool check_partition(map<uint32_t, uint32_t> Partition, unsigned int n)
{
  map<uint32_t, uint32_t>::iterator Part;
  uint32_t sum = 0;
  uint32_t rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += countSetBits((*Part).second);
    //cout << bitset<n>( (*Part).second ) << " \t";
  }
  //cout << bitset<n>(sum) << endl;

  return (countSetBits(sum) == rank);
}

/********************************************************************/
/*******    PRINT INFO on each PART of an MCM (= a partition)   *****/
/********************************************************************/
void PrintTerminal_MCM_Info(map<uint32_t, unsigned int > Kset, unsigned int N, unsigned int n, map<uint32_t, uint32_t> MCM_Partition)
{
  uint32_t Part = 0, m=0;
  double C_param=0, C_geom=0;
  Complexity_MCM(MCM_Partition, N, n, &C_param, &C_geom);
  double LogL = LogL_MCM(Kset, MCM_Partition, N, n);

  cout << "********** General Information about the MCM: **********" << endl; 
  cout << "Best MCM has " << MCM_Partition.size() << " partitions and the following properties:" << endl;
  cout << "\t LogL = " << LogL << endl;
  cout << " \t C_param = " << C_param << " \t \t C_geom = " << C_geom << endl;
  cout << " \t Total complexity = " << C_param + C_geom << endl;
  cout << " \t MDL = " << LogL - C_param - C_geom << endl;
  cout << "  \t LogE = " << LogE_MCM(Kset, MCM_Partition, N, n) << endl;

  cout << endl << "********** Information about each part of the MCM: **********";
  cout << endl << "\t (the total LogE of the model is the sum of the values for each part)";
  cout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  cout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  cout << "## 1:Part_int \t 2:Part_binary \t 3:LogL \t 4:C_param \t 5:C_geom \t 6:C_tot \t 7:LogE" << endl;

  for (map<uint32_t, uint32_t>::iterator i = MCM_Partition.begin(); i != MCM_Partition.end(); i++)
  {    
    Part = (*i).second;
    m = countSetBits(Part);  // rank of the part (i.e. rank of the SCM)
    C_param = ParamComplexity_SubCM(m, N);
    C_geom = GeomComplexity_SubCM(m);

    cout << " \t " << Part << " \t " << int_to_bstring(Part, n) << " \t";
    cout << LogL_SubCM(Kset, Part, N, n) << " \t";
    cout << C_param << " \t " << C_geom << " \t" << C_param + C_geom << " \t";
    cout << LogE_SubCM(Kset, Part, N, n) << endl;
  }
  cout << endl;
}

/********************************************************************/
/**************************    PRINT INFO    ************************/
/******    ON SUCCESSIVE INDEPENDENT MODELS IN THE NEW BASIS   ******/
/********************************************************************/
void PrintInfo_All_Indep_Models(map<uint32_t, unsigned int> Kset, unsigned int N, unsigned int n)
{
  map<uint32_t, uint32_t> Partition_Indep;  uint32_t Op = 1;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_Indep[i] = Op;
    cout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N, n) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N, n) << endl;    
    Op = Op << 1;
  }
  Partition_Indep.clear();
}

/********************************************************************/
/**************************    PRINT INFO    ************************/
/******    ON SUCCESSIVE SUB_COMPLETE MODELS IN THE NEW BASIS   *****/
/********************************************************************/
void PrintInfo_All_SubComplete_Models(map<uint32_t, unsigned int> Kset, unsigned int N, unsigned int n)
{
  map<uint32_t, uint32_t> Partition_SC;  uint32_t Op = 1;
  Partition_SC[0] = 0;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_SC[0] += Op;
    cout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_SC, N, n) << " \t LogL = " << LogL_MCM(Kset, Partition_SC, N, n) << endl;    
    Op = Op << 1;
  }
  Partition_SC.clear();
}

